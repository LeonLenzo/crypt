#!/usr/bin/env python3
"""
build.py — download selected assemblies and build the Kraken2 pathogen DB.

Scope: fungi + oomycetes selected by kraken/screen_refs.py.
  - Seed assemblies: up to 10 geographically/temporally diverse assemblies per
    PHI-base pathogen species (pan-genome coverage).
  - Genus fill-in: best annotated scaffold-plus assembly for each non-PHI-base
    species within detected PHI-base genera.
  - No BBDuk masking — Kraken2's LCA algorithm handles shared k-mers natively.
    Shared k-mers between species resolve to their LCA (genus/family), not to a
    false species-level hit. This allows a broad DB without sensitivity loss.

Host sequences are intentionally excluded — host identity is inferred from SRA
run metadata; including host CDS creates irreducible noise from k-mer similarity
between related plant genomes.

Run from crypt/ (requires kraken2 and datasets CLI in PATH):
    python kraken/build.py                         # default dirs
    python kraken/build.py [--db-dir PATH] [--genomes-dir PATH]   # Setonix scratch

Steps:
    1. Load assembly list from kraken/ref_screen.tsv (run screen_refs.py first)
    2. For each assembly: download FASTA via `datasets download genome accession`
         --include cds  (annotated assemblies)
         --include genome  (unannotated fallback)
    3. Tag FASTA headers with kraken:taxid|TAXID| (bypasses 15 GB accession2taxid)
    4. Add all FASTAs to Kraken2 library (--no-masking keeps low-complexity regions)
    5. Download taxdump + build the database with `kraken2-build --build`

Prerequisites:
    - kraken/ref_screen.tsv  (run: python kraken/screen_refs.py)
    - NCBI datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
    - Kraken2: https://github.com/DerrickWood/kraken2
    - ~5 GB disk for FASTAs; ~8 GB for built Kraken2 DB
    - Setonix: high-memory node (≥ 64 GB RAM) for --build step

Output:
    {db_dir}/           Kraken2 database directory
    {genomes_dir}/      downloaded FASTAs (one subdir per accession)
    kraken/output/build/logs/  build log
"""

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, make_log_dir, link_latest

REF_SCREEN   = Path("kraken/ref_screen.tsv")
OUT_DIR      = Path("kraken/output/build")
DEFAULT_DB   = Path("kraken/output/kraken_db_build/data/db")
DEFAULT_GEN  = Path("kraken/output/kraken_db_search/data/cds_from_genomic")

TAXDUMP_URL     = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
ACACIA_BUCKET   = "pawsey1168-llenzo-kraken-db"
ACACIA_ENDPOINT = "https://projects.pawsey.org.au"
ACACIA_PROFILE  = "acacia"

# Classification parameters (used at classify time; stored here for reference)
KRAKEN2_CONFIDENCE     = 0.15
KRAKEN2_MIN_HIT_GROUPS = 3


def load_ref_screen(tsv_path: Path) -> list[dict]:
    """Load ref_screen.tsv → list of assembly dicts, one row per selected assembly."""
    if not tsv_path.exists():
        print(f"WARNING: {tsv_path} not found — run kraken/screen_refs.py first",
              flush=True)
        return []
    rows = []
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("accession"):
                continue
            rows.append({
                "taxid":         int(row["taxid"]),
                "organism_name": row.get("organism_name", str(row["taxid"])),
                "kingdom":       row.get("kingdom", ""),
                "source":        row.get("source", "seed"),
                "accession":     row["accession"],
                "fasta_type":    row.get("fasta_type", "genome"),
                "has_annotation": row.get("has_annotation", "False") == "True",
            })
    n_seed = sum(1 for r in rows if r["source"] == "seed")
    n_fill = sum(1 for r in rows if r["source"] == "genus_fill")
    print(f"Assemblies in ref_screen.tsv: {len(rows)} "
          f"({n_seed} seed pan-genome, {n_fill} genus fill-in)", flush=True)
    return rows


def scan_existing_accessions(genomes_dir: Path) -> dict[str, Path]:
    """Scan genomes_dir for already-downloaded assemblies regardless of directory naming.

    Works with both old-style ({taxid}_{kingdom}/) and new-style ({accession}/)
    subdirectory layouts by reading assembly_data_report.jsonl files bundled with
    every NCBI datasets download.

    Returns {accession: dir_containing_fna_files}.
    """
    found: dict[str, Path] = {}
    for report in genomes_dir.glob("**/assembly_data_report.jsonl"):
        try:
            data = json.loads(report.read_text())
            acc = data.get("accession", "") or report.parent.name
            if not (acc.startswith("GCA_") or acc.startswith("GCF_")):
                continue
            fnas = list(report.parent.glob("**/*.fna")) or list(report.parent.parent.glob("*.fna"))
            if fnas:
                found[acc] = report.parent
        except Exception:
            pass
    return found


def read_assembly_taxid(taxon_dir: Path, seed_taxid: int) -> int:
    """Read organism taxid from assembly_data_report.jsonl shipped with the download.
    Falls back to seed_taxid if the report is absent or malformed — but the three
    PHI-base seeds with stale taxids (105487, 694573, 914237) need this to get the
    correct current NCBI taxid embedded in Kraken2 headers."""
    reports = list(taxon_dir.glob("**/assembly_data_report.jsonl"))
    if not reports:
        return seed_taxid
    try:
        d = json.loads(reports[0].read_text())
        taxid = d.get("organism", {}).get("taxId")
        if taxid:
            return int(taxid)
    except Exception:
        pass
    return seed_taxid


def tag_fasta_headers(fasta_path: Path, taxid: int) -> Path:
    """
    Rewrite FASTA headers to kraken:taxid|TAXID|original_header format.
    Kraken2 reads this directly, bypassing the 15 GB accession2taxid lookup.
    Returns the path of the tagged file (replaces original in-place).
    """
    tagged_path = fasta_path.with_suffix(".tagged.fna")
    with open(fasta_path) as fin, open(tagged_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip()
                fout.write(f">kraken:taxid|{taxid}|{header}\n")
            else:
                fout.write(line)
    fasta_path.unlink()
    tagged_path.rename(fasta_path)
    return fasta_path


def download_fasta(taxid: int, accession: str, fasta_type: str,
                   dest_dir: Path, name: str) -> list[Path]:
    """
    Download FASTA for a specific accession via NCBI datasets CLI.
    fasta_type: 'cds' (coding sequences) or 'genome' (genomic)
    Returns list of downloaded .fna paths (may be empty on failure).
    Skips download if dest_dir already contains tagged .fna files.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)

    existing = list(dest_dir.glob("*.fna")) + list(dest_dir.glob("**/*.fna"))
    if existing:
        print(f"  [{taxid}] {name[:50]}: already downloaded ({len(existing)} fna)",
              flush=True)
        return existing

    include_flag = "cds" if fasta_type == "cds" else "genome"
    print(f"  [{taxid}] {name[:50]}: downloading {accession} "
          f"(--include {include_flag}) …", flush=True)

    zip_path = dest_dir / "ncbi_dataset.zip"
    cmd = [
        "datasets", "download", "genome", "accession", accession,
        "--include", include_flag,
        "--filename", str(zip_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not zip_path.exists():
        print(f"    WARNING: download failed for {accession} ({name}): "
              f"{result.stderr.strip()[:200]}", flush=True)
        return []

    unzip = subprocess.run(
        ["unzip", "-q", "-o", str(zip_path), "-d", str(dest_dir)],
        capture_output=True
    )
    zip_path.unlink(missing_ok=True)
    if unzip.returncode != 0:
        print(f"    WARNING: unzip failed for {accession}", flush=True)
        return []

    fnas = list(dest_dir.glob("**/*.fna"))
    if not fnas:
        for gz in dest_dir.glob("**/*.fna.gz"):
            out = gz.with_suffix("")
            with gzip.open(gz, "rb") as fin, open(out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            gz.unlink()
        fnas = list(dest_dir.glob("**/*.fna"))

    print(f"    → {len(fnas)} fna file(s) ({include_flag})", flush=True)
    return fnas


def add_to_library(fna_path: Path, db_dir: Path) -> bool:
    """Add a FASTA to the Kraken2 library."""
    result = subprocess.run(
        ["kraken2-build", "--add-to-library", str(fna_path),
         "--db", str(db_dir), "--no-masking"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"    WARNING: --add-to-library failed for {fna_path.name}: "
              f"{result.stderr.strip()[:200]}", flush=True)
        return False
    return True




def upload_to_acacia(local_dir: Path, s3_prefix: str) -> None:
    """Sync a local directory to Acacia S3 using aws s3 sync."""
    s3_uri = f"s3://{ACACIA_BUCKET}/{s3_prefix}/"
    print(f"\nUploading {local_dir} → {s3_uri} …", flush=True)
    cmd = [
        "aws", "s3", "sync", str(local_dir), s3_uri,
        "--profile", ACACIA_PROFILE,
        "--endpoint-url", ACACIA_ENDPOINT,
    ]
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        print(f"WARNING: upload to Acacia failed (exit {result.returncode})", flush=True)
    else:
        print(f"Upload complete → {s3_uri}", flush=True)


def download_taxonomy(db_dir: Path) -> None:
    """Download and install NCBI taxonomy into the Kraken2 DB directory."""
    taxdump_path = db_dir / "taxdump.tar.gz"
    if (db_dir / "taxonomy" / "nodes.dmp").exists():
        print("Taxonomy already installed.", flush=True)
        return
    print(f"Downloading taxonomy from {TAXDUMP_URL} …", flush=True)
    (db_dir / "taxonomy").mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(TAXDUMP_URL, taxdump_path)
    subprocess.run(
        ["tar", "-xzf", str(taxdump_path), "-C", str(db_dir / "taxonomy")],
        check=True
    )
    taxdump_path.unlink(missing_ok=True)
    print("Taxonomy installed.", flush=True)


def build_database(db_dir: Path, threads: int = 32) -> None:
    """Run kraken2-build --build."""
    print(f"\nBuilding Kraken2 database in {db_dir} (threads={threads}) …", flush=True)
    result = subprocess.run(
        ["kraken2-build", "--build", "--db", str(db_dir), "--threads", str(threads)],
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"kraken2-build --build failed (exit {result.returncode})")
    print("Database build complete.", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-screen", default=str(REF_SCREEN),
                    help="Path to ref_screen.tsv (from kraken/screen_refs.py)")
    ap.add_argument("--db-dir", default=str(DEFAULT_DB),
                    help="Kraken2 database output directory")
    ap.add_argument("--genomes-dir", default=str(DEFAULT_GEN),
                    help="Directory for downloaded FASTAs")
    ap.add_argument("--threads", type=int, default=32,
                    help="Threads for kraken2-build --build (default: 32)")
    ap.add_argument("--download-only", action="store_true",
                    help="Download FASTAs only; skip kraken2-build steps")
    ap.add_argument("--build-only", action="store_true",
                    help="Skip downloads; add existing FASTAs to library and build")
    ap.add_argument("--upload-to-acacia", action="store_true",
                    help="After downloading, sync FASTAs to Acacia "
                         f"(s3://{ACACIA_BUCKET}/kraken_transcriptomes/)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_build.log")
    link_latest(logs_base, log_dir / "kraken_build.log")
    sys.stdout = log

    try:
        db_dir      = Path(args.db_dir)
        genomes_dir = Path(args.genomes_dir)
        if not args.download_only:
            db_dir.mkdir(parents=True, exist_ok=True)
        genomes_dir.mkdir(parents=True, exist_ok=True)

        ref = load_ref_screen(Path(args.ref_screen))
        if not ref:
            print("No assemblies to process. Run kraken/screen_refs.py first.", flush=True)
            return

        n_downloaded = 0
        n_cached     = 0
        n_skipped    = 0
        n_failed     = 0
        n_added      = 0
        seen_accessions: set[str] = set()

        # ── Phase 1: Download ─────────────────────────────────────────────────
        downloaded: list[tuple[int, list[Path]]] = []  # (taxid, fnas)

        # Pre-scan: find assemblies already on disk (old-style or new-style dirs)
        if not args.build_only:
            existing = scan_existing_accessions(genomes_dir)
            if existing:
                print(f"Pre-scan: {len(existing)} accessions already on disk in {genomes_dir}",
                      flush=True)
        else:
            existing = {}

        for i, asm in enumerate(ref, 1):
            taxid      = asm["taxid"]
            name       = asm["organism_name"]
            accession  = asm["accession"]
            fasta_type = asm["fasta_type"]
            source     = asm["source"]

            print(f"\n[{i}/{len(ref)}] {accession}  taxid={taxid}  ({source}) {name}",
                  flush=True)

            if fasta_type == "genome":
                print("  SKIP: no annotation (CDS-only DB)", flush=True)
                n_skipped += 1
                continue

            if accession in seen_accessions:
                print("  SKIP: duplicate accession", flush=True)
                n_skipped += 1
                continue

            # One subdir per accession — supports multiple assemblies per taxid
            taxon_dir = genomes_dir / accession

            if not args.build_only and accession in existing:
                # Already on disk from a previous run (any directory layout)
                cached_dir = existing[accession]
                fnas = list(cached_dir.glob("**/*.fna"))
                tag_taxid = read_assembly_taxid(cached_dir, taxid)
                print(f"  CACHED: {accession} found in {cached_dir.parent.name}/ "
                      f"({len(fnas)} fna)", flush=True)
                seen_accessions.add(accession)
                n_cached += 1
            elif not args.build_only:
                fnas = download_fasta(taxid, accession, fasta_type, taxon_dir, name)
                if fnas:
                    tag_taxid = read_assembly_taxid(taxon_dir, taxid)
                    if tag_taxid != taxid:
                        print(f"  NOTE: taxid corrected {taxid} → {tag_taxid} "
                              f"(stale taxid in ref_screen)", flush=True)
                    for fna in fnas:
                        tag_fasta_headers(fna, tag_taxid)
                    seen_accessions.add(accession)
                    n_downloaded += 1
            else:
                fnas = list(taxon_dir.glob("**/*.fna"))
                tag_taxid = taxid
                if fnas:
                    seen_accessions.add(accession)

            if not fnas:
                n_failed += 1
                continue

            downloaded.append((tag_taxid, fnas))

        print(f"\n── Download summary ──────────────────────────────────────────")
        print(f"  Assemblies total:  {len(ref)}")
        print(f"  Cached (reused):   {n_cached}")
        print(f"  Downloaded (new):  {n_downloaded}")
        print(f"  Skipped:           {n_skipped}")
        print(f"  Failed:            {n_failed}")

        if args.upload_to_acacia:
            upload_to_acacia(genomes_dir, "kraken_transcriptomes")

        if args.download_only:
            print("\n--download-only: skipped kraken2-build steps.")
            if args.upload_to_acacia:
                print(f"To retrieve on Setonix:")
                print(f"  aws s3 sync s3://{ACACIA_BUCKET}/kraken_transcriptomes/ "
                      f"{DEFAULT_GEN}/ "
                      f"--profile acacia --endpoint-url {ACACIA_ENDPOINT}")
            return

        # ── Phase 2: Add to library ───────────────────────────────────────────
        # No BBDuk masking — Kraken2 LCA handles shared k-mers natively.
        print(f"\n── Adding to Kraken2 library ────────────────────────────────────")

        lib_dir = db_dir / "library"
        if lib_dir.exists() and list(lib_dir.iterdir()):
            print(f"  Clearing existing library ({lib_dir}) …", flush=True)
            shutil.rmtree(lib_dir)

        for tag_taxid, fnas in downloaded:
            for fna in fnas:
                if add_to_library(fna, db_dir):
                    n_added += 1

        print(f"\n── Library summary ──────────────────────────────────────────────")
        print(f"  FASTAs added to library: {n_added}")

        # ── Phase 3: Download taxonomy + build ────────────────────────────────
        download_taxonomy(db_dir)
        build_database(db_dir, threads=args.threads)

        print(f"\nKraken2 DB ready: {db_dir}")
        print(f"Classify with:  kraken2 --db {db_dir} "
              f"--confidence {KRAKEN2_CONFIDENCE} "
              f"--minimum-hit-groups {KRAKEN2_MIN_HIT_GROUPS} "
              f"--threads N ...")

    finally:
        log.close()


if __name__ == "__main__":
    main()
