#!/usr/bin/env python3
"""
kraken_build.py — download PHI-base eukaryotic pathogen transcriptomes and
build a Kraken2 database for co-infection screening of RNA-seq data.

Scope: fungi + oomycetes from phibase_db.json.
Downloads CDS FASTAs (cds_from_genomic.fna) for all annotated assemblies.
NCBI generates this file for any assembly with annotation (NCBI or author-provided).
Unannotated assemblies are skipped — all have zero detections in the SRA data.
Uses specific accessions from kraken/ref_screen.tsv — run screen_refs.py first.

Run from crypt/ on Setonix (requires kraken2 and datasets CLI in PATH):
    python kraken/build.py [--db-dir /scratch/kraken_db] [--genomes-dir /scratch/genomes]

Steps:
    0. Load reference map from kraken/ref_screen.tsv (best accession per seed)
    1. Load seed taxids from phibase_db.json (fungal_to_seed + oomycete_to_seed)
    2. For each seed: download FASTA via `datasets download genome accession`
         --include rna    (transcriptome) when assembly has annotation
         --include genome (genomic)       when no annotation available
    3. Tag FASTA headers with kraken:taxid|TAXID| (bypasses 15 GB accession2taxid)
    4. Add each FASTA to the Kraken2 library with `kraken2-build --add-to-library`
    5. Download taxdump + build the database with `kraken2-build --build`

Prerequisites:
    - kraken/ref_screen.tsv  (run: python kraken/screen_refs.py)
    - NCBI datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
    - Kraken2: https://github.com/DerrickWood/kraken2
    - ~2 GB disk for CDS FASTAs; ~4 GB for built Kraken2 DB (estimate)
    - Setonix: request high-memory node (512 GB RAM) for --build step

Output:
    {db_dir}/           Kraken2 database directory
    {genomes_dir}/      downloaded FASTAs (one subdir per taxid)
    output/kraken_build/logs/  build log
"""

import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, make_log_dir, link_latest

DB_PATH      = Path("output/00_build/data/phibase_db.json")
REF_SCREEN   = Path("kraken/ref_screen.tsv")
OUT_DIR      = Path("output/kraken_build")
DEFAULT_DB   = Path("/scratch/leon/kraken_db")
DEFAULT_GEN  = Path("/scratch/leon/kraken_transcriptomes")

TAXDUMP_URL   = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
ACACIA_BUCKET = "pawsey1168-llenzo-kraken-db"
ACACIA_ENDPOINT = "https://projects.pawsey.org.au"
ACACIA_PROFILE  = "acacia"

# Confidence threshold used at classify time (stored here for reference)
KRAKEN2_CONFIDENCE = 0.1


def load_ref_screen(tsv_path: Path) -> dict[int, dict]:
    """Load ref_screen.tsv → {taxid: {accession, fasta_type, has_annotation, notes}}."""
    ref = {}
    if not tsv_path.exists():
        print(f"WARNING: {tsv_path} not found — run kraken/screen_refs.py first",
              flush=True)
        return ref
    with open(tsv_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(header, parts))
            taxid = int(row["taxid"])
            ref[taxid] = {
                "accession":      row.get("best_accession", ""),
                "fasta_type":     row.get("fasta_type", "genome"),
                "has_annotation": row.get("has_annotation", "False") == "True",
                "release_date":   row.get("release_date", ""),
                "genes":          row.get("protein_coding_genes", ""),
                "notes":          row.get("notes", ""),
            }
    return ref


def load_seed_taxids(db_path: Path) -> tuple[dict[int, str], dict[str, str]]:
    """Return {seed_taxid: kingdom} for fungi + oomycetes."""
    with open(db_path) as f:
        raw = json.load(f)

    seeds: dict[int, str] = {}
    for taxid in set(raw["fungal_to_seed"].values()):
        seeds[int(taxid)] = "fungi"
    for taxid in set(raw["oomycete_to_seed"].values()):
        seeds[int(taxid)] = "oomycete"

    t2n = raw["taxid_to_name"]
    print(f"Seed taxids: {len(seeds)} "
          f"({sum(1 for k in seeds.values() if k=='fungi')} fungi, "
          f"{sum(1 for k in seeds.values() if k=='oomycete')} oomycetes)",
          flush=True)
    return seeds, t2n


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
    pattern = re.compile(r"^(>.+)")
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
    fasta_type: 'rna' (transcriptome) or 'genome' (genomic)
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
        # Some transcriptome downloads use .fa or .fna.gz
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
        ["kraken2-build", "--add-to-library", str(fna_path), "--db", str(db_dir)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"    WARNING: --add-to-library failed for {fna_path.name}: "
              f"{result.stderr.strip()[:200]}", flush=True)
        return False
    return True


def upload_to_acacia(local_dir: Path, s3_prefix: str) -> None:
    """
    Sync a local directory to Acacia S3 using aws s3 sync.
    s3_prefix: path within ACACIA_BUCKET, e.g. 'kraken_transcriptomes'
    """
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
    ap.add_argument("--db-path", default=str(DB_PATH),
                    help="Path to phibase_db.json")
    ap.add_argument("--ref-screen", default=str(REF_SCREEN),
                    help="Path to ref_screen.tsv (from screen_references.py)")
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
        seeds, t2n = load_seed_taxids(Path(args.db_path))

        n_rna      = 0
        n_skipped  = 0
        n_failed   = 0
        n_added    = 0
        seen_accessions: dict[str, int] = {}  # accession → taxid used, for dedup

        for i, (taxid, kingdom) in enumerate(sorted(seeds.items()), 1):
            name = t2n.get(str(taxid), str(taxid))
            print(f"\n[{i}/{len(seeds)}] taxid={taxid} ({kingdom}) {name}", flush=True)

            info = ref.get(taxid, {})
            accession = info.get("accession", "")
            fasta_type = info.get("fasta_type", "genome")

            if not accession:
                print(f"  SKIP: no assembly in ref_screen.tsv", flush=True)
                n_skipped += 1
                continue

            if fasta_type == "genome":
                # Unannotated assembly — no CDS FASTA available.
                # All unannotated seeds have zero detections in the SRA data.
                print(f"  SKIP: no annotation — zero detections in data", flush=True)
                n_skipped += 1
                continue

            if accession in seen_accessions:
                # PHI-base duplicate seeds sharing the same assembly (e.g. 105487/694573)
                print(f"  SKIP: duplicate of {accession} "
                      f"(already added as taxid={seen_accessions[accession]})", flush=True)
                n_skipped += 1
                continue

            if info.get("notes"):
                print(f"  NOTE: {info['notes']}", flush=True)

            taxon_dir = genomes_dir / f"{taxid}_{kingdom}"

            if not args.build_only:
                fnas = download_fasta(taxid, accession, fasta_type, taxon_dir, name)
                if fnas:
                    # Use taxid from assembly metadata — corrects stale PHI-base taxids
                    # (e.g. 105487/694573 → 578113 Cytospora mali)
                    tag_taxid = read_assembly_taxid(taxon_dir, taxid)
                    if tag_taxid != taxid:
                        print(f"  NOTE: taxid corrected {taxid} → {tag_taxid} "
                              f"(stale PHI-base taxid)", flush=True)
                    for fna in fnas:
                        tag_fasta_headers(fna, tag_taxid)
                    seen_accessions[accession] = tag_taxid
                    n_rna += 1
            else:
                fnas = list(taxon_dir.glob("**/*.fna"))
                if fnas:
                    seen_accessions[accession] = taxid

            if not fnas:
                n_failed += 1
                continue

            if not args.download_only:
                for fna in fnas:
                    if add_to_library(fna, db_dir):
                        n_added += 1

        print(f"\n── Download summary ──")
        print(f"  Seeds total:                     {len(seeds)}")
        print(f"  CDS FASTA (--include cds):       {n_rna}")
        print(f"  Skipped (no annotation):         {n_skipped}")
        print(f"  Download failed:                 {n_failed}")
        print(f"  FASTAs added to library:         {n_added}")

        if args.upload_to_acacia:
            upload_to_acacia(genomes_dir, "kraken_transcriptomes")

        if not args.download_only:
            download_taxonomy(db_dir)
            build_database(db_dir, threads=args.threads)
            print(f"\nKraken2 DB ready: {db_dir}")
            print(f"Classify with:  kraken2 --db {db_dir} "
                  f"--confidence {KRAKEN2_CONFIDENCE} --threads N ...")
        else:
            print("\n--download-only: skipped kraken2-build steps.")
            if args.upload_to_acacia:
                print(f"To retrieve on Setonix:")
                print(f"  aws s3 sync s3://{ACACIA_BUCKET}/kraken_transcriptomes/ "
                      f"/scratch/leon/kraken_transcriptomes/ "
                      f"--profile acacia --endpoint-url {ACACIA_ENDPOINT}")

    finally:
        log.close()


if __name__ == "__main__":
    main()
