#!/usr/bin/env python3
"""
kraken_db_build.py — build the Kraken2 pathogen DB from BUSCO-selected assemblies.
Submodule 1, step 3 of 3 (search → busco → build).

Reads kraken_db_busco/data/busco_scores.tsv, takes every row with selected=True,
and reads its CDS FASTA from kraken_db_search/data/cds/pathogen/{accession}/ — this
script never downloads anything itself; kraken_db_search.py is the one place in
the pipeline that fetches assemblies.

Scope: fungi + oomycetes selected by kraken_db_search.py + kraken_db_busco.py.
  - Seed assemblies: geographically/temporally diverse assemblies per PHI-base
    pathogen species that pass BUSCO (pan-genome coverage).
  - Genus fill-in: best annotated scaffold-plus assembly for each non-PHI-base
    species within detected PHI-base genera.
  - No BBDuk masking — Kraken2's LCA algorithm handles shared k-mers natively.
    Shared k-mers between species resolve to their LCA (genus/family), not to a
    false species-level hit. This allows a broad DB without sensitivity loss.

Host sequences are intentionally excluded — host identity is inferred from SRA
run metadata; including host CDS creates irreducible noise from k-mer similarity
between related plant genomes.

Run from crypt/ (requires kraken2 and datasets CLI in PATH):
    python kraken/db/kraken_db_build.py                              # default dirs
    python kraken/db/kraken_db_build.py --db-dir PATH --genomes-dir PATH   # Setonix scratch

Steps:
    1. Load selected assemblies from kraken_db_busco/data/busco_scores.tsv
       (run kraken_db_search.py + kraken_db_busco.py first)
    2. Read pre-downloaded CDS FASTA from --genomes-dir (no fetching here)
    3. Tag FASTA headers with kraken:taxid|TAXID| (bypasses 15 GB accession2taxid)
    4. Add all FASTAs to Kraken2 library (--no-masking keeps low-complexity regions)
    5. Download taxdump + build the database with `kraken2-build --build`

Prerequisites:
    - kraken_db_busco/data/busco_scores.tsv (run kraken_db_search.py + kraken_db_busco.py)
    - NCBI datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
    - Kraken2: https://github.com/DerrickWood/kraken2
    - Setonix: high-memory node (>= 64 GB RAM) for --build step

Output:
    {db_dir}/                   Kraken2 database directory
    kraken/output/db/build/logs/  build log
"""

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _util import _Tee, make_log_dir, link_latest, upload_to_acacia

BUSCO_SCORES = Path("kraken/output/db/busco/data/busco_scores.tsv")
OUT_DIR      = Path("kraken/output/db/build")
DEFAULT_DB          = OUT_DIR / "data" / "db"
DEFAULT_GENOMES_DIR = Path("kraken/output/db/search/data/cds/pathogen")

TAXDUMP_URL   = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
ACACIA_BUCKET = "pawsey1168-llenzo-kraken-db"

# Classification parameters (used at classify time; stored here for reference)
KRAKEN2_CONFIDENCE     = 0.15
KRAKEN2_MIN_HIT_GROUPS = 3


def load_selected(tsv_path: Path) -> list:
    """Load busco_scores.tsv → list of selected assembly dicts (selected=True only)."""
    if not tsv_path.exists():
        print(f"WARNING: {tsv_path} not found — run kraken_db_search.py + "
              f"kraken_db_busco.py first", flush=True)
        return []
    rows = []
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("selected") != "True" or not row.get("accession"):
                continue
            rows.append({
                "taxid":         int(row["taxid"]),
                "organism_name": row.get("organism_name", str(row["taxid"])),
                "kingdom":       row.get("kingdom", ""),
                "source":        row.get("source", "seed"),
                "accession":     row["accession"],
                "fasta_type":    row.get("fasta_type", "genome"),
                "busco_status":  row.get("busco_status", ""),
            })
    n_seed = sum(1 for r in rows if r["source"] == "seed")
    n_fill = sum(1 for r in rows if r["source"] == "genus_fill")
    n_fallback = sum(1 for r in rows if "fallback" in r["busco_status"])
    print(f"Selected assemblies: {len(rows)} ({n_seed} seed, {n_fill} genus fill-in, "
          f"{n_fallback} below-threshold fallback)", flush=True)
    return rows


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
    """Rewrite FASTA headers to kraken:taxid|TAXID|original_header format.
    Kraken2 reads this directly, bypassing the 15 GB accession2taxid lookup.
    Returns the path of the tagged file (replaces original in-place). No-op if
    already tagged (idempotent across re-builds without re-downloading)."""
    with open(fasta_path, "rb") as f:
        head = f.read(30)
    if head.startswith(b">kraken:taxid|"):
        return fasta_path
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


def add_to_library(fna_path: Path, db_dir: Path) -> bool:
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


def download_taxonomy(db_dir: Path) -> None:
    taxdump_path = db_dir / "taxdump.tar.gz"
    if (db_dir / "taxonomy" / "nodes.dmp").exists():
        print("Taxonomy already installed.", flush=True)
        return
    print(f"Downloading taxonomy from {TAXDUMP_URL} …", flush=True)
    (db_dir / "taxonomy").mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(TAXDUMP_URL, taxdump_path)
    subprocess.run(["tar", "-xzf", str(taxdump_path), "-C", str(db_dir / "taxonomy")], check=True)
    taxdump_path.unlink(missing_ok=True)
    print("Taxonomy installed.", flush=True)


def build_database(db_dir: Path, threads: int = 32) -> None:
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
    ap.add_argument("--busco-scores", default=str(BUSCO_SCORES),
                    help="Path to busco_scores.tsv (from kraken_db_busco.py)")
    ap.add_argument("--db-dir", default=str(DEFAULT_DB),
                    help="Kraken2 database output directory")
    ap.add_argument("--genomes-dir", default=str(DEFAULT_GENOMES_DIR),
                    help="Directory with CDS FASTA already downloaded by kraken_db_search.py")
    ap.add_argument("--threads", type=int, default=32,
                    help="Threads for kraken2-build --build (default: 32)")
    ap.add_argument("--build-only", action="store_true",
                    help="Skip library add; taxonomy download + build only "
                         "(assumes library already populated)")
    ap.add_argument("--upload-to-acacia", action="store_true",
                    help="Sync --genomes-dir to Acacia "
                         f"(s3://{ACACIA_BUCKET}/kraken_transcriptomes/) before building")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_db_build.log")
    link_latest(logs_base, log_dir / "kraken_db_build.log")
    sys.stdout = log

    try:
        db_dir      = Path(args.db_dir)
        genomes_dir = Path(args.genomes_dir)
        db_dir.mkdir(parents=True, exist_ok=True)

        selected = load_selected(Path(args.busco_scores))
        if not selected:
            print("No assemblies to build with. Run kraken_db_search.py + "
                  "kraken_db_busco.py first.", flush=True)
            return

        if args.upload_to_acacia:
            upload_to_acacia(genomes_dir, "kraken_transcriptomes", bucket=ACACIA_BUCKET)

        n_found, n_missing, n_skipped, n_added = 0, 0, 0, 0
        seen_accessions: set = set()
        prepared: list = []  # (tag_taxid, fnas)

        if not args.build_only:
            print(f"\n── Preparing FASTAs from {genomes_dir} ──────────────────────────")
            for i, asm in enumerate(selected, 1):
                taxid, name = asm["taxid"], asm["organism_name"]
                accession, fasta_type, source = asm["accession"], asm["fasta_type"], asm["source"]
                print(f"\n[{i}/{len(selected)}] {accession}  taxid={taxid}  ({source}) {name}",
                      flush=True)

                if fasta_type != "cds":
                    print("  SKIP: no annotation (CDS-only DB)", flush=True)
                    n_skipped += 1
                    continue
                if accession in seen_accessions:
                    print("  SKIP: duplicate accession", flush=True)
                    n_skipped += 1
                    continue

                taxon_dir = genomes_dir / accession
                fnas = [f for f in taxon_dir.glob("**/*.fna")
                        if not f.name.endswith(("_clean.fna", "_combined.fna"))]
                if not fnas:
                    print(f"  MISSING: no CDS at {taxon_dir} — run "
                          f"kraken_db_search.py --download first", flush=True)
                    n_missing += 1
                    continue

                tag_taxid = read_assembly_taxid(taxon_dir, taxid)
                if tag_taxid != taxid:
                    print(f"  NOTE: taxid corrected {taxid} → {tag_taxid} "
                          f"(stale taxid in ref_candidates)", flush=True)
                for fna in fnas:
                    tag_fasta_headers(fna, tag_taxid)
                seen_accessions.add(accession)
                n_found += 1
                prepared.append((tag_taxid, fnas))

            print(f"\n── Prepare summary ──────────────────────────────────────────────")
            print(f"  Selected total: {len(selected)}")
            print(f"  Found on disk:  {n_found}")
            print(f"  Skipped:        {n_skipped}")
            print(f"  Missing:        {n_missing}")

            # ── Add to library — no BBDuk masking, Kraken2 LCA handles shared k-mers ──
            print(f"\n── Adding to Kraken2 library ────────────────────────────────────")
            lib_dir = db_dir / "library"
            if lib_dir.exists() and list(lib_dir.iterdir()):
                print(f"  Clearing existing library ({lib_dir}) …", flush=True)
                shutil.rmtree(lib_dir)
            for tag_taxid, fnas in prepared:
                for fna in fnas:
                    if add_to_library(fna, db_dir):
                        n_added += 1
            print(f"  FASTAs added to library: {n_added}")

        download_taxonomy(db_dir)
        build_database(db_dir, threads=args.threads)

        print(f"\nKraken2 DB ready: {db_dir}")
        print(f"Classify with:  kraken2 --db {db_dir} "
              f"--confidence {KRAKEN2_CONFIDENCE} "
              f"--minimum-hit-groups {KRAKEN2_MIN_HIT_GROUPS} --threads N ...")

    finally:
        log.close()


if __name__ == "__main__":
    main()
