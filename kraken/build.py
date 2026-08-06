#!/usr/bin/env python3
"""
kraken_build.py — download PHI-base eukaryotic pathogen transcriptomes and
build a Kraken2 database for co-infection screening of RNA-seq data.

Scope: fungi + oomycetes from phibase_db.json.
Downloads CDS FASTAs (cds_from_genomic.fna) for all annotated assemblies.
NCBI generates this file for any assembly with annotation (NCBI or author-provided).
Unannotated assemblies are skipped — all have zero detections in the SRA data.
Uses specific accessions from kraken/ref_screen.tsv — run screen_refs.py first.

Run from crypt/ (requires kraken2, datasets CLI, and bbduk.sh in PATH):
    python kraken/build.py                         # uses kraken/db + kraken/cds_from_genomic
    python kraken/build.py [--db-dir PATH] [--genomes-dir PATH]   # override for Setonix scratch

Steps:
    0. Load reference map from kraken/ref_screen.tsv (best accession per seed)
    1. Load seed taxids from phibase_db.json (fungal_to_seed + oomycete_to_seed)
    2. For each seed: download FASTA via `datasets download genome accession`
         --include cds (coding sequences) when assembly has annotation
         --include genome (genomic)       when no annotation available
    3. Tag FASTA headers with kraken:taxid|TAXID| (bypasses 15 GB accession2taxid)
    4. BBDuk masking: remove k-mers from host CDS that also appear in pathogen CDS.
         Prevents conserved eukaryotic k-mers (ribosomes, histones, metabolic genes)
         from matching non-host reads and producing spurious host assignments.
    5. Add each FASTA to the Kraken2 library with `kraken2-build --add-to-library`
         (host FASTAs added as a single BBDuk-masked combined file)
    6. Download taxdump + build the database with `kraken2-build --build`

Prerequisites:
    - kraken/ref_screen.tsv  (run: python kraken/screen_refs.py)
    - NCBI datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
    - Kraken2: https://github.com/DerrickWood/kraken2
    - BBDuk (BBTools ≥ 38.0): https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
    - ~2 GB disk for CDS FASTAs; ~10 GB for built Kraken2 DB
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
DEFAULT_DB   = Path("kraken/db")
DEFAULT_GEN  = Path("kraken/cds_from_genomic")

TAXDUMP_URL   = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
ACACIA_BUCKET = "pawsey1168-llenzo-kraken-db"
ACACIA_ENDPOINT = "https://projects.pawsey.org.au"
ACACIA_PROFILE  = "acacia"

# Classification parameters (used at classify time; stored here for reference)
KRAKEN2_CONFIDENCE     = 0.15
KRAKEN2_MIN_HIT_GROUPS = 3

# BBDuk k-mer length — must match Kraken2 kmer-len (default 35)
BBDUK_KMER_LEN = 35


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
    """Return {seed_taxid: kingdom} for fungi + oomycetes + hosts (if in ref_screen.tsv)."""
    with open(db_path) as f:
        raw = json.load(f)

    seeds: dict[int, str] = {}
    for taxid in set(raw["fungal_to_seed"].values()):
        seeds[int(taxid)] = "fungi"
    for taxid in set(raw["oomycete_to_seed"].values()):
        seeds[int(taxid)] = "oomycete"
    for taxid in set(raw["host_to_seed"].values()):
        seeds[int(taxid)] = "host"

    t2n = raw["taxid_to_name"]
    n_fungi    = sum(1 for k in seeds.values() if k == "fungi")
    n_oomycete = sum(1 for k in seeds.values() if k == "oomycete")
    n_host     = sum(1 for k in seeds.values() if k == "host")
    print(f"Seed taxids: {len(seeds)} "
          f"({n_fungi} fungi, {n_oomycete} oomycetes, {n_host} hosts)",
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
        ["kraken2-build", "--add-to-library", str(fna_path), "--db", str(db_dir)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"    WARNING: --add-to-library failed for {fna_path.name}: "
              f"{result.stderr.strip()[:200]}", flush=True)
        return False
    return True


def _concat_fnas(fnas: list[Path], dest: Path) -> None:
    """Concatenate FASTA files into dest."""
    with open(dest, "wb") as out:
        for p in fnas:
            with open(p, "rb") as f:
                shutil.copyfileobj(f, out)


def _kmer_shared(combined: Path, out: Path, threads: int) -> Path | None:
    """
    Extract k-mers appearing ≥2 times in combined.fna using kmercountexact.sh.
    These are k-mers shared between ≥2 sequences (cross-species) or duplicated
    within a single genome (repetitive CDS — minor over-masking, acceptable).
    Output is a FASTA of k-mer sequences usable as a BBDuk reference.
    Returns out path on success, None on failure.
    """
    if out.exists():
        print(f"    shared k-mers already computed ({out.name})", flush=True)
        return out

    sz_gb = combined.stat().st_size / 1e9
    print(f"    kmercountexact: counting k={BBDUK_KMER_LEN}-mers in "
          f"{combined.name} ({sz_gb:.2f} GB) …", flush=True)
    result = subprocess.run([
        "kmercountexact.sh",
        f"in={combined}",
        f"out={out}",
        f"k={BBDUK_KMER_LEN}",
        "fastadump=t",   # output as FASTA (BBDuk-readable reference)
        "mincount=2",    # only k-mers shared between ≥2 locations
        f"threads={threads}",
        "overwrite=t",
    ], capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        print(f"    kmercountexact ERROR (exit {result.returncode}):", flush=True)
        print(result.stderr[-1000:], flush=True)
        return None

    n_shared = sum(1 for l in result.stdout.splitlines() if l.startswith(">"))
    print(f"    → {n_shared:,} shared k-mers written to {out.name}", flush=True)
    return out


def _run_bbduk(in_fna: Path, out_fna: Path, ref_paths: list[Path], threads: int) -> bool:
    """Mask k-mers in in_fna that appear in any ref_paths file. kmask=N."""
    if out_fna.exists():
        print(f"    already masked: {out_fna.name}", flush=True)
        return True

    ref_str = ",".join(str(p) for p in ref_paths)
    result = subprocess.run([
        "bbduk.sh",
        f"in={in_fna}",
        f"out={out_fna}",
        f"ref={ref_str}",
        f"k={BBDUK_KMER_LEN}",
        "hdist=0",
        "kmask=N",
        f"threads={threads}",
        "overwrite=t",
    ], capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        print(f"    BBDuk ERROR (exit {result.returncode}) for {in_fna.name}:", flush=True)
        print(result.stderr[-500:], flush=True)
        return False

    for line in result.stderr.splitlines():
        s = line.strip()
        if s and any(x in s for x in ["Reads", "Bases", "Mask", "Result"]):
            print(f"    {s}", flush=True)
    return True


def bbduk_mask_sequences(genomes_dir: Path, threads: int = 8) -> tuple[Path | None, Path | None]:
    """
    Reduce all CDS sequences to species-diagnostic k-mers using a two-pass approach:

    Pass 1 (pathogen): identify k-mers shared between ≥2 pathogen sequences
      (kmercountexact.sh mincount=2), then mask pathogens against those shared
      k-mers + all host CDS. Eliminates same-order cross-hits (e.g. Melampsora
      in PST runs) and cross-kingdom hits (plant genes in fungal libraries).

    Pass 2 (host): identify k-mers shared between ≥2 host sequences, then mask
      hosts against those shared k-mers + all pathogen CDS. Eliminates same-family
      cross-hits (e.g. Sorghum in maize runs) and cross-kingdom hits.

    Kraken2's LCA algorithm already promotes k-mers shared between species IN THE
    DATABASE to ancestor nodes. This masking extends that principle to handle cases
    where one reference has sequences that match reads from a different species
    (due to reference contamination or divergent sequences not captured in both).

    Input : {genomes_dir}/{taxid}_{kingdom}/**/*.fna (already tagged)
    Output: (_pathogen_masked.fna, _host_masked.fna)  ready for kraken2-build
    Intermediates kept in genomes_dir for inspection and resumability.
    """
    # Collect FASTAs
    pathogen_fnas = sorted(
        fna
        for d in sorted(genomes_dir.iterdir())
        if d.is_dir() and not d.name.startswith("_") and not d.name.endswith("_host")
        for fna in d.glob("**/*.fna")
    )
    host_fnas = sorted(
        fna
        for d in sorted(genomes_dir.iterdir())
        if d.is_dir() and d.name.endswith("_host")
        for fna in d.glob("**/*.fna")
    )

    print(f"  BBDuk: {len(pathogen_fnas)} pathogen FASTAs "
          f"({sum(p.stat().st_size for p in pathogen_fnas)/1e9:.2f} GB), "
          f"{len(host_fnas)} host FASTAs "
          f"({sum(h.stat().st_size for h in host_fnas)/1e9:.2f} GB)", flush=True)

    if not pathogen_fnas and not host_fnas:
        print("  BBDuk: no FASTAs found — skipping", flush=True)
        return None, None

    # ── Concatenate combined inputs ────────────────────────────────────────────
    combined_pathogens = genomes_dir / "_pathogen_combined.fna"
    combined_hosts     = genomes_dir / "_host_combined.fna"

    if pathogen_fnas and not combined_pathogens.exists():
        print(f"  Concatenating pathogens → {combined_pathogens.name} …", flush=True)
        _concat_fnas(pathogen_fnas, combined_pathogens)

    if host_fnas and not combined_hosts.exists():
        print(f"  Concatenating hosts → {combined_hosts.name} …", flush=True)
        _concat_fnas(host_fnas, combined_hosts)

    # ── Pass 1: Mask pathogens (intra-pathogen shared + cross-kingdom) ─────────
    masked_pathogens = genomes_dir / "_pathogen_masked.fna"
    if not masked_pathogens.exists() and pathogen_fnas:
        print("\n  ── Pass 1: pathogen masking ──────────────────────────────────",
              flush=True)
        shared_pathogen_kmers = genomes_dir / "_shared_pathogen_kmers.fna"
        pathogen_ref_files: list[Path] = []

        # k-mers shared between ≥2 pathogen sequences (same-order cross-hits)
        if len(pathogen_fnas) > 1:
            r = _kmer_shared(combined_pathogens, shared_pathogen_kmers, threads)
            if r:
                pathogen_ref_files.append(r)

        # Also mask against all host CDS (cross-kingdom hits)
        if combined_hosts.exists():
            pathogen_ref_files.append(combined_hosts)

        if pathogen_ref_files:
            print(f"    BBDuk: masking pathogens "
                  f"(ref={[p.name for p in pathogen_ref_files]}) …", flush=True)
            _run_bbduk(combined_pathogens, masked_pathogens, pathogen_ref_files, threads)
        else:
            print("    No reference available; copying pathogens unmasked.", flush=True)
            shutil.copy2(str(combined_pathogens), str(masked_pathogens))
    elif masked_pathogens.exists():
        print(f"  _pathogen_masked.fna already exists "
              f"({masked_pathogens.stat().st_size/1e9:.2f} GB)", flush=True)

    # ── Pass 2: Mask hosts (intra-host shared + cross-kingdom) ───────────────
    masked_hosts = genomes_dir / "_host_masked.fna"
    if not masked_hosts.exists() and host_fnas:
        print("\n  ── Pass 2: host masking ─────────────────────────────────────",
              flush=True)
        shared_host_kmers = genomes_dir / "_shared_host_kmers.fna"
        host_ref_files: list[Path] = []

        # k-mers shared between ≥2 host sequences (same-family cross-hits)
        if len(host_fnas) > 1:
            r = _kmer_shared(combined_hosts, shared_host_kmers, threads)
            if r:
                host_ref_files.append(r)

        # Also mask against all pathogen CDS (cross-kingdom hits)
        if combined_pathogens.exists():
            host_ref_files.append(combined_pathogens)

        if host_ref_files:
            print(f"    BBDuk: masking hosts "
                  f"(ref={[p.name for p in host_ref_files]}) …", flush=True)
            _run_bbduk(combined_hosts, masked_hosts, host_ref_files, threads)
        else:
            print("    No reference available; copying hosts unmasked.", flush=True)
            shutil.copy2(str(combined_hosts), str(masked_hosts))
    elif masked_hosts.exists():
        print(f"  _host_masked.fna already exists "
              f"({masked_hosts.stat().st_size/1e9:.2f} GB)", flush=True)

    return (
        masked_pathogens if masked_pathogens.exists() else None,
        masked_hosts     if masked_hosts.exists()     else None,
    )


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
    ap.add_argument("--db-path", default=str(DB_PATH),
                    help="Path to phibase_db.json")
    ap.add_argument("--ref-screen", default=str(REF_SCREEN),
                    help="Path to ref_screen.tsv (from screen_references.py)")
    ap.add_argument("--db-dir", default=str(DEFAULT_DB),
                    help="Kraken2 database output directory")
    ap.add_argument("--genomes-dir", default=str(DEFAULT_GEN),
                    help="Directory for downloaded FASTAs")
    ap.add_argument("--threads", type=int, default=32,
                    help="Threads for kraken2-build --build and BBDuk (default: 32)")
    ap.add_argument("--download-only", action="store_true",
                    help="Download FASTAs only; skip BBDuk masking and kraken2-build steps")
    ap.add_argument("--build-only", action="store_true",
                    help="Skip downloads; BBDuk-mask hosts then add to library and build")
    ap.add_argument("--skip-bbduk", action="store_true",
                    help="Skip BBDuk host masking step (not recommended; produces noisy host calls)")
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

        n_downloaded = 0
        n_skipped    = 0
        n_failed     = 0
        n_added      = 0
        seen_accessions: dict[str, int] = {}

        # ── Phase 1: Download (or scan if --build-only) ───────────────────────
        # Collect all FASTAs before masking so BBDuk has the full pathogen set.
        downloaded: list[tuple[int, str, list[Path]]] = []  # (tag_taxid, kingdom, fnas)

        for i, (taxid, kingdom) in enumerate(sorted(seeds.items()), 1):
            name = t2n.get(str(taxid), str(taxid))
            print(f"\n[{i}/{len(seeds)}] taxid={taxid} ({kingdom}) {name}", flush=True)

            info      = ref.get(taxid, {})
            accession = info.get("accession", "")
            fasta_type = info.get("fasta_type", "genome")

            if not accession:
                print("  SKIP: no assembly in ref_screen.tsv", flush=True)
                n_skipped += 1
                continue

            if fasta_type == "genome":
                print("  SKIP: no annotation — zero detections in data", flush=True)
                n_skipped += 1
                continue

            if accession in seen_accessions:
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
                    tag_taxid = read_assembly_taxid(taxon_dir, taxid)
                    if tag_taxid != taxid:
                        print(f"  NOTE: taxid corrected {taxid} → {tag_taxid} "
                              f"(stale PHI-base taxid)", flush=True)
                    for fna in fnas:
                        tag_fasta_headers(fna, tag_taxid)
                    seen_accessions[accession] = tag_taxid
                    n_downloaded += 1
            else:
                fnas = list(taxon_dir.glob("**/*.fna"))
                tag_taxid = taxid
                if fnas:
                    seen_accessions[accession] = taxid

            if not fnas:
                n_failed += 1
                continue

            downloaded.append((tag_taxid, kingdom, fnas))

        print(f"\n── Download summary ──────────────────────────────────────────")
        print(f"  Seeds total:     {len(seeds)}")
        print(f"  Downloaded:      {n_downloaded}")
        print(f"  Skipped:         {n_skipped}")
        print(f"  Failed:          {n_failed}")

        if args.upload_to_acacia:
            upload_to_acacia(genomes_dir, "kraken_transcriptomes")

        if args.download_only:
            print("\n--download-only: skipped BBDuk masking and kraken2-build steps.")
            if args.upload_to_acacia:
                print(f"To retrieve on Setonix:")
                print(f"  aws s3 sync s3://{ACACIA_BUCKET}/kraken_transcriptomes/ "
                      f"kraken/cds_from_genomic/ "
                      f"--profile acacia --endpoint-url {ACACIA_ENDPOINT}")
            return

        # ── Phase 2: BBDuk masking ─────────────────────────────────────────────
        masked_pathogen_fna: Path | None = None
        masked_host_fna:     Path | None = None
        if not args.skip_bbduk:
            print(f"\n── BBDuk bidirectional masking ───────────────────────────────")
            masked_pathogen_fna, masked_host_fna = bbduk_mask_sequences(
                genomes_dir,
                threads=min(args.threads, 16),  # BBDuk doesn't benefit from 32+ threads
            )
            # Report size reduction from masking
            for label, pre, masked in [
                ("pathogens", genomes_dir / "_pathogen_combined.fna", masked_pathogen_fna),
                ("hosts",     genomes_dir / "_host_combined.fna",     masked_host_fna),
            ]:
                if pre.exists() and masked and masked.exists():
                    pre_mb    = pre.stat().st_size / 1e6
                    masked_mb = masked.stat().st_size / 1e6
                    pct       = 100 * (1 - masked_mb / pre_mb) if pre_mb else 0
                    print(f"  {label}: {pre_mb:.0f} MB → {masked_mb:.0f} MB "
                          f"({pct:.1f}% masked)", flush=True)
        else:
            print("\n--skip-bbduk: sequences added unmasked (not recommended)", flush=True)

        # ── Phase 3: Clear library + add masked sequences ──────────────────────
        print(f"\n── Adding to Kraken2 library ────────────────────────────────────")

        # Clear existing library so old unmasked files don't persist across rebuilds
        lib_dir = db_dir / "library"
        if lib_dir.exists() and list(lib_dir.iterdir()):
            print(f"  Clearing existing library ({lib_dir}) …", flush=True)
            shutil.rmtree(lib_dir)

        # Add pathogen sequences: masked combined file (preferred) or individual (fallback)
        if masked_pathogen_fna and masked_pathogen_fna.exists():
            print(f"  Adding BBDuk-masked pathogens: {masked_pathogen_fna.name}", flush=True)
            if add_to_library(masked_pathogen_fna, db_dir):
                n_added += 1
        else:
            print("  Adding unmasked pathogen FASTAs individually …", flush=True)
            for tag_taxid, kingdom, fnas in downloaded:
                if kingdom == "host":
                    continue
                for fna in fnas:
                    if add_to_library(fna, db_dir):
                        n_added += 1

        # Add host sequences: masked combined file (preferred) or individual (fallback)
        if masked_host_fna and masked_host_fna.exists():
            print(f"  Adding BBDuk-masked hosts: {masked_host_fna.name}", flush=True)
            if add_to_library(masked_host_fna, db_dir):
                n_added += 1
        else:
            print("  Adding unmasked host FASTAs individually …", flush=True)
            for tag_taxid, kingdom, fnas in downloaded:
                if kingdom != "host":
                    continue
                for fna in fnas:
                    if add_to_library(fna, db_dir):
                        n_added += 1

        print(f"\n── Library summary ──────────────────────────────────────────────")
        print(f"  FASTAs added to library: {n_added}")

        # ── Phase 4: Download taxonomy + build ────────────────────────────────
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
