#!/usr/bin/env python3
"""
busco_screen.py — download CDS FASTA and run BUSCO in transcriptome mode
for every candidate assembly in kraken/ref_candidates.tsv.

Downloads go to the same directory used by build.py (--genomes-dir), so
build.py will reuse them and skip re-downloading.

Run on Setonix: sbatch kraken/slurm/busco_screen.slurm
Run from crypt/:
    python kraken/busco_screen.py \\
        --busco-db-path /scratch/pawsey1168/llenzo/busco_dbs \\
        [--genomes-dir /scratch/pawsey1168/llenzo/kraken/cds_v2] \\
        [--busco-out   /scratch/pawsey1168/llenzo/kraken/busco] \\
        [--workers 16] [--cpus-per-busco 8]

Outputs:
    kraken/busco_scores.tsv          (tracked — final scores per accession)
    {busco-out}/{accession}/         (BUSCO run directories)
    {genomes-dir}/{accession}/       (downloaded CDS FNA, reused by build.py)
"""

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, make_log_dir, link_latest

CANDIDATES_TSV = Path("kraken/ref_candidates.tsv")
SCORES_TSV     = Path("kraken/busco_scores.tsv")
OUT_DIR        = Path("kraken/output/busco_screen")

SCORES_COLS = [
    "accession", "organism_name", "kingdom", "busco_lineage",
    "complete_pct", "single_pct", "duplicated_pct",
    "fragmented_pct", "missing_pct", "n_markers",
    "status",   # pass | fail | no_cds | busco_error
]

# Regex for the BUSCO summary line: C:91.2%[S:88.4%,D:2.8%],F:3.1%,M:5.7%,n:758
_BUSCO_RE = re.compile(
    r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],"
    r"F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)"
)


# ── CDS download ──────────────────────────────────────────────────────────────

def download_cds(accession: str, dest_dir: Path) -> list[Path]:
    """
    Download CDS FASTA for accession to dest_dir via NCBI datasets CLI.
    Returns list of .fna paths (empty on failure).
    Skips if dest_dir already contains .fna files (resumable).
    """
    dest_dir.mkdir(parents=True, exist_ok=True)

    existing = list(dest_dir.glob("**/*.fna"))
    if existing:
        return existing

    zip_path = dest_dir / "ncbi_dataset.zip"
    cmd = [
        "datasets", "download", "genome", "accession", accession,
        "--include", "cds",
        "--filename", str(zip_path),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0 or not zip_path.exists():
        return []

    subprocess.run(["unzip", "-q", "-o", str(zip_path), "-d", str(dest_dir)],
                   capture_output=True)
    zip_path.unlink(missing_ok=True)

    fnas = list(dest_dir.glob("**/*.fna"))
    if not fnas:
        for gz in dest_dir.glob("**/*.fna.gz"):
            out = gz.with_suffix("")
            with gzip.open(gz, "rb") as fin, open(out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            gz.unlink()
        fnas = list(dest_dir.glob("**/*.fna"))

    return fnas


def concat_fnas(fnas: list[Path], dest: Path) -> Path:
    """Concatenate multiple .fna files into a single file for BUSCO."""
    if len(fnas) == 1:
        return fnas[0]
    with open(dest, "wb") as fout:
        for fna in sorted(fnas):
            with open(fna, "rb") as fin:
                shutil.copyfileobj(fin, fout)
    return dest


# ── BUSCO run + parse ─────────────────────────────────────────────────────────

def run_busco(fna: Path, lineage: str, out_name: str,
              out_path: Path, db_path: Path, cpus: int) -> dict | None:
    """
    Run BUSCO in transcriptome mode. Returns parsed dict or None on error.
    Skips if short_summary already exists (resumable).
    """
    run_dir = out_path / out_name
    # Check for existing short_summary
    existing = list(run_dir.glob("short_summary.specific.*.txt"))
    if not existing:
        existing = list(run_dir.glob("run_*/short_summary.txt"))

    if not existing:
        run_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "busco",
            "-m", "transcriptome",
            "-i", str(fna),
            "-l", lineage,
            "--out", out_name,
            "--out_path", str(out_path),
            "--cpu", str(cpus),
            "--offline",
            "--download_path", str(db_path),
            "--force",   # overwrite incomplete runs
        ]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if r.returncode != 0:
            return None
        existing = list(run_dir.glob("short_summary.specific.*.txt"))
        if not existing:
            existing = list(run_dir.glob("run_*/short_summary.txt"))

    if not existing:
        return None

    return parse_busco_summary(existing[0])


def parse_busco_summary(summary_path: Path) -> dict | None:
    text = summary_path.read_text()
    m = _BUSCO_RE.search(text)
    if not m:
        return None
    return {
        "complete_pct":    float(m.group(1)),
        "single_pct":      float(m.group(2)),
        "duplicated_pct":  float(m.group(3)),
        "fragmented_pct":  float(m.group(4)),
        "missing_pct":     float(m.group(5)),
        "n_markers":       int(m.group(6)),
    }


# ── per-accession worker ──────────────────────────────────────────────────────

def process_one(row: dict, genomes_dir: Path, busco_out: Path,
                busco_db: Path, cpus: int,
                thresholds: dict[str, float]) -> dict:
    acc      = row["accession"]
    name     = row["organism_name"]
    kingdom  = row["kingdom"]
    lineage  = row["busco_lineage"]
    ftype    = row["fasta_type"]

    result = {
        "accession":     acc,
        "organism_name": name,
        "kingdom":       kingdom,
        "busco_lineage": lineage,
        "complete_pct":  "",
        "single_pct":    "",
        "duplicated_pct": "",
        "fragmented_pct": "",
        "missing_pct":   "",
        "n_markers":     "",
        "status":        "",
    }

    if ftype != "cds":
        result["status"] = "no_cds"
        return result

    dest_dir = genomes_dir / acc
    fnas = download_cds(acc, dest_dir)
    if not fnas:
        result["status"] = "no_cds"
        return result

    combined = dest_dir / f"{acc}_combined.fna"
    input_fna = concat_fnas(fnas, combined)

    scores = run_busco(input_fna, lineage, acc, busco_out, busco_db, cpus)
    if scores is None:
        result["status"] = "busco_error"
        return result

    result.update(scores)
    threshold = thresholds.get(kingdom, 50.0)
    result["status"] = "pass" if scores["complete_pct"] >= threshold else "fail"
    return result


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes-dir",    default="/scratch/pawsey1168/llenzo/kraken/cds_v2",
                    help="CDS download directory (same as build.py --genomes-dir)")
    ap.add_argument("--busco-out",      default="/scratch/pawsey1168/llenzo/kraken/busco",
                    help="Directory for BUSCO run output")
    ap.add_argument("--busco-db-path",  default="/scratch/pawsey1168/llenzo/busco_dbs",
                    help="Path to pre-downloaded BUSCO lineage databases")
    ap.add_argument("--scores-tsv",     default=str(SCORES_TSV),
                    help="Output scores TSV (default: kraken/busco_scores.tsv)")
    ap.add_argument("--workers",        type=int, default=16,
                    help="Parallel download+BUSCO workers")
    ap.add_argument("--cpus-per-busco", type=int, default=8,
                    help="CPUs passed to each BUSCO call")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "busco_screen.log")
    link_latest(logs_base, log_dir / "busco_screen.log")
    sys.stdout = log

    genomes_dir = Path(args.genomes_dir)
    busco_out   = Path(args.busco_out)
    busco_db    = Path(args.busco_db_path)
    scores_tsv  = Path(args.scores_tsv)

    genomes_dir.mkdir(parents=True, exist_ok=True)
    busco_out.mkdir(parents=True, exist_ok=True)

    # BUSCO thresholds — used for status column; filter_refs.py applies them
    thresholds = {"fungal": 50.0, "oomycete": 65.0}

    # Load candidates
    candidates: list[dict] = []
    with open(CANDIDATES_TSV, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            candidates.append(row)
    print(f"Candidates: {len(candidates):,}", flush=True)

    # Load existing scores for resumability
    done_accs: set[str] = set()
    if scores_tsv.exists():
        with open(scores_tsv, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                done_accs.add(row["accession"])
    print(f"Already scored: {len(done_accs):,}  remaining: {len(candidates) - len(done_accs):,}",
          flush=True)

    todo = [r for r in candidates if r["accession"] not in done_accs]
    if not todo:
        print("All accessions already scored.", flush=True)
        return

    # Open scores TSV for appending
    write_header = not scores_tsv.exists() or scores_tsv.stat().st_size == 0
    scores_fh  = open(scores_tsv, "a", newline="")
    writer     = csv.DictWriter(scores_fh, fieldnames=SCORES_COLS, delimiter="\t",
                                extrasaction="ignore")
    if write_header:
        writer.writeheader()

    completed = 0
    failed    = 0
    no_cds    = 0

    def work(row):
        return process_one(row, genomes_dir, busco_out, busco_db,
                           args.cpus_per_busco, thresholds)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(work, row): row["accession"] for row in todo}
        for fut in as_completed(futs):
            acc = futs[fut]
            try:
                result = fut.result()
            except Exception as e:
                print(f"  ERROR {acc}: {e}", flush=True)
                failed += 1
                continue

            writer.writerow(result)
            scores_fh.flush()

            status = result["status"]
            pct    = result.get("complete_pct", "")
            pct_s  = f"{pct:.1f}%" if isinstance(pct, float) else "—"
            completed += 1
            if status == "no_cds":
                no_cds += 1
            elif status == "busco_error":
                failed += 1
            total_done = completed + failed
            print(f"  [{total_done}/{len(todo)}] {acc}  {status}  C={pct_s}",
                  flush=True)

    scores_fh.close()

    n_pass = sum(1 for r in csv.DictReader(open(scores_tsv), delimiter="\t")
                 if r.get("status") == "pass")
    print(f"\n── BUSCO screen summary ─────────────────────────────────────────")
    print(f"  Processed:   {completed:,}")
    print(f"  No CDS:      {no_cds:,}")
    print(f"  BUSCO error: {failed:,}")
    print(f"  Pass:        {n_pass:,}  (across all runs including prior)")
    print(f"\nScores: {scores_tsv}")
    print(f"Next:   python kraken/filter_refs.py")


if __name__ == "__main__":
    main()
