#!/usr/bin/env python3
"""
kraken_classify.py — stream reads from ENA FTP and classify with Kraken2.

For each run in {mode}_runs.json:
  1. Query ENA portal for FASTQ FTP URLs
  2. Stream N_READS from ENA (curl | gunzip | head) to a temp FASTQ
  3. Run kraken2 --report against the pre-built database
  4. Parse report: extract species-level detections + host %
  5. Append result to kraken_cache.jsonl (one line per run, resumable)

Run from crypt/ on Setonix (requires kraken2 in PATH):
    module load kraken2   # or equivalent on Setonix
    python kraken/classify.py --mode mal --db /scratch/kraken_db
    python kraken/classify.py --mode hal --db /scratch/kraken_db
    python kraken/classify.py --mode both --db /scratch/kraken_db

Scale: ~593k runs total (MAL ~48k + HAL ~560k). At 500k reads/run with 8
parallel workers, expect ~3-5 days on Setonix. Use Slurm array jobs for
production runs (see kraken_classify.slurm).

Reads from:  output/01_fetch_runs/data/{mode}_runs.json
Output:      output/kraken_classify/data/kraken_cache.jsonl
             output/kraken_classify/data/kraken_cache_index.txt
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, make_log_dir, link_latest

# ── Constants ─────────────────────────────────────────────────────────────────

N_READS        = 500_000    # reads to stream per run (per mate for paired-end)
CONFIDENCE     = 0.1        # kraken2 --confidence threshold
ENA_RATE       = 8          # max concurrent ENA FTP connections
KRAKEN_THREADS = 4          # kraken2 threads per worker (workers × threads ≤ total cores)
WORKERS        = 8          # parallel runs

# Temp dir: use Setonix scratch if available, else system tmp
SCRATCH = Path(os.environ.get("MYSCRATCH", tempfile.gettempdir())) / "kraken_tmp"

IN_DIR  = Path("output/01_fetch_runs/data")
OUT_DIR = Path("output/kraken_classify")

# ── ENA FTP helpers ───────────────────────────────────────────────────────────

def _ena_fastq_urls(run: str) -> list[str]:
    """Return ENA FTP FASTQ URLs for a run accession. Empty list if unavailable."""
    api = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={run}&result=read_run&fields=fastq_ftp&format=tsv")
    try:
        req = urllib.request.Request(api, headers={"User-Agent": "crypt/kraken_classify"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            body = resp.read().decode()
    except Exception:
        return []
    lines = [l for l in body.strip().split("\n") if l]
    if len(lines) < 2:
        return []
    ftp_field = lines[1].split("\t")[-1]
    return [f"ftp://{p.strip()}" for p in ftp_field.split(";") if p.strip()]


def _stream_fastq(url: str, n_reads: int, dest: Path) -> bool:
    """Stream n_reads from a gzipped ENA FTP URL into dest (plain FASTQ)."""
    cmd = ["bash", "-c",
           f'curl --silent --fail --max-time 300 "{url}" | gunzip -c | head -n {n_reads * 4}']
    try:
        with open(dest, "w") as out:
            r = subprocess.run(cmd, stdout=out, timeout=360)
        # head exits 141 (SIGPIPE) once it has enough lines — that's fine
        return dest.exists() and dest.stat().st_size > 0
    except Exception:
        return False


# ── Kraken2 runner ────────────────────────────────────────────────────────────

def _run_kraken2(db_dir: Path, reads: list[Path],
                 report: Path, threads: int) -> bool:
    """Run kraken2. reads is [r1] for SE or [r1, r2] for PE."""
    cmd = [
        "kraken2",
        "--db", str(db_dir),
        "--report", str(report),
        "--confidence", str(CONFIDENCE),
        "--threads", str(threads),
        "--output", "/dev/null",   # discard per-read output; we only need the report
    ]
    if len(reads) == 2:
        cmd += ["--paired", str(reads[0]), str(reads[1])]
    else:
        cmd.append(str(reads[0]))

    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        return r.returncode == 0
    except Exception:
        return False


def _parse_report(report: Path) -> dict:
    """
    Parse a Kraken2 report file.
    Returns {pct_classified, pct_unclassified, n_reads, species: [{taxid, name, pct, reads}]}.

    Report columns: pct  reads  taxReads  rank  taxid  name
    """
    species = []
    pct_unclassified = 0.0
    n_reads = 0

    try:
        with open(report) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                pct   = float(parts[0])
                reads = int(parts[1])
                rank  = parts[3].strip()
                taxid = int(parts[4].strip())
                name  = parts[5].strip()

                if rank == "U":                 # unclassified root
                    pct_unclassified = pct
                    n_reads = int(parts[1]) + int(parts[1])   # U + classified
                elif rank == "R" and taxid == 1:
                    # root = total classified; n_reads = classified + unclassified
                    # (we'll overwrite below with the correct total)
                    pass
                elif rank == "S" and pct > 0:
                    species.append({
                        "taxid": taxid,
                        "name":  name,
                        "pct":   round(pct, 4),
                        "reads": reads,
                    })
    except Exception:
        pass

    # Reconstruct total reads from unclassified line (reads = U count, total = U + classified)
    # Re-parse for n_reads robustly
    try:
        with open(report) as f:
            lines = f.readlines()
        u_line  = next((l for l in lines if "\tU\t" in l), None)
        r_line  = next((l for l in lines if "\tR\t1\t" in l or "\tR\t\t1\t" in l), None)
        u_reads = int(u_line.split("\t")[1]) if u_line else 0
        r_reads = int(r_line.split("\t")[1]) if r_line else 0
        n_reads = u_reads + r_reads
    except Exception:
        pass

    return {
        "pct_unclassified": round(pct_unclassified, 4),
        "pct_classified":   round(100.0 - pct_unclassified, 4),
        "n_reads":          n_reads,
        "species":          sorted(species, key=lambda x: -x["pct"]),
    }


# ── Per-run worker ────────────────────────────────────────────────────────────

def _gzip_file(src: Path, dest: Path) -> bool:
    """Gzip src → dest. Returns True on success."""
    try:
        result = subprocess.run(
            ["gzip", "-c", str(src)],
            stdout=open(dest, "wb"), timeout=120
        )
        return result.returncode == 0
    except Exception:
        return False


def _process_run(run_id: str, layout: str, db_dir: Path,
                 tmp_dir: Path, threads: int,
                 reads_dir: Path | None = None) -> dict:
    """Download, classify, and parse one run. Returns result dict.

    If reads_dir is set, gzipped FASTQs are moved there after classification
    instead of being deleted — for archival to Acacia.
    """
    result = {
        "run":    run_id,
        "layout": layout,
        "error":  None,
        "ts":     time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    urls = _ena_fastq_urls(run_id)
    if not urls:
        result["error"] = "no_ena_urls"
        return result

    run_tmp = tmp_dir / run_id
    run_tmp.mkdir(parents=True, exist_ok=True)

    try:
        is_paired = layout == "PAIRED" and len(urls) >= 2

        if is_paired:
            r1 = run_tmp / "r1.fastq"
            r2 = run_tmp / "r2.fastq"
            ok1 = _stream_fastq(urls[0], N_READS, r1)
            ok2 = _stream_fastq(urls[1], N_READS, r2)
            if not (ok1 and ok2):
                result["error"] = "stream_failed"
                return result
            reads = [r1, r2]
        else:
            r1 = run_tmp / "r1.fastq"
            ok = _stream_fastq(urls[0], N_READS, r1)
            if not ok:
                result["error"] = "stream_failed"
                return result
            reads = [r1]

        report = run_tmp / "report.txt"
        ok = _run_kraken2(db_dir, reads, report, threads)
        if not ok or not report.exists():
            result["error"] = "kraken2_failed"
            return result

        result.update(_parse_report(report))

        # archive reads if requested
        if reads_dir is not None:
            reads_dir.mkdir(parents=True, exist_ok=True)
            suffixes = ["_1.fastq.gz", "_2.fastq.gz"] if is_paired else [".fastq.gz"]
            for fastq, suffix in zip(reads, suffixes):
                dest = reads_dir / f"{run_id}{suffix}"
                if not dest.exists():
                    _gzip_file(fastq, dest)

    finally:
        for f in run_tmp.iterdir():
            f.unlink(missing_ok=True)
        run_tmp.rmdir()

    return result


# ── Cache helpers ─────────────────────────────────────────────────────────────

_cache_lock = threading.Lock()


def _load_cache_index(cache_dir: Path) -> set[str]:
    """Return set of already-processed run accessions."""
    idx = cache_dir / "kraken_cache_index.txt"
    if idx.exists():
        return set(idx.read_text().splitlines())
    return set()


def _append_cache(result: dict, cache_dir: Path) -> None:
    """Append one result to kraken_cache.jsonl and update index."""
    with _cache_lock:
        with open(cache_dir / "kraken_cache.jsonl", "a") as f:
            f.write(json.dumps(result) + "\n")
        with open(cache_dir / "kraken_cache_index.txt", "a") as f:
            f.write(result["run"] + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    global N_READS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal", "both"], default="both")
    ap.add_argument("--db", required=True, help="Path to Kraken2 database directory")
    ap.add_argument("--workers", type=int, default=WORKERS,
                    help=f"Parallel runs (default: {WORKERS})")
    ap.add_argument("--kraken-threads", type=int, default=KRAKEN_THREADS,
                    help=f"Kraken2 threads per worker (default: {KRAKEN_THREADS})")
    ap.add_argument("--n-reads", type=int, default=N_READS,
                    help=f"Reads to stream per run (default: {N_READS:,})")
    ap.add_argument("--tmp-dir", default=str(SCRATCH),
                    help=f"Temp directory for intermediate files (default: {SCRATCH})")
    ap.add_argument("--reads-dir", default=None,
                    help="If set, gzipped FASTQs are kept here after classification "
                         "(for archival to Acacia). Omit to discard reads immediately.")
    ap.add_argument("--limit", type=int, default=None,
                    help="Process at most N runs (useful for testing)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cache_dir = OUT_DIR / "data"
    cache_dir.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_classify.log")
    link_latest(logs_base, log_dir / "kraken_classify.log")
    sys.stdout = log

    db_dir  = Path(args.db)
    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    N_READS = args.n_reads

    reads_dir = Path(args.reads_dir) if args.reads_dir else None
    if reads_dir:
        reads_dir.mkdir(parents=True, exist_ok=True)
        print(f"Reads will be archived (gzipped) to: {reads_dir}")
        print(f"  Estimated storage: ~{608_000 * 55 // 1024 // 1024:.0f} TB for full run set")
        print(f"  Sync to Acacia with: aws s3 sync {reads_dir} s3://<bucket>/kraken_reads/")
    else:
        print("Reads will be discarded after classification (use --reads-dir to keep)")

    # ── Load runs ──────────────────────────────────────────────────────────────
    modes = ["mal", "hal"] if args.mode == "both" else [args.mode]
    all_runs: dict[str, str] = {}   # run_id → layout

    for mode in modes:
        runs_path = IN_DIR / f"{mode}_runs.json"
        if not runs_path.exists():
            print(f"WARNING: {runs_path} not found — skipping {mode}")
            continue
        with open(runs_path) as f:
            runs = json.load(f)
        for run_id, meta in runs.items():
            all_runs[run_id] = meta.get("LibraryLayout", "SINGLE")
        print(f"Loaded {len(runs):,} {mode.upper()} runs from {runs_path}")

    # ── Resume from cache ──────────────────────────────────────────────────────
    done = _load_cache_index(cache_dir)
    todo = [(run_id, layout) for run_id, layout in all_runs.items()
            if run_id not in done]
    if args.limit:
        todo = todo[:args.limit]
        print(f"--limit {args.limit}: processing {len(todo)} runs")

    print(f"\nTotal runs: {len(all_runs):,} | "
          f"Already done: {len(done):,} | "
          f"Remaining: {len(todo):,}")
    print(f"Settings: n_reads={N_READS:,}, confidence={CONFIDENCE}, "
          f"workers={args.workers}, kraken_threads={args.kraken_threads}")
    print(f"DB: {db_dir}")
    print(f"Tmp: {tmp_dir}\n")

    if not todo:
        print("All runs already classified. Done.")
        return

    # ── Classify ───────────────────────────────────────────────────────────────
    n_ok = n_err = 0
    t0 = time.time()

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(_process_run, run_id, layout, db_dir,
                        tmp_dir, args.kraken_threads, reads_dir): run_id
            for run_id, layout in todo
        }

        for i, fut in enumerate(as_completed(futures), 1):
            run_id = futures[fut]
            try:
                result = fut.result()
            except Exception as e:
                result = {"run": run_id, "error": str(e),
                          "ts": time.strftime("%Y-%m-%dT%H:%M:%S")}

            _append_cache(result, cache_dir)

            if result.get("error"):
                n_err += 1
                status = f"ERR({result['error']})"
            else:
                n_ok += 1
                status = f"{result.get('pct_classified', 0):.1f}% classified"

            elapsed = time.time() - t0
            rate = i / elapsed if elapsed > 0 else 0
            eta  = (len(todo) - i) / rate if rate > 0 else 0
            print(f"[{i:>7}/{len(todo):,}] {run_id}  {status}  "
                  f"| {rate:.1f}/s  ETA {eta/3600:.1f}h", flush=True)

    elapsed = time.time() - t0
    print(f"\n── Summary ────────────────────────────────────────")
    print(f"  Classified: {n_ok:,}")
    print(f"  Errors:     {n_err:,}")
    print(f"  Elapsed:    {elapsed/3600:.1f}h")
    print(f"  Cache:      {cache_dir / 'kraken_cache.jsonl'}")

    log.close()


if __name__ == "__main__":
    main()
