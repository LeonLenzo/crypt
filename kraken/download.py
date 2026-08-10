#!/usr/bin/env python3 -u
"""
download.py — batch-download FASTQ reads from ENA HTTPS for Kraken2 classification.

Saves R1 (and R2 for paired-end) as gzipped FASTQ to --reads-dir. Resumable via
download_index.txt. Use with classify.py --reads-dir for the separated workflow.

Output naming:
  {reads-dir}/{run}_1.fastq.gz   paired R1
  {reads-dir}/{run}_2.fastq.gz   paired R2
  {reads-dir}/{run}.fastq.gz     single-end

Usage:
  python kraken/download.py \\
      --run-list kraken/control/output/data/run_ids.txt \\
      --reads-dir /scratch/pawsey1168/llenzo/kraken_reads \\
      --workers 16 --n-reads 500000
"""

import argparse
import csv
import fcntl
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

IN_DIR  = Path("stat/output/fetch_runs/data")
OUT_DIR = Path("kraken/output/download")

N_READS = 500_000
WORKERS = 16

_index_lock = threading.Lock()


# ── ENA HTTPS helpers ─────────────────────────────────────────────────────────

def _ena_https_urls(run: str) -> list:
    """Return ENA HTTPS FASTQ URLs for a run. Empty list if unavailable."""
    api = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={run}&result=read_run&fields=fastq_ftp&format=tsv")
    req = urllib.request.Request(api, headers={"User-Agent": "crypt/download"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            body = r.read().decode()
    except Exception:
        return []
    lines = [l for l in body.strip().split("\n") if l]
    if len(lines) < 2:
        return []
    ftp_field = lines[1].split("\t")[-1]
    return [f"https://{p.strip()}" for p in ftp_field.split(";") if p.strip()]


# ── Download one run ──────────────────────────────────────────────────────────

def _download_run(run: str, n_reads: int, reads_dir: Path, tmp_dir: Path) -> dict:
    """Download R1 (and R2) for one run to reads_dir as gzipped FASTQ.
    Atomic: writes to a .tmp file then renames on success."""
    result = {"run": run, "error": None, "ts": time.strftime("%Y-%m-%dT%H:%M:%S")}

    urls = _ena_https_urls(run)
    if not urls:
        result["error"] = "no_urls"
        return result

    is_paired  = len(urls) >= 2
    suffixes   = ["_1.fastq.gz", "_2.fastq.gz"] if is_paired else [".fastq.gz"]
    n_files    = 0

    for url, suffix in zip(urls[:2], suffixes):
        dest = reads_dir / f"{run}{suffix}"
        if dest.exists() and dest.stat().st_size > 0:
            n_files += 1
            continue                        # already downloaded

        tmp = tmp_dir / f"{run}{suffix}.tmp"
        cmd = ["bash", "-c",
               f'curl --silent --fail --max-time 600 "{url}" | '
               f'gunzip -c | head -n {n_reads * 4} | gzip -c']
        try:
            with open(tmp, "wb") as out:
                subprocess.run(cmd, stdout=out, timeout=660, check=False)
            if tmp.exists() and tmp.stat().st_size > 0:
                tmp.rename(dest)
                n_files += 1
            else:
                tmp.unlink(missing_ok=True)
                result["error"] = "empty_output"
                return result
        except Exception as e:
            tmp.unlink(missing_ok=True)
            result["error"] = str(e)
            return result

    result["n_files"]  = n_files
    result["is_paired"] = is_paired
    return result


# ── Index helpers ─────────────────────────────────────────────────────────────

def _load_index(out_dir: Path) -> set:
    idx = out_dir / "download_index.txt"
    return set(idx.read_text().splitlines()) if idx.exists() else set()


def _append_index(run: str, out_dir: Path) -> None:
    with _index_lock:
        with open(out_dir / "download_index.txt", "a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                f.write(run + "\n")
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)


# ── Run list loading (same logic as classify.py) ──────────────────────────────

def _load_runs(args) -> list:
    run_ids = []
    if args.run_list:
        run_ids = [l.strip() for l in open(args.run_list) if l.strip()]
        print(f"Loaded {len(run_ids):,} runs from {args.run_list}")
    elif args.runs_tsv:
        with open(args.runs_tsv) as f:
            header = f.readline().strip().split("\t")
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if args.biosample_rep and parts[col["biosample_representative"]] != "True":
                    continue
                if args.hc and parts[col["same_genus_secondary"]] != "False":
                    continue
                run_ids.append(parts[col["Run"]])
        print(f"Loaded {len(run_ids):,} runs from {args.runs_tsv}")
    else:
        modes = ["mal", "hal"] if args.mode == "both" else [args.mode]
        for mode in modes:
            p = IN_DIR / f"{mode}_runs.json"
            if not p.exists():
                print(f"WARNING: {p} not found — skipping {mode}")
                continue
            import json
            with open(p) as f:
                runs = json.load(f)
            run_ids.extend(runs.keys())
            print(f"Loaded {len(runs):,} {mode.upper()} runs from {p}")
    return run_ids


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    global N_READS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode",         choices=["mal", "hal", "both"], default="both")
    ap.add_argument("--runs-tsv",     default=None, metavar="PATH")
    ap.add_argument("--run-list",     default=None, metavar="PATH",
                    help="Plain text file, one Run per line.")
    ap.add_argument("--biosample-rep", action="store_true")
    ap.add_argument("--hc",            action="store_true")
    ap.add_argument("--reads-dir",    required=True, metavar="PATH",
                    help="Directory to save gzipped FASTQs.")
    ap.add_argument("--out-dir",      default=str(OUT_DIR))
    ap.add_argument("--tmp-dir",      default=None,
                    help="Temp dir for partial downloads (default: reads-dir/.tmp)")
    ap.add_argument("--workers",      type=int, default=WORKERS,
                    help=f"Parallel downloads (default: {WORKERS})")
    ap.add_argument("--n-reads",      type=int, default=N_READS,
                    help=f"Reads to download per run (default: {N_READS:,})")
    ap.add_argument("--limit",        type=int, default=None)
    ap.add_argument("--array-id",     type=int, default=None)
    ap.add_argument("--array-count",  type=int, default=None)
    args = ap.parse_args()

    N_READS = args.n_reads

    reads_dir = Path(args.reads_dir)
    reads_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = Path(args.tmp_dir) if args.tmp_dir else reads_dir / ".tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    out_dir  = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir  = make_log_dir(out_dir / "logs")
    log      = _Tee(log_dir / "download.log")
    link_latest(out_dir / "logs", log_dir / "download.log")
    sys.stdout = log

    all_runs = _load_runs(args)
    done     = _load_index(reads_dir)
    todo     = [r for r in all_runs if r not in done]

    if args.array_id is not None and args.array_count:
        todo = todo[args.array_id::args.array_count]
        print(f"Array task {args.array_id}/{args.array_count}: {len(todo):,} runs assigned")
    if args.limit:
        todo = todo[:args.limit]
        print(f"--limit {args.limit}: processing {len(todo):,} runs")

    print(f"\nTotal: {len(all_runs):,} | Done: {len(done):,} | Remaining: {len(todo):,}")
    print(f"n_reads={N_READS:,}  workers={args.workers}")
    print(f"reads-dir: {reads_dir}\n")

    if not todo:
        print("All runs already downloaded.")
        log.close()
        return

    n_ok = n_err = 0
    t0 = time.time()

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(_download_run, run, N_READS, reads_dir, tmp_dir): run
            for run in todo
        }
        for i, fut in enumerate(as_completed(futures), 1):
            run = futures[fut]
            try:
                result = fut.result()
            except Exception as e:
                result = {"run": run, "error": str(e)}

            if result.get("error"):
                n_err += 1
                status = f"ERR({result['error']})"
            else:
                n_ok += 1
                pe_str = "PE" if result.get("is_paired") else "SE"
                status = pe_str
                _append_index(run, reads_dir)

            elapsed = time.time() - t0
            rate    = i / elapsed if elapsed > 0 else 0
            eta     = (len(todo) - i) / rate if rate > 0 else 0
            print(f"[{i:>7}/{len(todo):,}] {run}  {status}  "
                  f"| {rate:.1f}/s  ETA {eta/3600:.1f}h", flush=True)

    elapsed = time.time() - t0
    print(f"\n── Summary ─────────────────────────────")
    print(f"  Downloaded: {n_ok:,}")
    print(f"  Errors:     {n_err:,}")
    print(f"  Elapsed:    {elapsed/3600:.1f}h")
    print(f"  Reads dir:  {reads_dir}")
    log.close()


if __name__ == "__main__":
    main()
