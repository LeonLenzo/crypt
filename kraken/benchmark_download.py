#!/usr/bin/env python3 -u
"""
benchmark_download.py — measure read acquisition throughput for different download methods.

Tests methods in order to find the fastest strategy for streaming reads into Kraken2.
Run this on each target environment (Setonix, local) and compare results.

Methods tested:
  ena-ftp    ENA FTP via curl (current classify.py approach)
  ena-https  ENA HTTPS via curl (avoids FTP quirks; same data, different protocol)
  sra-s3     SRA open-access S3 via aws-cli  (requires: aws configure or IAM role)
  prefetch   NCBI prefetch + fasterq-dump    (requires: sra-toolkit in PATH)

Usage:
  python kraken/benchmark_download.py --runs 20 --n-reads 500000 --workers 1 4 8
  python kraken/benchmark_download.py --runs 20 --methods ena-ftp ena-https --workers 1 8
  python kraken/benchmark_download.py --run-list kraken/control/output/data/run_ids.txt \\
      --runs 20 --methods ena-ftp sra-s3
"""

import argparse
import csv
import json
import os
import random
import subprocess
import sys
import tempfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT        = Path(__file__).resolve().parent.parent
CONTROL_TSV = ROOT / "kraken/control/output/data/control_runs.tsv"
SCRATCH_DEF = Path(os.environ.get("MYSCRATCH", tempfile.gettempdir())) / "dl_bench"

RESULT_COLS = [
    "method", "workers", "run", "elapsed_s", "mb_downloaded",
    "mb_per_s", "records", "records_per_s", "success",
]

SEED = 42


# ── ENA metadata ──────────────────────────────────────────────────────────────

def _ena_urls(run: str) -> dict:
    """Return {'ftp': url, 'https': url} for R1, or empty dict."""
    api = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={run}&result=read_run&fields=fastq_ftp,fastq_bytes&format=tsv")
    req = urllib.request.Request(api, headers={"User-Agent": "crypt/benchmark"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            body = r.read().decode()
    except Exception:
        return {}
    lines = [l for l in body.strip().split("\n") if l]
    if len(lines) < 2:
        return {}
    fields = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    ftp_urls = [u.strip() for u in fields.get("fastq_ftp", "").split(";") if u.strip()]
    if not ftp_urls:
        return {}
    ftp = f"ftp://{ftp_urls[0]}"
    https = ftp.replace("ftp://ftp.sra.ebi.ac.uk", "https://ftp.sra.ebi.ac.uk", 1)
    return {"ftp": ftp, "https": https}


def _sra_s3_path(run: str) -> str:
    """Construct S3 URI for a run on NCBI SRA open-access bucket."""
    return f"s3://sra-pub-run-odp/sra/{run}/{run}"


# ── Download methods ──────────────────────────────────────────────────────────

def _ena_curl(url: str, n_reads: int, dest: Path, protocol: str) -> tuple[bool, float, int]:
    """Download via curl. Returns (success, mb_downloaded, records)."""
    t0 = time.time()
    cmd = ["bash", "-c",
           f'curl --silent --fail --max-time 600 "{url}" | gunzip -c | head -n {n_reads * 4}']
    try:
        with open(dest, "w") as out:
            subprocess.run(cmd, stdout=out, timeout=660, check=False)
        if not (dest.exists() and dest.stat().st_size > 0):
            return False, 0.0, 0
        mb = dest.stat().st_size / 1e6
        records = sum(1 for _ in open(dest)) // 4
        return True, mb, records
    except Exception:
        return False, 0.0, 0
    finally:
        elapsed = time.time() - t0
        return (dest.exists() and dest.stat().st_size > 0), \
               (dest.stat().st_size / 1e6 if dest.exists() else 0), \
               (sum(1 for _ in open(dest)) // 4 if dest.exists() else 0)


def _method_ena(run: str, n_reads: int, tmp: Path, protocol: str = "ftp") -> dict:
    urls = _ena_urls(run)
    url  = urls.get(protocol)
    if not url:
        return {"success": False, "mb_downloaded": 0, "records": 0, "elapsed_s": 0}
    dest = tmp / f"{run}_r1.fastq"
    t0 = time.time()
    cmd = ["bash", "-c",
           f'curl --silent --fail --max-time 600 "{url}" | gunzip -c | head -n {n_reads * 4}']
    try:
        with open(dest, "w") as out:
            subprocess.run(cmd, stdout=out, timeout=660, check=False)
    except Exception:
        pass
    elapsed = time.time() - t0
    ok = dest.exists() and dest.stat().st_size > 0
    mb = dest.stat().st_size / 1e6 if ok else 0
    records = sum(1 for _ in open(dest)) // 4 if ok else 0
    dest.unlink(missing_ok=True)
    return {"success": ok, "mb_downloaded": round(mb, 3),
            "records": records, "elapsed_s": round(elapsed, 2)}


def _method_sra_s3(run: str, n_reads: int, tmp: Path) -> dict:
    """Download from SRA S3 open-access bucket via aws-cli, then fasterq-dump or vdb-dump."""
    s3_path = _sra_s3_path(run)
    sra_file = tmp / f"{run}.sra"
    dest     = tmp / f"{run}_r1.fastq"
    t0 = time.time()
    try:
        # download SRA file from S3
        r = subprocess.run(
            ["aws", "s3", "cp", "--no-sign-request", s3_path, str(sra_file)],
            capture_output=True, text=True, timeout=300,
        )
        if r.returncode != 0 or not sra_file.exists():
            return {"success": False, "mb_downloaded": 0, "records": 0,
                    "elapsed_s": round(time.time() - t0, 2)}

        # convert to FASTQ, take first n_reads
        cmd = ["bash", "-c",
               f'fasterq-dump --stdout --split-spot -e 4 "{sra_file}" 2>/dev/null '
               f'| head -n {n_reads * 4}']
        with open(dest, "w") as out:
            subprocess.run(cmd, stdout=out, timeout=300, check=False)
    except FileNotFoundError:
        return {"success": False, "mb_downloaded": 0, "records": 0,
                "elapsed_s": round(time.time() - t0, 2),
                "_note": "aws or fasterq-dump not in PATH"}
    except Exception:
        pass
    finally:
        sra_file.unlink(missing_ok=True)

    elapsed = time.time() - t0
    ok = dest.exists() and dest.stat().st_size > 0
    mb = dest.stat().st_size / 1e6 if ok else 0
    records = sum(1 for _ in open(dest)) // 4 if ok else 0
    dest.unlink(missing_ok=True)
    return {"success": ok, "mb_downloaded": round(mb, 3),
            "records": records, "elapsed_s": round(elapsed, 2)}


def _method_prefetch(run: str, n_reads: int, tmp: Path) -> dict:
    """NCBI prefetch (SRA toolkit) + fasterq-dump."""
    t0 = time.time()
    try:
        r = subprocess.run(
            ["prefetch", "--max-size", "20g", "-O", str(tmp), run],
            capture_output=True, text=True, timeout=300,
        )
        sra_file = tmp / run / f"{run}.sra"
        if not sra_file.exists():
            sra_file = tmp / f"{run}.sra"

        dest = tmp / f"{run}_r1.fastq"
        cmd  = ["bash", "-c",
                f'fasterq-dump --stdout --split-spot -e 4 "{sra_file}" 2>/dev/null '
                f'| head -n {n_reads * 4}']
        with open(dest, "w") as out:
            subprocess.run(cmd, stdout=out, timeout=300, check=False)
    except FileNotFoundError:
        return {"success": False, "mb_downloaded": 0, "records": 0,
                "elapsed_s": round(time.time() - t0, 2),
                "_note": "prefetch or fasterq-dump not in PATH"}
    except Exception:
        pass
    finally:
        for p in tmp.glob(f"{run}*"):
            if p.is_file():
                p.unlink(missing_ok=True)

    elapsed = time.time() - t0
    ok = dest.exists() and dest.stat().st_size > 0 if "dest" in dir() else False
    mb = dest.stat().st_size / 1e6 if ok else 0
    records = sum(1 for _ in open(dest)) // 4 if ok else 0
    if ok:
        dest.unlink(missing_ok=True)
    return {"success": ok, "mb_downloaded": round(mb, 3),
            "records": records, "elapsed_s": round(elapsed, 2)}


METHODS = {
    "ena-ftp":   lambda run, n, tmp: _method_ena(run, n, tmp, "ftp"),
    "ena-https": lambda run, n, tmp: _method_ena(run, n, tmp, "https"),
    "sra-s3":    _method_sra_s3,
    "prefetch":  _method_prefetch,
}


# ── Worker ────────────────────────────────────────────────────────────────────

def _bench_run(method: str, run: str, n_reads: int, tmp_dir: Path) -> dict:
    fn  = METHODS[method]
    res = fn(run, n_reads, tmp_dir)
    elapsed = res["elapsed_s"]
    mb      = res["mb_downloaded"]
    records = res["records"]
    return {
        "method":        method,
        "run":           run,
        "elapsed_s":     elapsed,
        "mb_downloaded": mb,
        "mb_per_s":      round(mb / elapsed, 3) if elapsed > 0 else 0,
        "records":       records,
        "records_per_s": round(records / elapsed, 1) if elapsed > 0 else 0,
        "success":       res["success"],
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-list", metavar="PATH",
                    help="Flat file of run IDs to sample from (default: control_runs.tsv)")
    ap.add_argument("--runs",     type=int, default=10,
                    help="Number of runs to benchmark per method×workers combination (default 10)")
    ap.add_argument("--n-reads",  type=int, default=500_000,
                    help="Reads to acquire per run (default 500,000 — same as classify.py)")
    ap.add_argument("--methods",  nargs="+",
                    choices=list(METHODS), default=["ena-ftp", "ena-https"],
                    help="Methods to benchmark (default: ena-ftp ena-https)")
    ap.add_argument("--workers",  nargs="+", type=int, default=[1, 4, 8],
                    help="Parallel worker counts to test (default: 1 4 8)")
    ap.add_argument("--tmp-dir",  default=str(SCRATCH_DEF))
    ap.add_argument("--out",      default="kraken/output/benchmark_download.tsv",
                    help="Output TSV path (default: kraken/output/benchmark_download.tsv)")
    args = ap.parse_args()

    # load run IDs
    if args.run_list:
        runs = [l.strip() for l in open(args.run_list) if l.strip()]
    elif CONTROL_TSV.exists():
        with open(CONTROL_TSV) as f:
            runs = [row["Run"] for row in csv.DictReader(f, delimiter="\t")]
    else:
        sys.exit(f"ERROR: no --run-list and {CONTROL_TSV} not found")

    random.seed(SEED)
    random.shuffle(runs)
    test_runs = runs[:args.runs]
    print(f"Benchmark: {len(test_runs)} runs × {args.methods} × workers={args.workers}")
    print(f"n_reads={args.n_reads:,}  tmp={args.tmp_dir}\n")

    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_results = []

    for method in args.methods:
        for n_workers in args.workers:
            print(f"── {method}  workers={n_workers} ──────────────────────────")
            t_wall = time.time()
            run_results = []

            with ThreadPoolExecutor(max_workers=n_workers) as pool:
                futs = {
                    pool.submit(_bench_run, method, run, args.n_reads, tmp_dir): run
                    for run in test_runs
                }
                for i, fut in enumerate(as_completed(futs), 1):
                    res = fut.result()
                    res["workers"] = n_workers
                    run_results.append(res)
                    status = (f"{res['mb_per_s']:.1f} MB/s  {res['records_per_s']:.0f} rec/s"
                              if res["success"] else f"FAIL")
                    print(f"  [{i:>2}/{len(test_runs)}] {res['run']}  {res['elapsed_s']:.1f}s  {status}")

            wall = time.time() - t_wall
            ok   = [r for r in run_results if r["success"]]
            if ok:
                avg_mb  = sum(r["mb_per_s"] for r in ok) / len(ok)
                avg_rec = sum(r["records_per_s"] for r in ok) / len(ok)
                agg_thr = sum(r["mb_downloaded"] for r in ok) / wall
                print(f"\n  mean per-run: {avg_mb:.2f} MB/s  {avg_rec:.0f} rec/s")
                print(f"  aggregate (wall): {agg_thr:.2f} MB/s  ({n_workers} workers, {wall:.0f}s)\n")
            else:
                print("  no successful downloads\n")

            all_results.extend(run_results)

    # write TSV
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["workers"] + RESULT_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(all_results)
    print(f"Results written to {out_path}")

    # summary table
    print("\n── Summary (mean MB/s per method × workers, successful runs only) ──")
    print(f"{'method':<12}  {'workers':>7}  {'n_ok':>5}  {'mean_MB/s':>10}  {'agg_MB/s':>10}")
    from itertools import groupby
    keyfn = lambda r: (r["method"], r["workers"])
    for (method, nw), grp in groupby(sorted(all_results, key=keyfn), key=keyfn):
        rows = list(grp)
        ok   = [r for r in rows if r["success"]]
        if not ok:
            print(f"{method:<12}  {nw:>7}  {'0':>5}  {'—':>10}  {'—':>10}")
            continue
        mean_mb = sum(r["mb_per_s"] for r in ok) / len(ok)
        # aggregate not available here; skip
        print(f"{method:<12}  {nw:>7}  {len(ok):>5}  {mean_mb:>10.2f}")


if __name__ == "__main__":
    main()
