#!/usr/bin/env python3
"""Run one arm of the database ladder over the cohort reads on disk, a BioProject at a time.

    python release/03_classify/scripts/classify.py --arm db1p --list
    python release/03_classify/scripts/classify.py --arm db1p --project PRJEB39201 --threads 8
    python release/03_classify/scripts/classify.py --arm db2 --project PRJEB39201 --force
    python release/03_classify/scripts/classify.py --arm db1p --rescan

Reads  ../01_cohort/raw/reads/<BioProject>/<Run>*.fastq.gz   the cohort reads, by SCANNING
       ../02_reference/data/db/<arm>/                        the Kraken2 index for the arm
Writes data/raw/reports/<arm>/<BioProject>/<Run>.report.tsv  one report per run
       data/raw/reports_manifest_<arm>.tsv                   run -> status, bytes, rows, TRACKED

The settings are fixed and are the same in both arms:

    kraken2 --confidence 0.15 --minimum-hit-groups 3 --report-minimizer-data

`--report-minimizer-data` is not optional. Its `distinct` column is the numerator of both
coverage and saturation, and without it every detection metric is None. The failure mode is a
column-count change rather than an error, which is why `_reports.parse_report` tells the two
layouts apart by counting columns.

`--confidence 0.15 --minimum-hit-groups 3` is db_v3's published operating point, carried
unchanged, so a difference between arms is attributable to the database and not to the caller.

## What is classified: everything on disk, and nothing is selected

No depth filter and no selection of any kind. A detection is valid at any depth; only ABSENCE
is depth-dependent, and that caveat attaches to a denominator rather than to a call. A pathogen
recovered from a 0.019 Gkmer library had to be abundant to appear at all, which makes it a
stronger observation and not a weaker one. A filter here would be a researcher degree of
freedom applied before the data is seen.

## Presence is a NON-EMPTY REPORT IN ITS FINAL LOCATION

Never a name in a table, and never existence alone. Both halves were paid for:

  * A manifest is a snapshot. Deciding what to classify from `reads_manifest.tsv` once left a
    job reporting "no reads on disk" and exiting 0 against 120 runs sitting right there, because
    the manifest predated the fetch. A stale snapshot fails SILENTLY in the safe-looking
    direction: the job succeeds at doing nothing. The tree is therefore scanned, which costs one
    stat per file and cannot drift.
  * A kraken2 process killed by a timeout leaves a short, parseable, wrong file. 224 reports
    were truncated that way in 2026-09 and consumed without complaint, so `done()` tests size
    as well as existence.

Reports are written to a per-job temporary name and renamed on success. The per-job temporary
directory is a correctness requirement: the whole point of the resume check is that a job may
start while another is still running, and two jobs sharing one partial path would have two
kraken2 processes writing one file and both renaming it over the real report. The result is a
full-size report of interleaved output, which no size check can catch.

## Why one job per project and not a per-run array

The hash is 57.5 GB for db1p and 59.8 GB for db2, so the expensive operation is LOADING the
database, not classifying against it. A per-run array task pays one load per run. A long-lived
job stages the hash into /dev/shm once and every worker mmaps that single RAM-resident copy.

`--mmap` is correct ONLY against /dev/shm. Against Lustre it is catastrophic and measured: 13
minutes, zero reports, kraken2 pinned at 1.6% CPU, because every hash lookup becomes a
filesystem RPC. It is not defaulted on; the SLURM script decides.
"""

import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _reports import SPECIES_RANK, parse_report  # noqa: E402

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
RELEASE = MODULE.parent
READS = RELEASE / "01_cohort" / "raw" / "reads"
DB_DIR = RELEASE / "02_reference" / "data" / "db"
REPORTS = MODULE / "data" / "raw" / "reports"

ARMS = ("db1p", "db2")
CONFIDENCE = "0.15"
MIN_HIT_GROUPS = "3"

# A report smaller than this is a crash, not a sample: even a run with no fungal reads writes
# the unclassified and root lines.
MIN_REPORT_BYTES = 200

MANIFEST_COLS = ["run", "BioProject", "arm", "status", "seconds", "report_bytes",
                 "species_rows", "reads_classified"]


def scan_reads(reads: Path) -> dict:
    """{run: {'BioProject', 'files', 'bytes'}} from the filesystem. Never from a manifest."""
    tree = {}
    if not reads.is_dir():
        sys.exit(f"no reads tree at {reads}")
    for proj in sorted(p for p in reads.iterdir() if p.is_dir()):
        for f in sorted(proj.glob("*.fastq.gz")):
            if f.name.endswith(".part") or f.stat().st_size == 0:
                continue
            run = f.name.split("_")[0].split(".")[0]
            e = tree.setdefault(run, {"BioProject": proj.name, "files": [], "bytes": 0})
            e["files"].append(f)
            e["bytes"] += f.stat().st_size
    return tree


def done(report: Path) -> bool:
    """Presence is a non-empty report. See the module docstring."""
    return report.exists() and report.stat().st_size >= MIN_REPORT_BYTES


def report_stats(report: Path) -> tuple[int, int]:
    """Species rows and reads under them, read back off the report that was just written."""
    n_sp = reads = 0
    for rec in parse_report(report):
        if rec["rank"] == SPECIES_RANK and rec["reads_clade"]:
            n_sp += 1
            reads += rec["reads_clade"]
    return n_sp, reads


def classify_run(run: str, m: dict, db: Path, arm: str, threads: int, mmap: bool,
                 tmp_dir: Path, reports: Path) -> dict:
    proj = m["BioProject"]
    out_dir = reports / proj
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / f"{run}.report.tsv"
    row = {"run": run, "BioProject": proj, "arm": arm, "seconds": 0, "report_bytes": 0,
           "species_rows": "", "reads_classified": ""}
    if done(report):
        n_sp, reads = report_stats(report)
        return {**row, "status": "cached", "report_bytes": report.stat().st_size,
                "species_rows": n_sp, "reads_classified": reads}

    files = list(m["files"])
    # Re-checked rather than trusted even though the scan just ran: /scratch is under a 21-day
    # purge and a job that queues behind a long one can start hours after the scan.
    missing = [f.name for f in files if not f.exists()]
    if missing:
        return {**row, "status": f"missing:{missing[0]}"}

    tmp = tmp_dir / f"{run}.report.partial"
    cmd = ["kraken2", "--db", str(db), "--threads", str(threads),
           "--confidence", CONFIDENCE, "--minimum-hit-groups", MIN_HIT_GROUPS,
           "--report-minimizer-data", "--report", str(tmp), "--output", "/dev/null"]
    if mmap:
        cmd.append("--memory-mapping")
    if len(files) == 2:
        cmd.append("--paired")
    cmd += [str(f) for f in sorted(files)]

    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
    except Exception as e:                      # noqa: BLE001 - must never kill the pool
        return {**row, "status": f"error:{type(e).__name__}", "seconds": round(time.time() - t0)}
    el = round(time.time() - t0)
    if r.returncode != 0 or not tmp.exists() or tmp.stat().st_size < MIN_REPORT_BYTES:
        tmp.unlink(missing_ok=True)
        return {**row, "status": "failed", "seconds": el,
                "reads_classified": r.stderr.strip()[-200:] if r.stderr else ""}
    tmp.replace(report)

    n_sp, reads = report_stats(report)
    print(f"  [done ] {run} {el:5d}s  {n_sp:4d} species  {reads:>12,d} reads", flush=True)
    return {**row, "status": "ok", "seconds": el, "report_bytes": report.stat().st_size,
            "species_rows": n_sp, "reads_classified": reads}


def merge_manifest(path: Path, rows: list) -> int:
    """Update the manifest for the runs this job touched, keyed by run. Rewritten, not appended.

    Appending is what the research-tree version did, and it makes the table a job log rather
    than a statement about the tree: a re-run leaves two rows for one run with no rule saying
    which is current. One row per run, last write wins, so the file is a description of what is
    on disk and a diff on it is readable.
    """
    kept = {}
    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                kept[r["run"]] = r
    for r in rows:
        kept[r["run"]] = {k: r.get(k, "") for k in MANIFEST_COLS}
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tsv.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for run in sorted(kept):
            w.writerow(kept[run])
    tmp.replace(path)
    return len(kept)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", required=True, choices=ARMS, help="which arm of the ladder")
    ap.add_argument("--project", nargs="*", default=[],
                    help="BioProjects to classify; omit to list what is on disk")
    ap.add_argument("--list", action="store_true", help="what is on disk and what is done")
    ap.add_argument("--db", help="index directory; default the arm's, override with the "
                                 "/dev/shm copy the SLURM script stages")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workers", type=int, default=8, help="runs classified concurrently")
    ap.add_argument("--mmap", action="store_true",
                    help="--memory-mapping. Correct ONLY against /dev/shm, never Lustre")
    ap.add_argument("--force", action="store_true", help="reclassify runs that have a report")
    ap.add_argument("--rescan", action="store_true",
                    help="rebuild the manifest from the reports on disk and exit")
    ap.add_argument("--reads-dir", default=str(READS))
    ap.add_argument("--reports-dir", help="default data/raw/reports/<arm>")
    ap.add_argument("--tmp-dir", default=os.environ.get("MYSCRATCH", "/tmp") + "/kraken_tmp")
    args = ap.parse_args()

    db = Path(args.db) if args.db else DB_DIR / args.arm
    reports = Path(args.reports_dir) if args.reports_dir else REPORTS / args.arm
    manifest = MODULE / "data" / "raw" / f"reports_manifest_{args.arm}.tsv"

    if args.rescan:
        # The manifest is merged per job, so two jobs finishing seconds apart can each read it
        # before the other's commit lands and the later write can drop the earlier job's rows.
        # The reports are the truth, so this rebuilds the table from them. Run it once after a
        # sweep: `seconds` is lost for rows this pass did not time, which is the only thing
        # the file holds that a report does not.
        rows = []
        for proj_dir in sorted(p for p in reports.iterdir() if p.is_dir()):
            for report in sorted(proj_dir.glob("*.report.tsv")):
                run = report.name.split(".")[0]
                n_sp, reads = report_stats(report)
                rows.append({"run": run, "BioProject": proj_dir.name, "arm": args.arm,
                             "status": "ok" if done(report) else "short",
                             "seconds": "", "report_bytes": report.stat().st_size,
                             "species_rows": n_sp, "reads_classified": reads})
        kept = {r["run"]: r for r in rows}
        if manifest.exists():
            for r in csv.DictReader(open(manifest, newline="", encoding="utf-8"),
                                    delimiter="\t"):
                if r["run"] in kept and r.get("seconds"):
                    kept[r["run"]]["seconds"] = r["seconds"]
        manifest.unlink(missing_ok=True)
        n = merge_manifest(manifest, list(kept.values()))
        short = sum(1 for r in rows if r["status"] == "short")
        print(f"-> {manifest.relative_to(MODULE)}  {n} runs from {len(rows)} reports on disk"
              + (f", {short} SHORT and unusable" if short else ""))
        return

    if not (db / "hash.k2d").is_file():
        sys.exit(f"no Kraken2 index at {db} (hash.k2d missing)")

    tree = scan_reads(Path(args.reads_dir))
    by_proj = {}
    for run, m in tree.items():
        by_proj.setdefault(m["BioProject"], []).append(run)

    if args.list or not args.project:
        print(f"arm {args.arm}  db {db}")
        print(f"{'BioProject':14s} {'onDisk':>7s} {'done':>6s} {'todo':>6s} {'TiB':>7s}")
        tot = [0, 0, 0]
        for proj in sorted(by_proj, key=lambda p: -len(by_proj[p])):
            runs = by_proj[proj]
            d = sum(1 for r in runs if done(reports / proj / f"{r}.report.tsv"))
            b = sum(tree[r]["bytes"] for r in runs) / 2**40
            print(f"{proj:14s} {len(runs):7d} {d:6d} {len(runs)-d:6d} {b:7.3f}")
            tot[0] += len(runs)
            tot[1] += d
            tot[2] += len(runs) - d
        print(f"{'TOTAL':14s} {tot[0]:7d} {tot[1]:6d} {tot[2]:6d}")
        if not args.project:
            return

    unknown = [p for p in args.project if p not in by_proj]
    if unknown:
        sys.exit(f"no reads on disk for {', '.join(unknown)}")

    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    rows, rc = [], 0
    for proj in args.project:
        runs = sorted(by_proj[proj])
        if args.force:
            for r in runs:
                (reports / proj / f"{r}.report.tsv").unlink(missing_ok=True)
        todo = [r for r in runs if not done(reports / proj / f"{r}.report.tsv")]
        print(f"\n{proj}: {len(runs)} runs on disk, {len(todo)} to classify", flush=True)
        if not todo:
            continue
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(classify_run, r, tree[r], db, args.arm, args.threads,
                              args.mmap, tmp_dir, reports): r for r in todo}
            for i, fut in enumerate(as_completed(futs), 1):
                rows.append(fut.result())
                if i % 25 == 0:
                    el = time.time() - t0
                    print(f"  [{i}/{len(todo)}] {el/60:.1f} min, "
                          f"eta {el/i*(len(todo)-i)/60:.0f} min", flush=True)
        bad = [r for r in rows if r["BioProject"] == proj
               and r["status"] not in ("ok", "cached")]
        print(f"{proj}: {sum(1 for r in rows if r['BioProject'] == proj and r['status'] == 'ok')}"
              f" classified, {len(bad)} failed, {(time.time()-t0)/60:.1f} min")
        for r in bad[:10]:
            print(f"  FAILED {r['run']}: {r['status']} {str(r['reads_classified'])[:120]}")
        if bad:
            rc = 1

    if rows:
        n = merge_manifest(manifest, rows)
        print(f"\n-> {manifest.relative_to(MODULE)}  {n} runs recorded")
    sys.exit(rc)


if __name__ == "__main__":
    main()
