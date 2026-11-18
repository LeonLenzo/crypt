#!/usr/bin/env python3
"""freshen_reads.py - keep the cohort reads out of the /scratch 21-day purge.

    python release/01_cohort/scripts/freshen_reads.py --report
    python release/01_cohort/scripts/freshen_reads.py --all
    python release/01_cohort/scripts/freshen_reads.py --all --older-than 7 --dry-run

/scratch is Lustre mounted WITHOUT noatime, so the purge clock runs on atime and reading one
4 KB block refreshes a file exactly as well as reading all of it. A full pass over the cohort
is ~12 MB of I/O, not 5 TB.

## WHY THIS FILE EXISTS RATHER THAN research/02_classify/03_cohort/freshen.py

That script took its reads root from `_layout.READS`, which resolves to
`research/02_classify/data/bronze/reads`. The live cohort moved to
`release/01_cohort/raw/reads` when the release tree was built, so the old script now points at
8.24 TB of superseded reads that are MEANT to expire and protects none of the 4.90 TB that
matters. It would run, report a healthy tree, and exit 0.

That is the SECOND time the same bug has bitten: its own docstring describes the first, when
`migrate_reads.py` moved the cohort out from under a hardcoded path on 2026-10-20. A keep-alive
that silently covers the wrong files is worse than none, because the report looks healthy.

Two changes stop it recurring:

1. **The root is derived from this file's own location**, not imported from a layout module
   that something else is free to re-point. This script sits in `release/01_cohort/scripts/`
   and freshens `release/01_cohort/raw/reads/` - the same directory `fetch_reads.py` writes to.
   If the module moves, the script moves with it and the path follows.

2. **Coverage against the cohort is checked, and a low figure is a FAILURE, not a note.** The
   first bug was not an empty directory; it was a populated directory holding the wrong files.
   So "did I find files" cannot be the test. The test is what fraction of the runs in
   `screen_samples.tsv` were actually found and touched, and below `--min-coverage` the script
   exits non-zero and says what it found instead.
"""

import argparse
import csv
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
READS = MODULE / "raw" / "reads"
SCREEN = MODULE / "data" / "curation" / "screen_samples.tsv"

PURGE_DAYS = 21
BLOCK = 4096


def cohort_runs():
    """The runs the cohort expects, as the denominator for the coverage check."""
    with open(SCREEN, newline="", encoding="utf-8") as fh:
        return {r["run"] for r in csv.DictReader(fh, delimiter="\t") if r.get("run")}


def scan(root):
    """{run: {"BioProject", "files", "bytes", "atime_days"}} from the filesystem alone.

    Presence is a file on disk, never a row in a manifest - the same rule fetch_reads.py
    follows, and for the same reason: a stale manifest would hide exactly the drift this
    script exists to catch.
    """
    now = time.time()
    tree = {}
    for p in root.rglob("*.fastq.gz"):
        run = p.name.split("_")[0].replace(".fastq.gz", "")
        e = tree.setdefault(run, {"BioProject": p.parent.name, "files": [],
                                  "bytes": 0, "atime_days": 0.0})
        st = p.stat()
        e["files"].append(p)
        e["bytes"] += st.st_size
        e["atime_days"] = max(e["atime_days"], (now - st.st_atime) / 86400)
    return tree


def report(tree, expected):
    nfiles = sum(len(e["files"]) for e in tree.values())
    nbytes = sum(e["bytes"] for e in tree.values())
    print(f"tree: {READS}")
    print(f"{len(tree)} runs, {nfiles} files, {nbytes / 1024**4:.2f} TiB "
          f"({len(tree) / len(expected) * 100:.0f}% of the {len(expected)}-run cohort)\n")
    if not tree:
        print("NOTHING on disk at that path. Run this on Setonix, not locally.")
        return

    buckets = {}
    for e in tree.values():
        buckets.setdefault(int(e["atime_days"]), [0, 0])
        buckets[int(e["atime_days"])][0] += 1
        buckets[int(e["atime_days"])][1] += e["bytes"]
    print("  atime age (days)   runs      size   days until purge")
    for age in sorted(buckets):
        n, b = buckets[age]
        left = PURGE_DAYS - age
        flag = ("  <-- EXPIRED" if left <= 0
                else "  <-- within a week" if left <= 7 else "")
        print(f"  {age:14d} {n:6d}  {b / 1e12:6.2f} TB   {left:8d}{flag}")

    at_risk = [(e["atime_days"], r) for r, e in tree.items()
               if PURGE_DAYS - e["atime_days"] <= 7]
    if at_risk:
        tb = sum(tree[r]["bytes"] for _, r in at_risk) / 1e12
        print(f"\n{len(at_risk)} runs ({tb:.2f} TB) expire within 7 days.")


def freshen(tree, want, older_than, dry_run):
    """Read one block from each file of each wanted run older than `older_than` days."""
    now = time.time()
    touched = skipped = failed = 0
    nbytes = 0
    for run in sorted(want):
        for p in tree[run]["files"]:
            try:
                age = (now - p.stat().st_atime) / 86400
            except OSError as e:
                print(f"  FAILED stat {p.name}: {e}", file=sys.stderr)
                failed += 1
                continue
            if age < older_than:
                skipped += 1
                continue
            if dry_run:
                touched += 1
                continue
            try:
                with open(p, "rb") as fh:
                    nbytes += len(fh.read(BLOCK))
                touched += 1
            except OSError as e:
                print(f"  FAILED {p.name}: {e}", file=sys.stderr)
                failed += 1
    verb = "would freshen" if dry_run else "freshened"
    print(f"{verb} {touched} files over {len(want)} runs "
          f"({nbytes / 1e6:.1f} MB read, {skipped} already fresh, {failed} failed)")
    return failed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="every cohort run on disk")
    ap.add_argument("--older-than", type=float, default=0.0,
                    help="only freshen files whose atime exceeds this many days (default 0: "
                         "everything, so nothing expires on a staggered clock)")
    ap.add_argument("--min-coverage", type=float, default=50.0,
                    help="fail if fewer than this %% of cohort runs are found on disk "
                         "(default 50). This is the guard against freshening the wrong tree.")
    ap.add_argument("--report", action="store_true", help="show the age histogram and exit")
    ap.add_argument("--dry-run", action="store_true", help="count what would be read, read nothing")
    args = ap.parse_args()

    if not READS.is_dir():
        sys.exit(f"no reads tree at {READS} - this runs on Setonix, not locally")

    expected = cohort_runs()
    tree = scan(READS)
    report(tree, expected)

    if args.report:
        return

    if not args.all:
        print("\nNothing freshened. Pass --all to act.")
        return

    found = len(set(tree) & expected)
    coverage = found / len(expected) * 100
    print()
    rc = freshen(tree, set(tree), args.older_than, args.dry_run)

    # The coverage check is deliberately AFTER the work: if the tree is right but partial,
    # freshening what is there is still correct and we want it done before the complaint.
    print(f"\ncoverage: {found}/{len(expected)} cohort runs ({coverage:.0f}%)")
    if coverage < args.min_coverage:
        print(f"ABORT-WORTHY: coverage {coverage:.0f}% is below --min-coverage "
              f"{args.min_coverage:.0f}%. Either the fetch is far from done, or this is "
              f"pointing at the wrong tree. Check that {READS} is where fetch_reads.py "
              f"is writing.", file=sys.stderr)
        sys.exit(1)
    sys.exit(1 if rc else 0)


if __name__ == "__main__":
    main()
