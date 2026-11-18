#!/usr/bin/env python3
"""freshen_reads.py - keep the cohort reads out of the /scratch 21-day purge.

    python release/01_cohort/scripts/freshen_reads.py --report
    python release/01_cohort/scripts/freshen_reads.py --all
    python release/01_cohort/scripts/freshen_reads.py --all --older-than 7 --dry-run

/scratch is Lustre and the purge clock runs on atime, so the job is to move atime forward
without moving 5 TB. This sets it directly with `os.utime`, and then RE-STATS every file to
confirm it moved.

## READING A FILE DOES NOT RELIABLY UPDATE ATIME HERE, WHICH IS WHY THE OLD ONE FAILED

The previous script read a 4 KB block per file on the stated reasoning that Lustre is mounted
without noatime and a short read refreshes atime as well as a long one. Measured on 2026-11-18,
that is false at scale. A single-file 4 KB read does move atime; a job reading 2,950 of them
in succession does not. Job 53349626 reported `freshened 2950 files ... 0 failed` and exited 0,
and afterwards 2,596 of 2,962 files still carried their old atime. The client appears to
coalesce or drop the atime updates it has not been forced to flush, and nothing in the read
path makes it flush them.

`os.utime` is a metadata write to the MDS instead of a side effect of a data read, so there is
nothing to coalesce. It is also cheaper: ~3,000 MDS operations and no data path at all.

## AND IT VERIFIES, BECAUSE THE OLD ONE'S REPORT WAS THE PROBLEM

Both failures here were silent successes. The first pointed at the wrong tree and said it was
healthy; the second used a mechanism that did nothing and said it had worked. So this script
re-stats every file after writing and reports how many atimes ACTUALLY moved, and a file whose
atime did not move is a failure. A keep-alive that cannot prove it kept anything alive is the
bug, not the implementation detail underneath it.

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

## POLICY

This sets atime explicitly rather than as a side effect of use, which is a more overt way of
holding data on /scratch than reading it was. It is the same intent the old script had and
failed to achieve, but it is worth being deliberate about: the cohort is ~5 TB of public ENA
accessions that are re-fetchable, and the justification for holding them is that re-fetching
costs about a day of transfer, not that they are irreplaceable. If Pawsey would rather the
space came back, the answer is to let it purge and re-fetch, not to keep this running quietly.
"""

import argparse
import csv
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
READS = MODULE / "raw" / "reads"
SCREEN = MODULE / "data" / "curation" / "screen_samples.tsv"

PURGE_DAYS = 21


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
    """Set atime forward on every wanted file, then re-stat to prove it moved.

    mtime is preserved deliberately: it is the only record of when the file was fetched, and
    the manifest's size checks are the integrity story, not the timestamp. Passing the file's
    own st_mtime back to os.utime keeps it untouched.
    """
    now = time.time()
    want_files = []
    skipped = 0
    for run in sorted(want):
        for p in tree[run]["files"]:
            try:
                st = p.stat()
            except OSError as e:
                print(f"  FAILED stat {p.name}: {e}", file=sys.stderr)
                continue
            if (now - st.st_atime) / 86400 < older_than:
                skipped += 1
                continue
            want_files.append((p, st.st_mtime))

    if dry_run:
        print(f"would freshen {len(want_files)} files over {len(want)} runs "
              f"({skipped} already fresh)")
        return 0

    written = failed = 0
    for p, mtime in want_files:
        try:
            os.utime(p, (now, mtime))     # (atime, mtime): mtime passed back unchanged
            written += 1
        except OSError as e:
            print(f"  FAILED utime {p.name}: {e}", file=sys.stderr)
            failed += 1

    # THE VERIFICATION. Both previous failures of this keep-alive were silent successes, so
    # claiming the write happened is not enough: re-stat and count what actually moved.
    stale = []
    for p, _ in want_files:
        try:
            if (time.time() - p.stat().st_atime) / 86400 > 1.0:
                stale.append(p)
        except OSError:
            stale.append(p)

    print(f"freshened {written} files over {len(want)} runs "
          f"({skipped} already fresh, {failed} failed)")
    print(f"verified  {len(want_files) - len(stale)}/{len(want_files)} atimes now current")
    if stale:
        print(f"\nDID NOT MOVE: {len(stale)} files still carry an old atime. The keep-alive "
              f"is not working and these will purge. First few:", file=sys.stderr)
        for p in stale[:10]:
            print(f"  {p}", file=sys.stderr)
    return failed + len(stale)


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
