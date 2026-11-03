#!/usr/bin/env python3
"""freshen.py — keep the cohort reads out of the /scratch purge, cheaply.

Setonix /scratch is Lustre mounted WITHOUT noatime/relatime, so atime updates on any read and
the 21-day purge clock runs on atime. Reading 4 KB of a file refreshes it exactly as well as
reading all of it: a full pass over the 2,442 fastq now on disk is ~10 MB of I/O, not 4 TB.

## Why this exists rather than the 03_kraken one

`03_kraken/03_select/freshen.py` hardcoded the FLAT reads directory. `migrate_reads.py` then
moved the cohort into `reads/<BioProject>/` and deliberately preserved atime, so the purge
clock carried across unbroken. The old script still ran, still reported, and still exited 0 --
against a tree the cohort had left. The 4 TB that matters was protected only incidentally, by
classify jobs happening to read it.

That is the failure mode worth naming: a keep-alive that silently covers the wrong files is
worse than none, because the report looks healthy. This version takes its root from
`_layout.READS`, so the path cannot drift from the one the rest of the module writes to.

## Scope: the cohort tree only

The flat tree still holds the 2,041 non-cohort STAT-gated runs. They are NOT freshened here
and are meant to expire -- every one is a public SRA accession, re-pullable from
`run_list.tsv`. Keeping them alive was justified while the db_v3-vs-db_v4 comparison was open;
it is not open any more.

This is a keep-alive for data in active use, not a way to park terabytes indefinitely.

    python 02_classify/03_cohort/freshen.py --report
    python 02_classify/03_cohort/freshen.py --all --older-than 10
    python 02_classify/03_cohort/freshen.py --project PRJEB39201 --dry-run
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _layout import READS  # noqa: E402
from _reads import scan_tree  # noqa: E402

PURGE_DAYS = 21
BLOCK = 4096


def report(tree: dict) -> None:
    """Age histogram and the runs closest to the purge edge."""
    nfiles = sum(len(e["files"]) for e in tree.values())
    nbytes = sum(e["bytes"] for e in tree.values())
    print(f"{len(tree)} runs, {nfiles} files, {nbytes / 1024**4:.2f} TiB")
    print(f"tree: {READS}\n")
    if not tree:
        print("Nothing on disk at that path. Run this on Setonix, not locally.")
        return

    buckets: dict[int, int] = {}
    for e in tree.values():
        age = int(e["atime_days"])
        buckets[age] = buckets.get(age, 0) + 1
    print("  atime age (days)   runs   days until purge")
    for age in sorted(buckets):
        left = PURGE_DAYS - age
        flag = ("  <-- EXPIRED" if left <= 0
                else "  <-- within a week" if left <= 7 else "")
        print(f"  {age:14d}   {buckets[age]:6d}   {left:8d}{flag}")

    at_risk = sorted((int(e["atime_days"]), r) for r, e in tree.items()
                     if PURGE_DAYS - e["atime_days"] <= 7)
    if at_risk:
        print(f"\n{len(at_risk)} runs expire within 7 days:")
        for age, run in at_risk[:40]:
            print(f"  {run:14s} {tree[run]['BioProject']:14s} "
                  f"{age:3d} days old, {PURGE_DAYS - age:2d} left")
        if len(at_risk) > 40:
            print(f"  ... and {len(at_risk) - 40} more")


def freshen(tree: dict, want: set, older_than: float, dry_run: bool) -> int:
    """Read one block from each file of each wanted run older than `older_than` days."""
    now = time.time()
    touched = skipped = failed = 0
    nbytes = 0
    for run in sorted(want):
        for p in tree[run]["files"]:
            age = (now - p.stat().st_atime) / 86400
            if age < older_than:
                skipped += 1
                continue
            if dry_run:
                print(f"  would freshen {p.name} ({age:.1f} days)")
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
    print(f"\n{verb} {touched} files over {len(want)} runs "
          f"({nbytes / 1e6:.1f} MB read, {skipped} already fresh, {failed} failed)")
    return failed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--all", action="store_true", help="every run in the cohort tree")
    g.add_argument("--project", action="append", metavar="PRJ",
                   help="only these BioProjects; repeatable")
    g.add_argument("--runs", type=Path, help="file of run accessions, one per line")
    ap.add_argument("--older-than", type=float, default=10.0,
                    help="only freshen files whose atime exceeds this many days (default 10)")
    ap.add_argument("--report", action="store_true", help="show age histogram and exit")
    ap.add_argument("--dry-run", action="store_true",
                    help="list what would be read, read nothing")
    args = ap.parse_args()

    tree = scan_tree(READS)
    if args.report or not (args.all or args.project or args.runs):
        report(tree)
        if not args.report:
            print("\nNothing freshened. Pass --all, --project or --runs to act.")
        return

    if args.project:
        keep = set(args.project)
        want = {r for r, e in tree.items() if e["BioProject"] in keep}
        missing = keep - {e["BioProject"] for e in tree.values()}
        if missing:
            print(f"not on disk: {', '.join(sorted(missing))}", file=sys.stderr)
    elif args.runs:
        want = {l.strip() for l in open(args.runs) if l.strip()} & set(tree)
    else:
        want = set(tree)

    if not want:
        sys.exit("nothing to freshen: no requested run is on disk")
    sys.exit(1 if freshen(tree, want, args.older_than, args.dry_run) else 0)


if __name__ == "__main__":
    main()
