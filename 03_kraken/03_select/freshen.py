#!/usr/bin/env python3
"""freshen.py — keep reads in active use out of the /scratch purge, cheaply.

Setonix /scratch is Lustre mounted WITHOUT noatime/relatime, so atime updates on any read,
and the 21-day purge clock runs on atime. Reading 4 KB of a file refreshes it exactly as well
as reading all of it: a full pass over 6,015 fastq is ~24 MB of I/O, not 14 TB.

This is a keep-alive for data in active use, not a way to park 14 TB indefinitely. Scope it
to the runs a pending analysis actually needs (--strata), and let the rest expire; every read
here is a public SRA accession, re-pullable from `run_list.tsv`.

    python 03_kraken/03_select/freshen.py --report
    python 03_kraken/03_select/freshen.py --strata --older-than 10
    python 03_kraken/03_select/freshen.py --all --dry-run
"""
import argparse, csv, collections, re, sys, time
from pathlib import Path

ROOT   = Path(__file__).resolve().parents[2]
READS  = Path("/scratch/pawsey1168/llenzo/crypt/03_kraken/03_select/data/reads")
STRATA = ROOT / "04_kallisto/01_select/data/strata.tsv"
PURGE_DAYS = 21
BLOCK = 4096
RUN_RE = re.compile(r"([A-Z]RR\d+)")


def read_files():
    """run -> [paths] for everything in the reads directory."""
    out = collections.defaultdict(list)
    if not READS.is_dir():
        sys.exit(f"reads directory not visible from here: {READS}\n"
                 "Run this on Setonix, not locally.")
    for p in READS.iterdir():
        if not p.is_file():
            continue
        m = RUN_RE.match(p.name)
        if m:
            out[m.group(1)].append(p)
    return out


def strata_runs():
    if not STRATA.exists():
        sys.exit(f"no stratified set yet: {STRATA}\n"
                 "Run: python 04_kallisto/01_select/kallisto_select.py")
    return {r["run"] for r in csv.DictReader(open(STRATA), delimiter="\t")}


def report(files):
    """Age histogram and the runs closest to the purge edge."""
    now = time.time()
    ages = {}
    for run, ps in files.items():
        ages[run] = min(int((now - p.stat().st_atime) / 86400) for p in ps)
    hist = collections.Counter(ages.values())
    print(f"{len(files)} runs, {sum(len(v) for v in files.values())} files\n")
    print("  atime age (days)   runs   days until purge")
    for age in sorted(hist):
        left = PURGE_DAYS - age
        flag = "  <-- EXPIRED" if left <= 0 else "  <-- within a week" if left <= 7 else ""
        print(f"  {age:14d}   {hist[age]:6d}   {left:8d}{flag}")
    at_risk = sorted((a, r) for r, a in ages.items() if PURGE_DAYS - a <= 7)
    if at_risk:
        print(f"\n{len(at_risk)} runs expire within 7 days:")
        for a, r in at_risk:
            print(f"  {r:14s} {a:3d} days old, {PURGE_DAYS - a:2d} left")
    return ages


def freshen(files, want, older_than, dry_run):
    """Read one block from each file of each wanted run whose atime exceeds older_than days."""
    now = time.time()
    touched = skipped = failed = 0
    nbytes = 0
    for run in sorted(want):
        for p in files.get(run, []):
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--strata", action="store_true",
                   help="only the runs in 04_kallisto/01_select/data/strata.tsv (the pending analysis)")
    g.add_argument("--all", action="store_true", help="every run present (14 TB, 6k files)")
    g.add_argument("--runs", type=Path, help="file of run accessions, one per line")
    ap.add_argument("--older-than", type=float, default=10.0,
                    help="only freshen files whose atime exceeds this many days (default 10)")
    ap.add_argument("--report", action="store_true", help="show age histogram and exit")
    ap.add_argument("--dry-run", action="store_true", help="list what would be read, read nothing")
    args = ap.parse_args()

    files = read_files()
    if args.report or not (args.strata or args.all or args.runs):
        report(files)
        if not args.report:
            print("\nNothing freshened. Pass --strata, --all or --runs to act.")
        return

    if args.strata:
        want = strata_runs() & set(files)
    elif args.runs:
        want = {l.strip() for l in open(args.runs) if l.strip()} & set(files)
    else:
        want = set(files)
    sys.exit(1 if freshen(files, want, args.older_than, args.dry_run) else 0)


if __name__ == "__main__":
    main()
