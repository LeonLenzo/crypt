#!/usr/bin/env python3
"""migrate_reads.py — move the flat 03_kraken reads tree into per-BioProject directories.

RUN THIS ON SETONIX. It is a no-op anywhere else, because the source tree only exists there.

    python 02_classify/03_cohort/migrate_reads.py              # plan only, changes nothing
    python 02_classify/03_cohort/migrate_reads.py --apply
    python 02_classify/03_cohort/migrate_reads.py --apply --clean-orphans

## What moves, and what deliberately does not

    03_kraken/03_select/data/reads/<Run>_1.fastq.gz            5,989 files, one directory
        -> 02_classify/data/bronze/reads/<BioProject>/<Run>_1.fastq.gz

**Only the 1,184 cohort runs move.** The other 2,041 runs on scratch are the old STAT-gated
selection and they scatter across 219 BioProjects, so moving them would turn 41 directories
into 251, most holding one or two runs, for 9.09 TiB of data that is meant to expire under the
21-day purge. They stay exactly where they are. Nothing in 02_classify reads the flat tree, so
leaving them costs nothing but the disk they were already occupying.

## The move is a rename, and that is the whole reason this is cheap

Source and destination are on the same Lustre mount (verified: both report fsid
8f127ee00000000), so `os.rename` is a metadata operation. 4.05 TiB of cohort reads relocate in
seconds, no bytes are read or written, and **the file's atime is not touched** -- which matters
because atime is what Setonix's 21-day purge clock reads, so a freshly-moved file keeps exactly
the time it had left. That claim is not taken on trust: `--apply` records each file's atime
before the rename and compares it afterwards, and aborts the whole migration on the first
disagreement rather than discovering it 3,000 files later.

The README's "No data moves on Setonix" line meant no re-download and no copy. That is still
true. The line is amended rather than left standing.

## Safety

  * **Dry run is the default.** `--apply` is required to touch anything.
  * **Never overwrites.** A destination that already exists aborts the run. Two files claiming
    one path means the plan is wrong, and picking a winner would hide it.
  * **Verifies every rename**: size, inode and atime must match what the source had. Inode
    equality is the proof it was a rename rather than a copy.
  * **Orphans are reported, not removed**, unless `--clean-orphans` is passed. There are 39
    `_<RUN>_sra` prefetch directories in the flat tree from downloads that died midway; each
    one is evidence a run failed, which is worth seeing before it is deleted.
"""

import argparse
import csv
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _layout import READS, READS_MANIFEST, ROOT, RUN_LIST  # noqa: E402
from _reads import scan_flat, scan_tree, write_manifest  # noqa: E402

# The old flat tree. Named here and nowhere else; after this script has run once it is of
# historical interest only, and 02_classify never reads it.
FLAT = ROOT / "superseded/03_kraken/03_select/data/reads"


def plan() -> tuple[list, dict, list, dict]:
    """(moves, by_project, orphans, skipped) without touching anything."""
    if not RUN_LIST.exists():
        sys.exit(f"no run list: {RUN_LIST}\n"
                 "Run: python 02_classify/03_cohort/resolve.py\n"
                 "(it is tracked, so it should have arrived with the checkout)")
    proj_of = {r["Run"]: r["BioProject"]
               for r in csv.DictReader(open(RUN_LIST), delimiter="\t")}

    present, orphans = scan_flat(FLAT)
    moves, by_project, skipped = [], defaultdict(lambda: [0, 0]), defaultdict(int)
    for run, files in sorted(present.items()):
        proj = proj_of.get(run)
        if not proj:
            # The old STAT-gated set. Left in place on purpose; see the docstring.
            skipped["not in the cohort"] += 1
            continue
        for f in files:
            moves.append((f, READS / proj / f.name))
        by_project[proj][0] += len(files)
        by_project[proj][1] += sum(f.stat().st_size for f in files)
    return moves, dict(by_project), orphans, dict(skipped)


def apply(moves: list) -> None:
    """Rename each file, verifying size, inode and atime across the move.

    Aborts on the FIRST failure. A partial migration is recoverable because every move is an
    independent rename and the manifest is rebuilt by scanning; a migration that carried on
    past a surprise would not be.
    """
    done = 0
    for src, dst in moves:
        if dst.exists():
            sys.exit(f"ABORT: destination already exists: {dst}\n"
                     f"  source: {src}\n"
                     f"Two files claim one path, so the plan is wrong. {done} files were "
                     f"already moved; they are fine, rerun to continue.")
        st = src.stat()
        before = (st.st_size, st.st_ino, st.st_atime)
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.rename(src, dst)
        ds = dst.stat()
        after = (ds.st_size, ds.st_ino, ds.st_atime)
        if after != before:
            sys.exit(f"ABORT after {done} moves: {dst} does not match its source.\n"
                     f"  size/inode/atime before {before}\n"
                     f"  size/inode/atime after  {after}\n"
                     "A changed inode means this was a COPY, not a rename, so the two trees "
                     "are not on one filesystem. A changed atime means the purge clock was "
                     "reset. Stop and investigate before moving anything else.")
        done += 1
        if done % 500 == 0:
            print(f"  {done}/{len(moves)} moved", flush=True)
    print(f"  {done}/{len(moves)} moved")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="actually move the files")
    ap.add_argument("--clean-orphans", action="store_true",
                    help="also delete the `_<RUN>_sra` prefetch debris. Requires --apply")
    args = ap.parse_args()

    if not FLAT.is_dir():
        sys.exit(f"the flat reads tree is not visible from here:\n  {FLAT}\n"
                 "Run this on Setonix.")

    moves, by_project, orphans, skipped = plan()

    print(f"source      {FLAT}")
    print(f"destination {READS}/<BioProject>/\n")
    tb = sum(v[1] for v in by_project.values()) / 2**40
    print(f"{len(moves)} files over {len(by_project)} BioProjects, {tb:.2f} TiB\n")
    print(f"  {'BioProject':14s} {'files':>6s} {'TiB':>7s}")
    for proj in sorted(by_project, key=lambda p: -by_project[p][1]):
        n, b = by_project[proj]
        print(f"  {proj:14s} {n:6d} {b / 2**40:7.3f}")
    for reason, n in sorted(skipped.items()):
        print(f"\n{n} runs left in the flat tree: {reason}")
    if orphans:
        print(f"\n{len(orphans)} orphaned prefetch directories (a download that died midway):")
        for p in orphans[:6]:
            print(f"  {p.name}")
        if len(orphans) > 6:
            print(f"  ... and {len(orphans) - 6} more")
        print("  Each is evidence a run failed. Pass --clean-orphans to remove them.")

    if not args.apply:
        print("\nDRY RUN. Nothing was changed. Add --apply.")
        return
    if not moves:
        print("\nnothing to move")
    else:
        print(f"\nmoving {len(moves)} files...")
        apply(moves)

    if args.clean_orphans and orphans:
        print(f"\nremoving {len(orphans)} orphaned prefetch directories")
        for p in orphans:
            shutil.rmtree(p)

    tree = scan_tree(READS)
    write_manifest(tree, READS_MANIFEST)
    nb = sum(e["bytes"] for e in tree.values())
    print(f"\n-> {READS_MANIFEST.name}: {len(tree)} runs present, {nb / 2**40:.2f} TiB")
    single = sum(1 for e in tree.values() if len(e["files"]) == 1)
    print(f"   {single} runs have ONE file (single-end, or a fasterq-dump dead second mate)")


if __name__ == "__main__":
    main()
