#!/usr/bin/env python3
"""Cohort counts, and the materialised gold/cohort.tsv.

The predicate itself lives in `_cohort.py` and is not restated here; see its docstring for
why it may only exist in one place.

Run:  python 01a_Literature/report/cohort.py
      python 01a_Literature/report/cohort.py --write
"""

from __future__ import annotations

import argparse, collections, csv, sys
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _cohort import all_runs, cohort
from _layout import COHORT_TSV, GOLD


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true", help="materialise gold/cohort.tsv")
    args = ap.parse_args()

    runs = all_runs()
    coh = cohort(runs)
    locs = {r["location"] for r in coh if r["location"]}
    dated = sum(1 for r in coh if r["collection_date"])
    print(f"runs.tsv          {len(runs):,}")
    print(f"cohort            {len(coh):,}")
    print(f"  projects        {len({r['BioProject'] for r in coh})}")
    print(f"  localities      {len(locs)}")
    print(f"  dated           {dated:,} ({dated / len(coh):.1%})")
    print(f"  unassessed      {sum(1 for r in runs if not r['setting']):,}")
    print("  setting         " + "  ".join(
        f"{k}:{v}" for k, v in collections.Counter(r["setting"] for r in coh).most_common()))
    print("  selection       " + "  ".join(
        f"{k or '(none)'}:{v}" for k, v in
        collections.Counter(r["sampling_selection"] for r in coh).most_common()))
    if args.write:
        GOLD.mkdir(parents=True, exist_ok=True)
        with COHORT_TSV.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(coh[0]), delimiter="\t")
            w.writeheader()
            w.writerows(coh)
        print(f"\nwrote {COHORT_TSV}  ({len(coh):,} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
