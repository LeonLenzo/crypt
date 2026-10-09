#!/usr/bin/env python3
"""Cohort counts, and the materialised gold/cohort.tsv.

The predicate itself lives in `_cohort.py` and is not restated here; see its docstring for
why it may only exist in one place.

Run:  python 01a_Literature/report/cohort.py
      python 01a_Literature/report/cohort.py --write
"""

from __future__ import annotations

import argparse, collections, csv, re, sys
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
    # Crop x selection, because it is the structural fact the chapter turns on and no other
    # report shows it. Measured 2026-10-09: wheat is 96% disease-selected, rice and sorghum
    # are 100% unselected. The strata do not overlap, so a crop comparison is not estimable
    # and the two evidence standards are the only honest framing. Printed as a cross-tab
    # rather than two separate distributions precisely so the empty cells are visible.
    hosts = GOLD / "cohort_hosts.tsv"
    if hosts.exists():
        crop_of = {}
        for r in csv.DictReader(hosts.open(), delimiter="\t"):
            for pat, name in ((r"^Triticum", "Wheat"), (r"^Zea mays($| subsp\.? mays| \()", "Maize"),
                              (r"^Oryza sativa", "Rice"), (r"^Sorghum bicolor", "Sorghum")):
                if re.match(pat, r["host"] or ""):
                    crop_of[r["Run"]] = name
                    break
        tab = collections.defaultdict(collections.Counter)
        for r in coh:
            c = crop_of.get(r["Run"])
            if c:
                tab[c][r["sampling_selection"] or "(none)"] += 1
        sels = [s for s in ("unselected", "disease-selected", "inoculated", "symptom-avoided",
                            "fungicide-treated", "(none)")
                if any(t[s] for t in tab.values())]
        if tab:
            print("\n  crop x sampling_selection   (cereal field cohort; empty cells are the point)")
            print("    " + f"{'':<9}" + "".join(f"{s[:17]:>19}" for s in sels))
            for c in ("Wheat", "Maize", "Rice", "Sorghum"):
                if not tab[c]:
                    continue
                tot = sum(tab[c].values())
                print("    " + f"{c:<9}" + "".join(
                    f"{tab[c][s]:>12} ({100*tab[c][s]/tot:>3.0f}%)" for s in sels))

    if args.write:
        GOLD.mkdir(parents=True, exist_ok=True)
        with COHORT_TSV.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(coh[0]), delimiter="\t")
            w.writeheader()
            w.writerows(coh)
        print(f"\nwrote {COHORT_TSV}  ({len(coh):,} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
