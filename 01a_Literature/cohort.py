#!/usr/bin/env python3
"""THE definition of the cohort. One place, so two answers are impossible.

Added 2026-10-06. The three predicates below were being retyped into every ad-hoc query, and
on that day alone they were written out four times. One copy differed, which produced a
locality count of 425 in one script and 424 in another and cost three tool calls to chase a
discrepancy that meant nothing. A definition used in four places belongs in one.

A run is in the cohort when all three hold:

    setting starts with `field`      the TISSUE was field-grown. A `field*` setting covers
                                    excised-and-forced and transferred-outdoors variants,
                                    which are field samples because the organisms in them
                                    were acquired in the field (leon 2026-10-06).
    library_representative != "no"   one row per biological sample. Only assessed where
                                    duplication was actually found, so most rows are blank
                                    and blank means "representative".
    tissue is not root-like          aerial tissue only. Soil community makes a fungal hit
                                    in root material as likely a saprotroph or endophyte as
                                    a pathogen. Also drops the 3 insect runs.

What is deliberately NOT a criterion: `sampling_selection`. Every stratum is in the cohort;
the flag decides what counts as a finding, not whether the sample counts at all. See
triage.py's sampling_selection comment and crypt-two-evidence-standards.

Run:  python 01a_Literature/cohort.py            # the counts
      python 01a_Literature/cohort.py --write     # also materialise gold/cohort.tsv
"""

from __future__ import annotations

import argparse, collections, csv, re, sys
from pathlib import Path

from _layout import COHORT_TSV, GOLD, RUNS

OUT = COHORT_TSV

NON_AERIAL = re.compile(r"root|rhizo|tuber|nodule|whole.?(plant|seedling)|insect", re.I)


def is_cohort(r: dict) -> bool:
    return (r.get("setting", "").startswith("field")
            and r.get("library_representative", "") != "no"
            and not NON_AERIAL.search(r.get("tissue", "") or ""))


def all_runs(path: Path = RUNS) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cohort(runs: list[dict] | None = None) -> list[dict]:
    return [r for r in (runs if runs is not None else all_runs()) if is_cohort(r)]


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
        with OUT.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(coh[0]), delimiter="\t")
            w.writeheader()
            w.writerows(coh)
        print(f"\nwrote {OUT}  ({len(coh):,} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
