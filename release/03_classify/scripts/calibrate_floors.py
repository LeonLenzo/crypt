#!/usr/bin/env python3
"""Measure the off-host null for one arm, which is what the coverage floor is anchored to.

    python release/03_classify/scripts/calibrate_floors.py --arm db1p
    python release/03_classify/scripts/calibrate_floors.py --arm db2

Reads  data/raw/candidates_<arm>.tsv              every scored species row, from detect.py
       data/inputs/null_biotroph_orders.tsv       which orders are obligate biotrophs
       ../02_reference/data/outputs/supp_data_s3_reference.tsv   order and on_wheat per taxon
       ../01_cohort/data/outputs/supp_data_s2_cohort.tsv         the run's curated host
Writes data/curation/null_offhost_<arm>.tsv       every null row, with its metrics

## What the null is

Obligate biotrophs, scored in *Triticum aestivum* runs, that USDA has never recorded on wheat.
A rust or a powdery mildew of another host cannot be infecting the plant the library came from,
yet reads assigned to it are real: they are cross-mapping from a congener that IS there, plus
whatever was deposited on the leaf. So whatever coverage this population reaches is the ceiling
of the artefact process, measured on real data rather than simulated.

**The floor goes where the null is exhausted, not at a percentile ratio.** Under db_v4 the same
construction put it at 2%, where zero of 10,314 null rows cleared it. A floor that merely sits
"30x the null's 75th percentile" still leaves null rows above it, which is a categorically
weaker sentence.

**This script measures. It does not set the floor.** The decision is written into
`data/inputs/thresholds.tsv` by a person, with this table as its basis. A script that set the
floor from the data it is about to filter would be a floor tuned to a result.

## Three limits of this null, all of which belong in the methods

**It is wheat only.** Restricting to *T. aestivum* runs is what makes the host test exact, and
it drops Triticale on purpose: *x Triticosecale* matches no USDA genus token, so genuine wheat
rusts on it score as off-host and contaminate the tail (max 53.9% rather than 1.97% when that
was last measured).

**It is built from biotrophs**, which have large, well represented genomes. It may not
generalise to small or thinly represented taxa; the axenic controls are what would.

**It bounds cross-mapping, not deposition.** Deposited material is real DNA off a real genome
and spreads k-mers normally, so no abundance or shape statistic separates it from infection.
The null measures how high the artefact process reaches, which is a different question from
whether a high-coverage call is an infection.
"""

import argparse
import collections
import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
RELEASE = MODULE.parent
ORDERS = MODULE / "data" / "inputs" / "null_biotroph_orders.tsv"
THRESHOLDS = MODULE / "data" / "inputs" / "thresholds.tsv"
COHORT = RELEASE / "01_cohort" / "data" / "outputs" / "supp_data_s2_cohort.tsv"
REFERENCE = RELEASE / "02_reference" / "data" / "outputs" / "supp_data_s3_reference.tsv"

ARMS = ("db1p", "db2")
NULL_HOST = "Triticum aestivum"
OUT_COLS = ["run", "BioProject", "host", "arm", "taxid", "taxon", "order", "reads",
            "distinct", "ceiling", "coverage", "redundancy", "saturation", "judgeable"]


def read(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(values: list, q: float) -> float:
    """The q-th percentile by nearest rank, on a sorted list. No numpy dependency."""
    if not values:
        return float("nan")
    i = min(len(values) - 1, max(0, round(q / 100 * (len(values) - 1))))
    return values[i]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", required=True, choices=ARMS)
    ap.add_argument("--candidates", help="default data/raw/candidates_<arm>.tsv")
    args = ap.parse_args()

    arm = args.arm
    cand_path = (Path(args.candidates) if args.candidates
                 else MODULE / "data" / "raw" / f"candidates_{arm}.tsv")
    if not cand_path.is_file():
        sys.exit(f"no candidate table at {cand_path}; run detect.py --arm {arm} first")

    null_orders = {r["order"] for r in read(ORDERS)}
    ref = {int(r["ncbi_taxid"]): r for r in read(REFERENCE) if r["ncbi_taxid"].isdigit()}
    hosts = {r["run"]: r.get("host", "") for r in read(COHORT)}
    floors = {r["parameter"]: float(r["value"]) for r in read(THRESHOLDS)
              if r["arm"] in ("all", arm)}

    rows = []
    for c in read(cand_path):
        if hosts.get(c["run"], "").split("(")[0].strip() != NULL_HOST:
            continue
        rr = ref.get(int(c["taxid"])) if c["taxid"].isdigit() else None
        if not rr or rr.get("order") not in null_orders:
            continue
        if rr.get("on_wheat") in ("yes", "1", "true"):
            continue
        if c["call"] == "competitor_only" or not c["coverage"]:
            continue
        rows.append({"run": c["run"], "BioProject": c["BioProject"],
                     "host": hosts[c["run"]], "arm": arm, "taxid": c["taxid"],
                     "taxon": c["taxon"], "order": rr["order"], "reads": c["reads"],
                     "distinct": c["distinct"], "ceiling": c["ceiling"],
                     "coverage": c["coverage"], "redundancy": c["redundancy"],
                     "saturation": c["saturation"], "judgeable": c["judgeable"]})

    out = MODULE / "data" / "curation" / f"null_offhost_{arm}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    if not rows:
        print(f"-> {out.relative_to(MODULE)}  0 rows. The null is EMPTY, which is not a result: "
              "check that wheat runs are classified and that S3 carries the order column.")
        return

    covs = sorted(float(r["coverage"]) for r in rows)
    # The absolute floor is part of the call, so the null is reported both ways: every row, and
    # only the rows that could have been called at all. The second is the number the coverage
    # floor actually has to clear.
    gated = sorted(float(r["coverage"]) for r in rows
                   if r["distinct"] and int(r["distinct"]) >= floors["distinct_floor"])
    taxa = collections.Counter(r["taxon"] for r in rows)

    print(f"arm {arm}: {len(rows):,} null rows over {len({r['run'] for r in rows}):,} "
          f"wheat runs, {len(taxa)} off-host biotroph species")
    for label, vals in (("all null rows", covs),
                        (f"clearing distinct >= {floors['distinct_floor']:.0f}", gated)):
        if not vals:
            print(f"  {label}: none")
            continue
        print(f"  {label}: n={len(vals):,}  median {statistics.median(vals):.5%}  "
              f"p75 {pct(vals, 75):.5%}  p95 {pct(vals, 95):.5%}  "
              f"p99 {pct(vals, 99):.5%}  max {max(vals):.5%}")
    live = floors["coverage_floor"]
    over = [v for v in gated if v >= live]
    print(f"\n  the floor in thresholds.tsv is {live:.2%}, and {len(over)} of {len(gated):,} "
          f"gated null rows clear it")
    if over:
        print("  THE NULL IS NOT EXHAUSTED AT THIS FLOOR. The exhaustion point is just above "
              f"{max(gated):.5%}; a floor at or above that is the defensible one.")
    else:
        print("  the null is exhausted at this floor, which is the sentence the floor needs")
    print(f"\n  top contributing species: "
          f"{', '.join(f'{t} {n}' for t, n in taxa.most_common(5))}")
    print(f"\n-> {out.relative_to(MODULE)}  {len(rows):,} rows")
    print("This script measures. The floor is decided in data/inputs/thresholds.tsv.")


if __name__ == "__main__":
    main()
