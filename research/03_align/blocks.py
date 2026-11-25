#!/usr/bin/env python3
"""The single-versus-multi block inventory for question 3, computed rather than assumed.

    python research/03_align/blocks.py

Reads  release/03_classify/data/outputs/supp_data_s5_detections.tsv   the db2 calls
       release/03_classify/data/outputs/supp_data_s6_runs.tsv         one row per run
       release/01_cohort/data/outputs/supp_data_s2_cohort.tsv         lat / lon per run
Writes nothing. Prints the table the README quotes.

The design constraint is that site is confounded with cultivar, sowing date, protocol and
depth, so the contrast has to be WITHIN a site and site is a blocking factor rather than a
result. This script answers the only question that decides whether that design exists: how
many sites hold enough of both classes.

`single` is exactly one call with is_wheat_pathogen = yes, `multi` is two or more, and runs
with none are counted and then set aside: a run with no call carries a resident community
rather than nothing, which is a third state and not a zero. Resident taxa are deliberately
not in the definition.

COORDINATES COME FROM THE COHORT TABLE, NOT FROM A FIGURE'S COLUMN TABLE. The samples figure
fills a drawing position for every run it draws; a run whose only location is a country name
would then sit at that country's centroid and invent a block shared with every other run from
that country. Only runs with lat / lon in supp_data_s2 are eligible, which is 1,062 of 1,379.

The grid sweep is the check that the inventory is not an artefact of where the rounding falls:
if 100 m, 1 km and 10 km disagree, the blocks are being made by the rounding.
"""

import collections
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLS = ROOT / "release/03_classify/data/outputs"
COH = ROOT / "release/01_cohort/data/outputs/supp_data_s2_cohort.tsv"
ARM = "db2"
MIN_PER_CLASS = 2


def rd(path):
    with open(path, newline="", encoding="utf-8") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main():
    wheat_pathogen_calls = collections.Counter()
    for r in rd(CLS / "supp_data_s5_detections.tsv"):
        if r["arm"] == ARM and r["is_wheat_pathogen"] == "yes":
            wheat_pathogen_calls[r["run"]] += 1

    runs = [r for r in rd(CLS / "supp_data_s6_runs.tsv")
            if r["host"] == "Triticum aestivum" and r["setting"] == "field"]
    coh = {r["run"]: r for r in rd(COH)}

    def klass(run):
        n = wheat_pathogen_calls[run]
        return "none" if n == 0 else ("single" if n == 1 else "multi")

    counts = collections.Counter(klass(r["run"]) for r in runs)
    placed = [r for r in runs if coh.get(r["run"], {}).get("lat")]

    print(f"field T. aestivum runs                     {len(runs):>6}")
    for k in ("single", "multi", "none"):
        print(f"  {k:<40} {counts[k]:>6}")
    print(f"runs carrying real coordinates             {len(placed):>6}")

    def blocks(precision):
        out = collections.defaultdict(collections.Counter)
        for r in placed:
            c = coh[r["run"]]
            key = (r["BioProject"], round(float(c["lat"]), precision),
                   round(float(c["lon"]), precision))
            out[key][klass(r["run"])] += 1
        return out

    for precision, label in ((3, "~100 m"), (2, "~1 km"), (1, "~10 km")):
        b = blocks(precision)
        informative = {k: v for k, v in b.items()
                       if v["single"] >= MIN_PER_CLASS and v["multi"] >= MIN_PER_CLASS}
        sizes = sorted(v["single"] + v["multi"] for v in informative.values())
        median = sizes[len(sizes) // 2] if sizes else 0
        print(f"grid {label:<7} {len(b):>4} blocks, {len(informative):>3} informative, "
              f"{sum(sizes):>4} runs in them, median {median}")

    b = blocks(2)
    per = collections.Counter()
    per_runs = collections.Counter()
    for (project, _, _), v in b.items():
        if v["single"] >= MIN_PER_CLASS and v["multi"] >= MIN_PER_CLASS:
            per[project] += 1
            per_runs[project] += v["single"] + v["multi"]
    print("\ninformative blocks by project, at ~1 km:")
    for project, n in per.most_common():
        note = "  INOCULATED, the manipulated check, out of the estimate" \
            if project == "PRJNA718488" else ""
        print(f"  {project:<12} {n:>2} blocks  {per_runs[project]:>4} runs{note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
