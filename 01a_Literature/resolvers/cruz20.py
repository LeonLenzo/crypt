#!/usr/bin/env python3
"""Pick one representative library per plant for the Zwijnaarde single-plant maize field trial.

    Cruz DF et al. 2020. "Using single-plant-omics in the field to link maize genes to
    functions and phenotypes." Molecular Systems Biology 16:e9667.
    doi:10.15252/msb.20209667   PMC7751767   open access.   PRJEB37824

This is the cohort's first genuinely field cereal project of any size, and it is unusual in a
way that matters: the samples are INDIVIDUAL PLANTS, not plot pools. Every other project
curated on 2026-10-06/07 pooled - three plants per leaf replicate, up to 32 ears per plot,
five heads per pot - so a co-detection there is per-plot. Here it is per-plant.

    "During the summer of 2015, 560 B104 maize inbred plants were grown under 'uncontrolled'
     field conditions at a site in Zwijnaarde, Belgium (51 00'35.2\"N, 3 42'56.5\"E)"
    "In total, 200 non-border plants that exhibited a primary ear at leaf 16 were harvested at
     the VT (tasseling) stage... plants were harvested on two different dates, 2015-08-25
     (164 plants) and 2015-09-02 (36 plants)"
    "Sixty of the 200 leaf samples for individual plants were randomly selected for
     RNA-sequencing"

## Why a resolver rather than a predicate

SRA holds 92 runs for those 60 plants, because "Sequencing was performed in paired-end mode
... in two batches". The plant is recoverable from `LibraryName`, which is the only place it
lives: `plant175_4393_3_p`. Counting them gives 29 plants with one library, 30 with two and
one with three, which is 92 and confirms the reading.

Without a representative flag the 30 doubly-sequenced plants would be counted twice, which is
the same error the Sato cluster produced at much larger scale. The choice of which library
keeps the flag is DETERMINISTIC and shallow: the deepest run per plant, with the accession as
a tie-break. The paper gives no basis for preferring one batch over the other, and depth is
what a quantification would prefer; recording the rule matters more than the rule itself.

Run:  python 01a_Literature/resolvers/cruz20.py
"""

from __future__ import annotations

import collections, csv, re, sys
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _resolver import finish
from _layout import RUNS, STUDIES

BIOPROJECT = "PRJEB37824"
OUT = STUDIES / "doi_10.15252_msb.20209667" / "cruz20_runs.csv"

N_PLANTS = 60        # "Sixty of the 200 leaf samples ... were randomly selected"
PLANT = re.compile(r"^plant(\d+)_")


def main() -> None:
    with open(RUNS) as fh:
        ours = [r for r in csv.DictReader(fh, delimiter="\t") if r["BioProject"] == BIOPROJECT]
    if not ours:
        sys.exit(f"REFUSED: no runs for {BIOPROJECT}")

    by_plant = collections.defaultdict(list)
    for r in ours:
        m = PLANT.match(r["LibraryName"] or "")
        if not m:
            sys.exit(f"REFUSED: {r['Run']} has LibraryName {r['LibraryName']!r}, which carries "
                     f"no plant number. The plant is the only thing that makes these runs "
                     f"dedupable and it lives nowhere else.")
        by_plant[m.group(1)].append(r)

    if len(by_plant) != N_PLANTS:
        sys.exit(f"REFUSED: {len(by_plant)} distinct plants, but the paper sequenced "
                 f"{N_PLANTS}. Either the parse is wrong or this is not that experiment.")

    rows = []
    for plant, runs in sorted(by_plant.items(), key=lambda kv: int(kv[0])):
        # Deepest first, accession as the tie-break so the choice is reproducible.
        runs.sort(key=lambda r: (-int(r["spots"] or 0), r["Run"]))
        for i, r in enumerate(runs):
            rows.append(dict(
                Run=r["Run"], Plant=f"plant{plant}", LibraryName=r["LibraryName"],
                Spots=r["spots"], NLibraries=len(runs),
                LibraryRepresentative="yes" if i == 0 else "no"))

    keep = sum(1 for r in rows if r["LibraryRepresentative"] == "yes")
    if keep != N_PLANTS:
        sys.exit(f"REFUSED: flagged {keep} representatives for {N_PLANTS} plants")

    print(f"{len(rows)} runs across {len(by_plant)} plants", file=sys.stderr)
    print(f"  libraries per plant: "
          f"{dict(sorted(collections.Counter(len(v) for v in by_plant.values()).items()))}",
          file=sys.stderr)
    print(f"  representatives: {keep}, dropped as duplicates: {len(rows) - keep}",
          file=sys.stderr)
    finish(rows, OUT, key_col="Run", bioprojects=[BIOPROJECT])


if __name__ == "__main__":
    main()
