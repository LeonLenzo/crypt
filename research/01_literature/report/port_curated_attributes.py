#!/usr/bin/env python3
"""Port the hand- and supplement-curated run attributes into the release's cohort module.

    python research/01_literature/report/port_curated_attributes.py

Reads  research/01_literature/release/cohort_samples.tsv   the prior curation
       release/01_cohort/raw/runs.tsv                       the new cohort's runs
Writes release/01_cohort/curation/curated_attributes.tsv    one row per overlapping run

This is a RESEARCH-TREE tool, like freeze_cohort_inputs.py: it may read the research tree,
the file it writes may not. The prior curation resolved host, tissue, setting, location and
date for the old cohort by reading the papers and their supplementary tables, and 76% of the
new pathogenomics runs are the same runs. That resolution is human work, not derivable from
the archive, so it enters the release as a curation input rather than being redone. Each value
keeps the id of the rule that set it, so the evidence chain survives the port.

Runs the prior curation never saw are simply absent here; the archive and BioSample attributes
fill those, and whatever remains is the genuinely new curation task.
"""

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
ROOT = MODULE.parent
REPO = ROOT.parent
PRIOR = MODULE / "release" / "cohort_samples.tsv"
RUNS = REPO / "release" / "01_cohort" / "raw" / "runs.tsv"
OUT = REPO / "release" / "01_cohort" / "curation" / "curated_attributes.tsv"

COLS = ["run", "host", "tissue", "setting", "location", "collection_date", "lat", "lon",
        "location_precision", "tissue_rule", "setting_rule", "location_rule",
        "collection_date_rule"]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    if not RUNS.exists():
        sys.exit("run release/01_cohort/fetch_runs.py first")
    newruns = {r["run"] for r in read(RUNS)}
    prior = {r["Run"]: r for r in read(PRIOR)}

    rows = []
    for run in sorted(newruns & set(prior)):
        p = prior[run]
        rows.append({c: p.get(c if c != "run" else "Run", "") for c in COLS})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"-> {OUT.relative_to(REPO)}  {len(rows)} runs ported of {len(newruns)} new runs "
          f"({100 * len(rows) / len(newruns):.0f}%)")
    for c in ("host", "tissue", "setting", "location", "collection_date", "lat"):
        n = sum(1 for r in rows if r[c])
        print(f"   {c:<16}{n:>5}  ({100 * n / max(1, len(rows)):.0f}%)")


if __name__ == "__main__":
    main()
