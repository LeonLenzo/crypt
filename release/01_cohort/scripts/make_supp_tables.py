#!/usr/bin/env python3
"""Build the manuscript's cohort supplementary data tables from the curated cohort.

    python release/01_cohort/scripts/make_supp_tables.py

Reads  data/curation/screen_fulltext.tsv   the candidate studies and their BioProjects
       data/curation/screen_samples.tsv     the per-run cohort, each value with its source
Writes data/outputs/supp_data_s1_studies.tsv   one row per candidate study
       data/outputs/supp_data_s2_cohort.tsv     one row per cohort run, values + sources

S1 lists every study the full-text screen assessed, its verdict, the facts read from it, and the
BioProjects it references with how many cohort runs each contributed. S2 is the run-level cohort
exactly as curated: every attribute beside the source that set it (curated / supplement /
biosample / archive), which is what the manuscript cites as the provenance of each value.

This is a projection of the curation, not new curation: it adds nothing a reader could not derive
from the data/curation/ tables, it only arranges them as the paper cites them.
"""

import collections
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
DATA = MODULE / "data"
FULLTEXT = DATA / "curation" / "screen_fulltext.tsv"
SAMPLES = DATA / "curation" / "screen_samples.tsv"
OUTDIR = DATA / "outputs"

S1_COLS = ["doi", "verdict", "host", "field", "tissue", "years", "target_pathogens",
           "n_pathogen_genera", "bioprojects", "bioprojects_own", "n_cohort_runs"]
S2_COLS = ["run", "BioProject", "sample", "own_deposit", "library_strategy", "platform",
           "model", "read_count", "host", "host_src", "tissue", "tissue_src",
           "setting", "setting_src", "location", "location_src",
           "collection_date", "collection_date_src", "lat", "lon"]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    cohort = read(SAMPLES)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # runs per BioProject, to report each study's contribution to the cohort
    runs_by_proj = collections.Counter(r["BioProject"] for r in cohort)

    # ---- S1: candidate studies ----
    s1 = []
    for r in read(FULLTEXT):
        projs = [p for p in (r.get("bioprojects") or "").split(";") if p]
        s1.append({
            "doi": r["doi"], "verdict": r.get("verdict", ""),
            "host": r.get("host", ""), "field": r.get("field", ""),
            "tissue": r.get("tissue", ""), "years": r.get("years", ""),
            "target_pathogens": r.get("target_pathogens", ""),
            "n_pathogen_genera": r.get("n_pathogen_genera_named", ""),
            "bioprojects": r.get("bioprojects", ""),
            "bioprojects_own": r.get("bioprojects_own", ""),
            "n_cohort_runs": sum(runs_by_proj.get(p, 0) for p in projs),
        })
    with open(OUTDIR / "supp_data_s1_studies.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=S1_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(s1)

    # ---- S2: per-run cohort ----
    with open(OUTDIR / "supp_data_s2_cohort.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=S2_COLS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(cohort)

    kept = sum(1 for r in s1 if r["verdict"] == "keep")
    print(f"-> data/outputs/supp_data_s1_studies.tsv   {len(s1)} candidate studies ({kept} kept)")
    print(f"-> data/outputs/supp_data_s2_cohort.tsv     {len(cohort)} runs, "
          f"{len(runs_by_proj)} BioProjects, {len({r['sample'] for r in cohort})} samples")


if __name__ == "__main__":
    main()
