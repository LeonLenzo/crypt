#!/usr/bin/env python3
"""fetch_run_titles.py — the submitter's own description of each run, from ENA.

    python release/01_cohort/scripts/fetch_run_titles.py
    python release/01_cohort/scripts/fetch_run_titles.py --projects PRJEB12497,PRJNA396589

Reads  data/outputs/supp_data_s2_cohort.tsv   which projects the cohort holds
Writes data/raw/run_titles.tsv                run -> experiment_title, sample_title, library_name

## Why this exists

The BioSample attributes we already cache describe the ISOLATE and the site. They do not say
what the LIBRARY was made from, and for the lab runs that is the only thing that matters: a
germinated-spore library and a 9-dpi infected leaf are different controls and the attribute
record calls both of them `host = wheat`. That is the isolate-vs-sample trap again
(crypt-isolate-vs-sample-trap), and it hid a 16-run pure-fungal control set inside 86 runs
the cohort was carrying as one undifferentiated "lab" bucket.

The run title is the submitter's own words about the library -- "PST - germinating spores
Rep1", "Haustoria-Percoll-3", "RNAseq of Pst 104 E137 A-: 9 dpi" -- so it is evidence about
the material rather than an inference from a field. It is still the ARCHIVE, so it is a
triage signal and the paper is the authority where the two disagree.

ONE REQUEST PER PROJECT, not per run (crypt-ena-portal-fetch). EBI shares a 2 req/s limiter
across its services, so the sleep is not optional politeness.
"""

import argparse
import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
COHORT = HERE / "data/outputs/supp_data_s2_cohort.tsv"
OUT = HERE / "data/raw/run_titles.tsv"
FIELDS = "run_accession,experiment_title,sample_title,library_name,sample_alias"
COLS = ["run", "BioProject", "experiment_title", "sample_title", "library_name",
        "sample_alias"]


def fetch(project, retries=3):
    q = urllib.parse.urlencode({"accession": project, "result": "read_run",
                                "fields": FIELDS, "format": "tsv", "limit": "0"})
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + q
    for n in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return list(csv.DictReader(r.read().decode().splitlines(), delimiter="\t"))
        except Exception as exc:                       # noqa: BLE001
            if n == retries - 1:
                print(f"  {project}: FAILED {exc}", file=sys.stderr)
                return []
            time.sleep(2 ** n)
    return []


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--projects", help="comma-separated; default is every cohort project")
    args = ap.parse_args()

    with open(COHORT) as fh:
        cohort = list(csv.DictReader(fh, delimiter="\t"))
    want = {r["run"]: r["BioProject"] for r in cohort}
    projects = (args.projects.split(",") if args.projects
                else sorted({r["BioProject"] for r in cohort}))

    rows, missing = [], 0
    for p in projects:
        got = fetch(p)
        time.sleep(0.6)
        kept = 0
        for r in got:
            if r["run_accession"] in want:
                rows.append({"run": r["run_accession"], "BioProject": p,
                             **{c: " ".join((r.get(c) or "").split())
                                for c in COLS[2:]}})
                kept += 1
        missing += sum(1 for k, v in want.items() if v == p) - kept
        print(f"  {p}: {kept} of {sum(1 for v in want.values() if v == p)} cohort runs")

    rows.sort(key=lambda r: (r["BioProject"], r["run"]))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    titled = sum(1 for r in rows if r["experiment_title"] or r["sample_title"])
    print(f"-> {OUT}  {len(rows)} runs, {titled} with a title, {missing} not returned")


if __name__ == "__main__":
    main()
