#!/usr/bin/env python3
"""Resolve every BioProject the kept papers name into its sequencing runs.

    python release/01_cohort/fetch_runs.py

Reads  curation/accessions.tsv       accessions found in each paper, with type + section
       curation/screen_fulltext.tsv  the verdict, to use only keep papers
Writes raw/runs.tsv                  one row per run, archive metadata

One ENA portal request per BioProject, not per sample: the portal returns every run of a
project in a single read_run report, which is ~150x fewer requests than walking samples.

The per-run library_strategy and library_source TRAVEL rather than being filtered here, so
the RNA-seq / WGS distinction is visible in the data: a rust survey's field RNA-seq and the
reference genome it cites sit in the same report, and only the strategy column tells them
apart. `named_in_deposit_section` carries whether the paper claimed the project as its own
deposit or merely cited it, the own-vs-cited call a later step needs. Nothing is dropped here;
the cohort filter is a separate, readable step.
"""

import collections
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ACC = HERE / "curation" / "accessions.tsv"
SCREEN = HERE / "curation" / "screen_fulltext.tsv"
OUT = HERE / "raw" / "runs.tsv"
PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
PAUSE = 0.3
FIELDS = ["run_accession", "sample_accession", "experiment_accession", "study_accession",
          "tax_id", "scientific_name", "library_strategy", "library_source",
          "library_selection", "library_layout", "instrument_platform", "instrument_model",
          "read_count", "base_count", "first_public", "country", "location",
          "collection_date", "isolate", "host", "cultivar", "tissue_type", "dev_stage"]
COLS = ["run", "BioProject", "sample", "experiment", "tax_id", "scientific_name",
        "library_strategy", "library_source", "library_selection", "library_layout",
        "platform", "model", "read_count", "base_count", "first_public", "country",
        "location", "collection_date", "isolate", "host", "cultivar", "tissue_type",
        "dev_stage", "named_by", "named_in_deposit_section"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def runs_for(bioproject):
    url = PORTAL + "?" + urllib.parse.urlencode(
        {"accession": bioproject, "result": "read_run",
         "fields": ",".join(FIELDS), "format": "json", "limit": "0"})
    with urllib.request.urlopen(url, timeout=90) as fh:
        return json.load(fh)


def main():
    verdict = {r["doi"]: r["verdict"] for r in read_tsv(SCREEN)}
    keep = {d for d, v in verdict.items() if v == "keep"}

    # BioProject -> which keep papers named it, and whether any named it as its own deposit
    named_by = collections.defaultdict(set)
    deposit = collections.defaultdict(bool)
    for a in read_tsv(ACC):
        if a["type"] == "bioproject" and a["doi"] in keep:
            named_by[a["accession"]].add(a["doi"])
            if a["in_deposit_section"] == "1":
                deposit[a["accession"]] = True
    projects = sorted(named_by)
    print(f"{len(projects)} BioProjects named by {len(keep)} keep papers")

    rows, failed = [], []
    for i, bp in enumerate(projects, 1):
        try:
            recs = runs_for(bp)
        except Exception as e:
            failed.append((bp, type(e).__name__))
            recs = []
        for r in recs:
            rows.append({
                "run": r.get("run_accession", ""), "BioProject": bp,
                "sample": r.get("sample_accession", ""),
                "experiment": r.get("experiment_accession", ""),
                "tax_id": r.get("tax_id", ""), "scientific_name": r.get("scientific_name", ""),
                "library_strategy": r.get("library_strategy", ""),
                "library_source": r.get("library_source", ""),
                "library_selection": r.get("library_selection", ""),
                "library_layout": r.get("library_layout", ""),
                "platform": r.get("instrument_platform", ""),
                "model": r.get("instrument_model", ""),
                "read_count": r.get("read_count", ""), "base_count": r.get("base_count", ""),
                "first_public": r.get("first_public", ""), "country": r.get("country", ""),
                "location": r.get("location", ""),
                "collection_date": r.get("collection_date", ""),
                "isolate": r.get("isolate", ""), "host": r.get("host", ""),
                "cultivar": r.get("cultivar", ""), "tissue_type": r.get("tissue_type", ""),
                "dev_stage": r.get("dev_stage", ""),
                "named_by": ";".join(sorted(named_by[bp])),
                "named_in_deposit_section": int(deposit[bp]),
            })
        if i % 10 == 0 or i == len(projects):
            print(f"  {i}/{len(projects)}  {len(rows)} runs so far")
        time.sleep(PAUSE)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    strat = collections.Counter(r["library_strategy"] or "?" for r in rows)
    rna = [r for r in rows if r["library_strategy"] in ("RNA-Seq", "ncRNA-Seq")]
    own = {bp for bp in projects if deposit[bp]}
    rna_own = [r for r in rna if r["named_in_deposit_section"] == 1]
    print(f"\n-> {OUT.name}  {len(rows)} runs over {len({r['BioProject'] for r in rows})} projects")
    print(f"   library_strategy: {dict(strat.most_common())}")
    print(f"   RNA-seq runs: {len(rna)}  over {len({r['BioProject'] for r in rna})} projects")
    print(f"   ... from a project named as the paper's own deposit: {len(rna_own)} "
          f"over {len({r['BioProject'] for r in rna_own})} projects")
    if failed:
        print(f"   FAILED to resolve: {failed}")


if __name__ == "__main__":
    main()
