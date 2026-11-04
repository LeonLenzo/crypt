#!/usr/bin/env python3
"""Assemble the run-level cohort, filling each attribute from the best available source.

    python release/01_cohort/build_cohort.py

Reads  raw/samples.tsv                   every run + its BioSample metadata, one table
       curation/curated_attributes.tsv   prior hand/supplement curation, ported
       curation/screen_fulltext.tsv      the keep/drop verdict per paper
Writes curation/screen_samples.tsv       one row per RNA-seq run in a kept project

PRECEDENCE, highest first, recorded per value in a `<field>_src` column so nothing is taken
on trust: curated (a human read the paper or its supplement) > biosample (the submitter's
attribute record) > archive (the run report's fixed field). A value and where it came from
travel together, because "leaf, from the paper" and "leaf, guessed from a sparse archive
field" are not the same evidence and a reviewer must be able to tell them apart.

Scope: RNA-seq / ncRNA-Seq runs only (WGS reference genomes and AMPLICON assays are in
runs.tsv, labelled, and excluded here), from BioProjects named by a paper the screen kept.
Nothing about the foliar or depth gates happens here; those belong to 02_classify, which
this feeds.
"""

import collections
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
SAMPLES = HERE / "raw" / "samples.tsv"
CURATED = HERE / "curation" / "curated_attributes.tsv"
SCREEN = HERE / "curation" / "screen_fulltext.tsv"
OUT = HERE / "curation" / "screen_samples.tsv"

RNA = {"RNA-Seq", "ncRNA-Seq"}
FIELDS = ["host", "tissue", "setting", "location", "collection_date", "lat", "lon"]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    keep = {r["doi"] for r in read(SCREEN) if r["verdict"] == "keep"}
    curated = {r["run"]: r for r in read(CURATED)}
    # samples.tsv carries BioSample-derived host/tissue/location/date and the archive
    # country/date on each run row; `setting` has no archive source (only the paper says it).
    arch_map = {"host": "host", "tissue": "tissue", "setting": None,
                "location": "location", "collection_date": "collection_date",
                "lat": None, "lon": None}

    rows, src_tally = [], collections.Counter()
    for r in read(SAMPLES):
        if r["library_strategy"] not in RNA:
            continue
        # a run belongs if any kept paper named its project; runs.tsv carries named_by
        if not (set(r.get("named_by", "").split(";")) & keep):
            continue
        cur = curated.get(r["run"], {})
        bs = r  # samples.tsv row already holds this sample's attributes
        out = {"run": r["run"], "BioProject": r["BioProject"], "sample": r["sample"],
               "library_strategy": r["library_strategy"], "platform": r["platform"],
               "model": r["model"], "read_count": r["read_count"],
               "named_by": r["named_by"],
               "own_deposit": r["named_in_deposit_section"]}
        for f in FIELDS:
            val, src = "", ""
            if cur.get(f):
                val, src = cur[f], "curated"
            elif bs.get(f):
                val, src = bs[f], "biosample"
            elif arch_map.get(f) and r.get(arch_map[f]):
                val, src = r[arch_map[f]], "archive"
            out[f] = val
            out[f"{f}_src"] = src
            src_tally[(f, src or "none")] += 1
        rows.append(out)

    cols = (["run", "BioProject", "sample", "library_strategy", "platform", "model",
             "read_count", "named_by", "own_deposit"]
            + [c for f in FIELDS for c in (f, f"{f}_src")])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"-> {OUT.name}  {len(rows)} RNA-seq runs over "
          f"{len({r['BioProject'] for r in rows})} BioProjects, "
          f"{len({r['sample'] for r in rows})} BioSamples")
    print(f"   from a paper's own deposit: {sum(1 for r in rows if r['own_deposit']=='1')}")
    print("\n   attribute coverage, and from where:")
    for f in FIELDS:
        c = {s: src_tally[(f, s)] for s in ("curated", "biosample", "archive", "none")}
        filled = len(rows) - c["none"]
        print(f"     {f:<16}{filled:>5}/{len(rows)} ({100*filled/len(rows):>3.0f}%)   "
              f"curated {c['curated']}, biosample {c['biosample']}, archive {c['archive']}")


if __name__ == "__main__":
    main()
