#!/usr/bin/env python3
"""Fetch every run and its BioSample metadata, for the BioProjects the kept papers name.

    python release/01_cohort/scripts/fetch_samples.py

Reads  data/curation/screen_fulltext.tsv   keep papers and their BioProjects
Writes data/raw/samples.tsv                one row per RUN, with its BioSample attributes merged
       data/raw/sample_attrs.jsonl         every BioSample attribute verbatim (gitignored)

ONE ROW PER RUN, not per sample: reads are classified per run downstream, and a few samples
carry several runs (one has 15), so sample-keying would bury run accessions. Each run row
still carries the full attributes of its sample.

Two archive calls per project: the ENA portal read_run report gives every run at once
(run, sample, library, platform, read counts, and the sample fields ENA promotes to columns),
then one ENA BioSample XML per distinct sample gives the submitter's full free-text attributes,
where tissue and locality actually live. The verbatim attributes go to the jsonl so nothing is
lost; the common ones are flattened onto each run row.

library_strategy travels unfiltered: a survey's field RNA-seq and the reference genome it cites
sit in the same report, and only that column tells them apart. The RNA-seq filter is
screen_samples.py's job, not this one's.
"""

import collections
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
DATA = MODULE / "data"
SCREEN = DATA / "curation" / "screen_fulltext.tsv"
OUT = DATA / "raw" / "samples.tsv"
ATTRS = DATA / "raw" / "sample_attrs.jsonl"
PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
PAUSE = 0.25
RUN_FIELDS = ["run_accession", "sample_accession", "experiment_accession", "study_accession",
              "tax_id", "scientific_name", "library_strategy", "library_source",
              "library_selection", "library_layout", "instrument_platform",
              "instrument_model", "read_count", "base_count", "first_public", "country",
              "collection_date"]
# submitter synonyms -> the attribute we want, first match wins
SYN = {"host": ["host", "host scientific name", "specific host"],
       "tissue": ["isolation_source", "isolation source", "tissue", "tissue_type",
                  "organism part", "plant structure", "sample type"],
       "location": ["geo_loc_name", "geographic location", "country", "location"],
       "lat_lon": ["lat_lon", "latitude and longitude"],
       "collection_date": ["collection_date", "collection date", "date"],
       "cultivar": ["cultivar", "variety", "genotype", "breeding line"],
       "dev_stage": ["dev_stage", "developmental stage", "growth stage"],
       "isolate": ["isolate", "strain"]}
COLS = (["run", "BioProject", "sample", "experiment", "tax_id", "scientific_name",
         "library_strategy", "library_source", "library_selection", "library_layout",
         "platform", "model", "read_count", "base_count", "first_public",
         "named_by", "named_in_deposit_section"]
        + list(SYN) + ["bs_country", "bs_collection_date", "n_attributes"])


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def get(url):
    with urllib.request.urlopen(url, timeout=90) as fh:
        return fh.read()


def runs_of(bp):
    url = PORTAL + "?" + urllib.parse.urlencode(
        {"accession": bp, "result": "read_run", "fields": ",".join(RUN_FIELDS),
         "format": "json", "limit": "0"})
    return json.loads(get(url))


def attrs_of(sample):
    xml = get(XML + sample).decode(errors="ignore")
    out = {}
    for tag, val in re.findall(r"<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>", xml, re.S):
        out[tag.strip()] = val.strip()
    return out


def pick(attrs, names):
    low = {k.lower(): v for k, v in attrs.items()}
    return next((low[n] for n in names if low.get(n)), "")


def main():
    screen = read(SCREEN)
    keep = {r["doi"] for r in screen if r["verdict"] == "keep"}
    named_by, deposit = collections.defaultdict(set), collections.defaultdict(bool)
    for r in screen:
        if r["doi"] not in keep:
            continue
        for bp in (r.get("bioprojects") or "").split(";"):
            if bp:
                named_by[bp].add(r["doi"])
        for bp in (r.get("bioprojects_own") or "").split(";"):
            if bp:
                deposit[bp] = True
    projects = sorted(named_by)
    print(f"{len(projects)} BioProjects named by {len(keep)} keep papers")

    cached = {}
    if ATTRS.exists():
        with open(ATTRS, encoding="utf-8") as fh:
            cached = {json.loads(l)["sample"]: json.loads(l)["attrs"]
                      for l in fh if l.strip()}

    runs, samples_seen = [], set()
    for i, bp in enumerate(projects, 1):
        try:
            recs = runs_of(bp)
        except Exception as e:
            print(f"  {bp}: portal {type(e).__name__}")
            recs = []
        for r in recs:
            runs.append((bp, r))
            samples_seen.add(r.get("sample_accession", ""))
        if i % 10 == 0 or i == len(projects):
            print(f"  projects {i}/{len(projects)}, {len(runs)} runs")
        time.sleep(PAUSE)

    todo = [s for s in samples_seen if s and s not in cached]
    print(f"{len(samples_seen)} samples, {len(cached)} cached, {len(todo)} to fetch")
    with open(ATTRS, "a", encoding="utf-8") as jl:
        for i, s in enumerate(todo, 1):
            try:
                a = attrs_of(s)
            except Exception:
                a = {}
            cached[s] = a
            jl.write(json.dumps({"sample": s, "attrs": a}) + "\n")
            if i % 100 == 0 or i == len(todo):
                print(f"  samples {i}/{len(todo)}")
            time.sleep(PAUSE)

    rows = []
    for bp, r in runs:
        s = r.get("sample_accession", "")
        a = cached.get(s, {})
        row = {"run": r.get("run_accession", ""), "BioProject": bp, "sample": s,
               "experiment": r.get("experiment_accession", ""),
               "tax_id": r.get("tax_id", ""), "scientific_name": r.get("scientific_name", ""),
               "library_strategy": r.get("library_strategy", ""),
               "library_source": r.get("library_source", ""),
               "library_selection": r.get("library_selection", ""),
               "library_layout": r.get("library_layout", ""),
               "platform": r.get("instrument_platform", ""),
               "model": r.get("instrument_model", ""),
               "read_count": r.get("read_count", ""), "base_count": r.get("base_count", ""),
               "first_public": r.get("first_public", ""),
               "named_by": ";".join(sorted(named_by[bp])),
               "named_in_deposit_section": int(deposit[bp]),
               "bs_country": r.get("country", ""),
               "bs_collection_date": r.get("collection_date", ""),
               "n_attributes": len(a)}
        for f, names in SYN.items():
            row[f] = pick(a, names)
        rows.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    strat = collections.Counter(r["library_strategy"] or "?" for r in rows)
    rna = [r for r in rows if r["library_strategy"] in ("RNA-Seq", "ncRNA-Seq")]
    print(f"\n-> {OUT.name}  {len(rows)} runs over {len({r['BioProject'] for r in rows})} "
          f"projects, {len({r['sample'] for r in rows})} samples")
    print(f"   library_strategy: {dict(strat.most_common())}")
    print(f"   RNA-seq: {len(rna)} runs over {len({r['BioProject'] for r in rna})} projects")
    for f in ("host", "tissue", "location"):
        n = sum(1 for r in rna if r[f])
        print(f"   RNA-seq {f:<10}{n:>5}/{len(rna)} ({100*n/max(1,len(rna)):.0f}%) from BioSample")


if __name__ == "__main__":
    main()
