#!/usr/bin/env python3
"""Fetch the full BioSample attribute record for every sample in raw/runs.tsv.

    python release/01_cohort/fetch_biosamples.py

Reads  raw/runs.tsv          the resolved runs (uses their distinct sample accessions)
Writes raw/biosamples.tsv    one row per BioSample, a curated set of common columns
       raw/biosample_attrs.jsonl   one JSON object per BioSample: EVERY attribute, verbatim

WHY BOTH. The run report promotes only a fixed set of sample fields to columns, and they are
thinly populated: tissue reached 10% of runs, location 4%. The BioSample XML carries the full
attribute set the submitter actually provided, and submitters name the same thing many ways
(`isolation_source`, `tissue`, `plant structure`, `organism part`). So the verbatim attributes
go to a JSONL where nothing is lost, and a best-effort flattening of the common ones goes to a
TSV for quick use. Curation reads the JSONL when the TSV's guess is wrong.

One ENA browser XML request per BioSample. Re-running skips samples already in the JSONL.
"""

import collections
import csv
import json
import re
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "raw" / "runs.tsv"
OUT_TSV = HERE / "raw" / "biosamples.tsv"
OUT_JSONL = HERE / "raw" / "biosample_attrs.jsonl"
XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
PAUSE = 0.25

# submitter synonyms -> the column we want. First match wins, checked in order.
MAP = {
    "host": ["host", "host scientific name", "specific host"],
    "tissue": ["isolation_source", "isolation source", "tissue", "tissue_type",
               "organism part", "plant structure", "plant_structure", "sample type"],
    "location": ["geo_loc_name", "geographic location", "geo loc name", "country",
                 "collection location", "location"],
    "lat_lon": ["lat_lon", "latitude and longitude", "lat lon"],
    "collection_date": ["collection_date", "collection date", "date", "sampling date"],
    "cultivar": ["cultivar", "variety", "genotype", "breeding line", "host cultivar"],
    "dev_stage": ["dev_stage", "developmental stage", "growth stage", "plant growth stage",
                  "phenological stage"],
    "isolate": ["isolate", "strain", "pathogen isolate"],
}
COLS = ["sample", "organism", "host", "tissue", "location", "lat_lon", "collection_date",
        "cultivar", "dev_stage", "isolate", "n_attributes"]


def attrs_of(sample):
    with urllib.request.urlopen(XML + sample, timeout=60) as fh:
        xml = fh.read().decode(errors="ignore")
    out = {}
    for tag, val in re.findall(r"<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>", xml, re.S):
        out[tag.strip()] = val.strip()
    org = re.search(r"<SCIENTIFIC_NAME>(.*?)</SCIENTIFIC_NAME>", xml)
    if org:
        out.setdefault("organism", org.group(1).strip())
    return out


def pick(attrs, names):
    low = {k.lower(): v for k, v in attrs.items()}
    for n in names:
        if low.get(n):
            return low[n]
    return ""


def main():
    with open(RUNS, newline="", encoding="utf-8") as fh:
        samples = sorted({r["sample"] for r in csv.DictReader(fh, delimiter="\t")
                          if r["sample"]})
    done = set()
    if OUT_JSONL.exists():
        with open(OUT_JSONL, encoding="utf-8") as fh:
            done = {json.loads(l)["sample"] for l in fh if l.strip()}
    todo = [s for s in samples if s not in done]
    print(f"{len(samples)} BioSamples, {len(done)} cached, {len(todo)} to fetch")

    failed = []
    with open(OUT_JSONL, "a", encoding="utf-8") as jl:
        for i, s in enumerate(todo, 1):
            try:
                a = attrs_of(s)
            except Exception as e:
                failed.append((s, type(e).__name__))
                a = {}
            jl.write(json.dumps({"sample": s, "attrs": a}) + "\n")
            if i % 100 == 0 or i == len(todo):
                print(f"  {i}/{len(todo)}")
            time.sleep(PAUSE)

    # flatten the JSONL (all of it, cached + new) to the TSV
    rows = []
    with open(OUT_JSONL, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            a = rec["attrs"]
            rows.append({"sample": rec["sample"], "organism": a.get("organism", ""),
                         "host": pick(a, MAP["host"]), "tissue": pick(a, MAP["tissue"]),
                         "location": pick(a, MAP["location"]), "lat_lon": pick(a, MAP["lat_lon"]),
                         "collection_date": pick(a, MAP["collection_date"]),
                         "cultivar": pick(a, MAP["cultivar"]),
                         "dev_stage": pick(a, MAP["dev_stage"]),
                         "isolate": pick(a, MAP["isolate"]), "n_attributes": len(a)})
    rows.sort(key=lambda x: x["sample"])
    with open(OUT_TSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\n-> {OUT_TSV.name}  {len(rows)} BioSamples")
    for c in ("host", "tissue", "location", "collection_date", "cultivar", "dev_stage"):
        n = sum(1 for r in rows if r[c])
        print(f"   {c:<16}{n:>5}  ({100 * n / max(1, len(rows)):.0f}%)")
    if failed:
        print(f"   FAILED: {len(failed)} e.g. {failed[:3]}")


if __name__ == "__main__":
    main()
