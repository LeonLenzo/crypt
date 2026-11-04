#!/usr/bin/env python3
"""Assemble the run-level cohort, filling each attribute from the best available source.

    python release/01_cohort/scripts/build_cohort.py

Reads  data/raw/samples.tsv                   every run + its BioSample metadata, one table
       data/curation/curated_attributes.tsv   prior hand/supplement curation, ported
       data/curation/supplements.tsv          per-accession provenance lifted from paper supplements
       data/curation/screen_fulltext.tsv      the keep/drop verdict per paper
Writes data/curation/screen_samples.tsv       one row per RNA-seq run in a kept project

PRECEDENCE, highest first, recorded per value in a `<field>_src` column so nothing is taken
on trust: curated > supplement > biosample (the submitter's attribute record) > archive (the
run report's fixed field). `curated` and `supplement` are both a human reading the paper;
curated ranks first only because it is the established port and the two agree wherever they
overlap, so supplement serves to fill the runs curation never reached. A value and where it
came from travel together, because "leaf, from the paper" and "leaf, guessed from a sparse
archive field" are not the same evidence and a reviewer must be able to tell them apart.

Scope: RNA-seq / ncRNA-Seq runs only (WGS reference genomes and AMPLICON assays are in
runs.tsv, labelled, and excluded here), from BioProjects named by a paper the screen kept.
Nothing about the foliar or depth gates happens here; those belong to 02_classify, which
this feeds.
"""

import collections
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
DATA = MODULE / "data"
SAMPLES = DATA / "raw" / "samples.tsv"
CURATED = DATA / "curation" / "curated_attributes.tsv"
SUPPLEMENTS = DATA / "curation" / "supplements.tsv"
SCREEN = DATA / "curation" / "screen_fulltext.tsv"
OUT = DATA / "curation" / "screen_samples.tsv"

RNA = {"RNA-Seq", "ncRNA-Seq"}
FIELDS = ["host", "tissue", "setting", "location", "collection_date", "lat", "lon"]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def normalize_host(h):
    """Fold the free-text host a submitter typed to a clean binomial.

    BioSample host is uncontrolled text: "wheat", "Soft red winter wheat: IL11-28222",
    "Triticum aestivum (bread wheat)", "not collected", "PDA" (a culture medium, not a host).
    The cohort is a few cereals, so the mapping is small and explicit rather than a taxonomy
    lookup; the cultivar after a ':'/';' is dropped (S2 does not carry it) and non-hosts blank.
    """
    h = (h or "").strip()
    low = h.lower()
    if low in ("", "none", "not collected", "missing", "na", "n/a", "pda"):
        return ""
    if low == "wheat" or ("soft re" in low and "winter wheat" in low):
        return "Triticum aestivum"
    if low == "perenial grass":
        return "perennial grass"
    if low.startswith("triticale"):
        return "Triticosecale"
    if low.startswith("triticum durum"):
        return "Triticum turgidum"
    h = re.sub(r"\s*\(.*?\)", "", h)   # drop a common-name gloss: "(bread wheat)", "(rye)"
    return re.split(r"[;:]", h)[0].strip()


def split_latlon(s):
    """Decimal (lat, lon) from a BioSample lat_lon string, e.g. '42.5 N 71.1 W' or '42.5,-71.1'."""
    hemi = re.findall(r"(-?\d+\.\d+)\s*([NSEW])", s or "")
    if len(hemi) == 2:
        d = {c: float(v) for v, c in hemi}
        lat = d["N"] if "N" in d else (-d["S"] if "S" in d else None)
        lon = d["E"] if "E" in d else (-d["W"] if "W" in d else None)
        return lat, lon
    nums = re.findall(r"-?\d+\.\d+", s or "")
    return (float(nums[0]), float(nums[1])) if len(nums) == 2 else (None, None)


def main():
    keep = {r["doi"] for r in read(SCREEN) if r["verdict"] == "keep"}
    curated = {r["run"]: r for r in read(CURATED)}
    supp = {r["run"]: r for r in read(SUPPLEMENTS)} if SUPPLEMENTS.exists() else {}
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
        sup = supp.get(r["run"], {})
        bs = dict(r)  # samples.tsv row already holds this sample's attributes
        # the archive stores coordinates as one combined lat_lon field; split it so the
        # biosample tier can fill lat/lon like any other attribute
        la, lo = split_latlon(r.get("lat_lon", ""))
        bs["lat"] = f"{la:.5f}" if la is not None else ""
        bs["lon"] = f"{lo:.5f}" if lo is not None else ""
        out = {"run": r["run"], "BioProject": r["BioProject"], "sample": r["sample"],
               "library_strategy": r["library_strategy"], "platform": r["platform"],
               "model": r["model"], "read_count": r["read_count"],
               "named_by": r["named_by"],
               "own_deposit": r["named_in_deposit_section"]}
        for f in FIELDS:
            val, src = "", ""
            if cur.get(f):
                val, src = cur[f], "curated"
            elif sup.get(f):
                val, src = sup[f], "supplement"
            elif bs.get(f):
                val, src = bs[f], "biosample"
            elif arch_map.get(f) and r.get(arch_map[f]):
                val, src = r[arch_map[f]], "archive"
            if f == "host":
                val = normalize_host(val)
                if not val:          # a non-host like "PDA" or "not collected": no value, no source
                    src = ""
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
        c = {s: src_tally[(f, s)]
             for s in ("curated", "supplement", "biosample", "archive", "none")}
        filled = len(rows) - c["none"]
        print(f"     {f:<16}{filled:>5}/{len(rows)} ({100*filled/len(rows):>3.0f}%)   "
              f"curated {c['curated']}, supplement {c['supplement']}, "
              f"biosample {c['biosample']}, archive {c['archive']}")


if __name__ == "__main__":
    main()
