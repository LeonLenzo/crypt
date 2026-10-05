#!/usr/bin/env python3
"""supp_provenance.py — mine per-sample provenance out of papers' supplementary tables.

The problem this solves: BioSample attributes carry locality for only ~64% of the cohort, and
`llm_geographic_location` is extracted per BIOPROJECT, so a multi-site survey gets one copied
value for every sample. Papers' supplementary tables often carry the real per-sample metadata,
and for surveillance corpora they carry it for hundreds of samples at once.

Worked example: Adams et al. 2021 (BMC Genomics) supplement covers 1,025 samples across 12
studies, and took PRJEB31334 from a single "Ethiopia" label to 75 distinct localities.

Accession flavours are the fiddly part. Supplements cite whatever the journal asked for: ENA
sample accessions (ERS), NCBI run accessions (SRR/ERR) or BioSamples (SAMEA/SAMN). Only the
last is what the pipeline keys on, so everything is resolved through ENA and the map is cached,
because the lookup is one HTTP request per accession.

Adding a supplement is a CONFIG entry, not new code. Header rows are found by looking for a row
containing an accession-like column, so the usual two or three rows of title and section
banners above the real header need no per-file offset.

**Not every supplement describes the sequenced samples, and merging the wrong one is worse than
merging nothing.** `pnas.2601719123.sd07.xlsx` (PRJNA1217477) tabulates the *Botrytis cinerea*
ISOLATES used as inoculum: Isolate, Origin, Country, Host, GPS, Source. Its geography is where
the pathogen cultures were collected, not where the plants were grown, and the plants were all
inoculated at UC Davis. Merging it would stamp Californian vineyard GPS coordinates onto
glasshouse Arabidopsis, and because it carries coordinates it would look like the best
provenance in the dataset while being entirely wrong. Deliberately left unconfigured, so it is
skipped. Before adding a CONFIG entry, check the table keys on a SAMPLE accession and not on an
isolate or strain name.

    python 02_literature/02_text/supp_provenance.py
    python 02_literature/02_text/supp_provenance.py --no-cache   # re-resolve accessions
"""
import argparse
import collections
import csv
import io
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent
SUPP_DIR = HERE / "data/supp_data"
OUT = HERE / "data/supp_provenance.tsv"
CACHE = HERE / "data/supp_accession_map.json"

# Which sheet and columns to read from each supplement. Column names are matched
# case-insensitively after stripping, so trailing spaces in the spreadsheet do not matter.
CONFIG = {
    "12864_2021_7488_MOESM1_ESM.xlsx": {
        "doi": "10.1186/s12864-021-07488-3",
        "sheet": "Sheet1",
        "accession": "Accession Code",
        "country": "Country", "location": "Location",
        "date": "Date Isolated1", "host": "Host Species",
        "sample_type": "Sample Type", "study": "Study",
    },
    "12864_2025_12230_MOESM1_ESM.xlsx": {
        "doi": "10.1186/s12864-025-12230-4",
        "sheet": "Table S1",
        "accession": "Accession",
        "country": "Country", "location": None,   # country and year only, no locality
        "date": "Year", "host": "Host",
        "sample_type": None, "study": "Reference",
    },
    "12915_2019_684_MOESM1_ESM.xlsx": {
        "doi": "10.1186/s12915-019-0684-y",
        "sheet": "Table S1",
        "accession": "Accession number",
        "country": "Country of Origin", "location": "Location",
        "date": "Date Collected", "host": "Host",
        "sample_type": None, "study": "Source",
    },
}

ACC_RE = re.compile(r"^(ERS|SRS|DRS|SAMEA|SAMN|SAMD|ERR|SRR|DRR)\d+$", re.I)
ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run"
       "&fields=run_accession,sample_accession,secondary_sample_accession,study_accession"
       "&format=tsv&limit=0")


def norm(s):
    return re.sub(r"\s+", " ", str(s or "").strip()).lower()


def read_sheet(path, cfg):
    """Rows of one supplement as dicts, with the header row located by content."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[cfg["sheet"]] if cfg["sheet"] in wb.sheetnames else wb[wb.sheetnames[0]]
    rows = [r for r in ws.iter_rows(values_only=True)]
    want = norm(cfg["accession"])
    hdr_i = next((i for i, r in enumerate(rows[:10])
                  if any(norm(c) == want for c in r)), None)
    if hdr_i is None:
        # fall back to any row naming an accession at all
        hdr_i = next((i for i, r in enumerate(rows[:10])
                      if any("accession" in norm(c) for c in r)), 0)
    hdr = [norm(c) for c in rows[hdr_i]]
    out = []
    for r in rows[hdr_i + 1:]:
        if not r or all(c is None for c in r):
            continue
        out.append(dict(zip(hdr, r)))
    return out, hdr


def pick(row, name):
    return "" if not name else str(row.get(norm(name)) or "").strip()


def resolve(accessions, use_cache=True):
    """accession -> BioSample, via ENA, cached on disk."""
    cache = {}
    if use_cache and CACHE.exists():
        cache = json.load(open(CACHE))
    todo = [a for a in accessions if a and a not in cache]
    if todo:
        print(f"  resolving {len(todo)} accessions through ENA "
              f"({len(cache)} already cached)", flush=True)
    for i, a in enumerate(todo, 1):
        try:
            req = urllib.request.Request(ENA.format(a),
                                         headers={"User-Agent": "crypt/supp_provenance"})
            with urllib.request.urlopen(req, timeout=40) as h:
                rr = list(csv.DictReader(io.StringIO(h.read().decode()), delimiter="\t"))
        except Exception:
            rr = []
        for x in rr:
            sam = x.get("sample_accession", "")
            if not sam:
                continue
            for k in (a, x.get("secondary_sample_accession", ""), x.get("run_accession", "")):
                if k:
                    cache[k] = sam
        cache.setdefault(a, "")      # remember the misses too, so a rerun is cheap
        if i % 100 == 0:
            print(f"    {i}/{len(todo)}", flush=True)
        time.sleep(0.12)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    json.dump(cache, open(CACHE, "w"))
    return cache


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-cache", action="store_true", help="re-resolve every accession")
    ap.add_argument("--supp-dir", type=Path, default=SUPP_DIR)
    args = ap.parse_args()

    files = sorted(p for p in args.supp_dir.glob("*.xlsx") if not p.name.startswith("~$"))
    if not files:
        sys.exit(f"no .xlsx found in {args.supp_dir}")

    parsed, unconfigured = [], []
    for p in files:
        cfg = CONFIG.get(p.name)
        if not cfg:
            unconfigured.append(p.name)
            continue
        rows, hdr = read_sheet(p, cfg)
        kept = 0
        for r in rows:
            acc = pick(r, cfg["accession"])
            if not ACC_RE.match(acc):
                continue
            parsed.append({
                "accession": acc, "source_file": p.name, "doi": cfg["doi"],
                "country": pick(r, cfg["country"]), "location": pick(r, cfg["location"]),
                "date": pick(r, cfg["date"])[:10], "host": pick(r, cfg["host"]),
                "sample_type": pick(r, cfg["sample_type"]), "study": pick(r, cfg["study"]),
            })
            kept += 1
        print(f"{p.name}: {kept} rows with a usable accession (of {len(rows)})")
    for n in unconfigured:
        print(f"{n}: SKIPPED, no CONFIG entry (add one naming its sheet and columns)")

    amap = resolve([r["accession"] for r in parsed], not args.no_cache)
    out_rows, unmapped = [], 0
    for r in parsed:
        sam = amap.get(r["accession"], "")
        if not sam:
            unmapped += 1
            continue
        out_rows.append({"BioSample": sam, **r})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["BioSample", "accession", "country", "location", "date", "host",
            "sample_type", "study", "doi", "source_file"]
    with open(OUT, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    print(f"\n{OUT}: {len(out_rows)} rows for {len({r['BioSample'] for r in out_rows})} "
          f"BioSamples ({unmapped} accessions unresolvable)")
    print(f"  countries: {len({r['country'] for r in out_rows if r['country']})}, "
          f"localities: {len({r['location'] for r in out_rows if r['location']})}")
    st = collections.Counter(r["sample_type"] for r in out_rows if r["sample_type"])
    if st:
        print(f"  author-declared sample type: {dict(st)}")


if __name__ == "__main__":
    main()
