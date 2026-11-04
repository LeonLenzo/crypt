#!/usr/bin/env python3
"""Seed curation/supplements.tsv from the two rust-aggregation supplement tables.

    python research/01_literature/retrieve/seed_supplement_rows.py

Reads  release/01_cohort/raw/supplements/*.xlsx   the fetched Additional-file tables
       release/01_cohort/curation/screen_samples.tsv   our runs and their BioSamples
Writes release/01_cohort/curation/supplements.tsv  run -> host/setting/location/date, cited

WHY THIS LIVES IN research/ AND RUNS ONCE. Two wheat yellow-rust papers re-analyse public
SRA and publish a per-accession provenance table (country, year, host, field-vs-lab) that the
archive does not carry. Joining those tables to our runs by accession is deterministic, so a
script does it; but the result is a curation record a human then owns and extends by hand for
the single-project papers, exactly like curated_attributes.tsv was ported once and then kept.
Re-running overwrites only the rows it produced (source begins "supp:"); hand-added rows with
any other source are preserved.

The s12915 table keys on run (SRR/ERR); the s12864 table keys on the old ENA secondary sample
accession (ERS), mapped to our run via one ENA portal call per rust BioProject.
"""

import csv
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
COHORT = REPO / "release" / "01_cohort"
SUPP = COHORT / "raw" / "supplements"
SCREEN = COHORT / "curation" / "screen_samples.tsv"
OUT = COHORT / "curation" / "supplements.tsv"
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
UA = {"User-Agent": "crypt/1.0 (+mailto:leon.lenzo@curtin.edu.au)"}
COLS = ["run", "host", "tissue", "setting", "location", "collection_date", "source"]

# Each spec maps one supplement sheet onto our fields. `key` is the column holding an
# accession we can match to a run (directly, or via ERS). `setting` normalises the paper's
# own Field/Lab label; where the table has none, the body established it (see notes below).
# ORDER MATTERS: the first spec to claim a run wins. s12864 is first because it carries the
# finer provenance (Country + Location + Host + an explicit Field/Lab label); s12915's single
# `Origin` column is coarser and, for isolates both papers list, sometimes wrong (it put
# ERS3526855, a Boshoff 2020 Hopetown field sample, in Zimbabwe), so it only fills what
# s12864 does not reach.
SPECS = [
    {"doi": "10.1186/s12864-021-07488-3", "file": "12864_2021_7488_MOESM1_ESM.xlsx",
     "table": "S1", "key": "Accession Code", "key_kind": "run_or_ers",
     "host": "Host Species", "location": "Country", "date": "Date Isolated1",
     "setting_col": "Sample Type", "setting_fixed": None},
    {"doi": "10.1186/s12915-023-01717-9", "file": "12915_2023_1717_MOESM1_ESM.xlsx",
     "table": "S1", "key": "SRR", "key_kind": "run",
     "host": None, "location": "Origin", "date": "Year of collection", "setting_col": None,
     "setting_fixed": "field"},  # body: field samples, 30 countries, 2013-2018
]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sheet_rows(path):
    """The supplement's first sheet as dicts, finding the header row (some have a title row)."""
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    grid = [[("" if c is None else str(c)).strip() for c in row]
            for row in wb.worksheets[0].iter_rows(values_only=True)]
    # the header is the first row mentioning an accession column we key on
    keys = {"SRR", "Accession Code"}
    h = next(i for i, row in enumerate(grid) if keys & set(row))
    hdr = grid[h]
    return [dict(zip(hdr, row)) for row in grid[h + 1:]]


def run_to_ers(bioprojects):
    """run -> ERS for the given projects, one ENA portal call each."""
    m = {}
    for bp in sorted(bioprojects):
        url = (f"{ENA}?accession={bp}&result=read_run"
               "&fields=run_accession,secondary_sample_accession&format=tsv")
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=90) as fh:
            for r in csv.DictReader(fh.read().decode().splitlines(), delimiter="\t"):
                if r.get("secondary_sample_accession"):
                    m[r["run_accession"]] = r["secondary_sample_accession"]
    return m


def norm_setting(label):
    s = (label or "").strip().lower()
    if s.startswith("field"):
        return "field"
    if s.startswith("lab") or "greenhouse" in s or "glasshouse" in s:
        return "lab"
    return ""


def main():
    runs = {r["run"]: r for r in read_tsv(SCREEN)}
    rust_doi = {s["doi"] for s in SPECS}
    rust_proj = {r["BioProject"] for r in runs.values()
                 if set(r.get("named_by", "").split(";")) & rust_doi}
    ers2run = {ers: run for run, ers in run_to_ers(rust_proj).items()}

    seeded, claimed = {}, 0
    for spec in SPECS:
        for row in sheet_rows(SUPP / spec["file"]):
            acc = row.get(spec["key"], "").strip()
            run = acc if acc in runs else ers2run.get(acc)
            if not run or run not in runs or run in seeded:
                continue
            setting = spec["setting_fixed"] or norm_setting(row.get(spec["setting_col"]))
            seeded[run] = {
                "run": run,
                "host": (row.get(spec["host"]) or "").strip() if spec["host"] else "",
                "tissue": "",  # the tables name no tissue; left for a hand read
                "setting": setting,
                "location": (row.get(spec["location"]) or "").strip(),
                "collection_date": str(row.get(spec["date"]) or "").strip(),
                "source": f"supp:{spec['doi']}#{spec['table']}",
            }
            claimed += 1

    # preserve any hand-added rows (source not beginning "supp:")
    kept = [r for r in (read_tsv(OUT) if OUT.exists() else [])
            if not r.get("source", "").startswith("supp:")]
    rows = list(seeded.values()) + kept
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    import collections
    print(f"-> {OUT.relative_to(REPO)}  {len(rows)} rows ({claimed} from supplements, "
          f"{len(kept)} hand-kept)")
    print(f"   settings: {dict(collections.Counter(r['setting'] for r in seeded.values()))}")
    print(f"   projects: {dict(collections.Counter(runs[r]['BioProject'] for r in seeded))}")


if __name__ == "__main__":
    main()
