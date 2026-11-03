#!/usr/bin/env python3
"""Per-isolate collection locality and date for PRJEB47693, from Lewis et al. 2022 Table S1.

    Lewis CM, Persoons A, Bebber DP, et al. 2022. "Wheat stem rust recorded for the first
    time in decades in Ireland ..." Plant Pathology 71:890.
    doi:10.1111/ppa.13532, Table S1 (PPA-71-890-s005.xlsx), pulled from Europe PMC
    PMC9303354 on 2026-10-12.

23 of the project's 34 runs sat on country centroids: ENA carries `country` and nothing
finer. Table S1 is the isolate table the article points at and it has the collection site,
down to the Irish county and the French commune.

## The key

Table S1's `Accession` is the ENA SECONDARY sample accession (`ERS7651702`), the same form
Tsushima 2025 uses. `sources.py ena` gained `secondary_sample_accession` on 2026-10-12;
refetch any project cached before that date or this resolver cannot run.

## Shape of the sheet

One row per accession, but the isolate-level columns are filled ONLY on the first row of
each isolate's block and left blank on its continuation rows. They must be forward-filled or
most rows read as having no country at all. The header is on row 2, under a title line.

`Latitude and longitude` is a column but is `-` on every row, as in this group's other
supplements; the place name is all there is, and the shared geocoder handles it.

## What it yields

Of the 34 runs, 7 are German genomic DNA libraries and 27 are the transcriptomes. 25 of the
27 gain a locality; `IT-05` and `FR-01` give `-` in Table S1 itself and keep their country.
`Louville la Chenard` is Louville-la-Chenard in Eure-et-Loir, already in the cohort as
"France: Louville" from PRJEB39201, the same surveillance network.

Collection dates come too, and they are finer than what we hold: Table S1 gives the day
(2020-07-28) where the archive gives the year.

    python 01_literature/resolvers/lewis22.py
"""

from __future__ import annotations

import collections, csv, re, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import ENA, STUDIES
from _resolver import finish

BIOPROJECT = "PRJEB47693"
STUDY = STUDIES / "doi_10.1111_ppa.13532"
XLSX = STUDY / "PPA-71-890-s005.xlsx"
SHEET = "Table S1"
OUT = STUDY / "lewis22_runs.csv"

N_EXPECTED = 34          # every run in the project appears in Table S1
N_WITH_LOCALITY = 32     # IT-05 and FR-01 give "-"

# Blank on a continuation row means "same as the isolate above", not "unknown".
FILL = ["Isolate ID", "Original isolate ID", "Date collected", "Country", "Location",
        "Latitude and longitude", "Host", "Host variety"]


def iso_date(v: str) -> str:
    """Table S1 mixes three date forms: 2019, 16/07/2013, and 2020-06-30 00:00:00."""
    v = v.strip()
    if m := re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})(?: .*)?", v):
        return f"{m.group(1)}-{m.group(2)}-{m.group(3)}"
    if m := re.fullmatch(r"(\d{2})/(\d{2})/(\d{4})", v):
        return f"{m.group(3)}-{m.group(2)}-{m.group(1)}"
    if re.fullmatch(r"\d{4}", v):
        return v
    return ""


def main() -> None:
    import openpyxl
    if not XLSX.exists():
        sys.exit(f"missing {XLSX}; Table S1 of doi:10.1111/ppa.13532 "
                 f"(Europe PMC {('PMC9303354')} supplementaryFiles)")

    rows = [["" if c is None else str(c).strip() for c in r]
            for r in openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
            [SHEET].iter_rows(values_only=True)]
    hdr = next((r for r in rows if "Location" in r and "Accession" in r), None)
    if hdr is None:
        sys.exit(f"REFUSED: no header row with both Location and Accession in {SHEET}")

    carried, s1 = {}, {}
    for r in rows[rows.index(hdr) + 1:]:
        d = dict(zip(hdr, r))
        for k in FILL:
            if d.get(k):
                carried[k] = d[k]
            else:
                d[k] = carried.get(k, "")
        if d.get("Accession"):
            s1[d["Accession"]] = d

    ena_tsv = ENA / f"{BIOPROJECT}.tsv"
    if not ena_tsv.exists():
        sys.exit(f"missing {ena_tsv}; run: sources.py ena {BIOPROJECT}")
    ena = list(csv.DictReader(ena_tsv.open(), delimiter="\t"))
    if not ena or "secondary_sample_accession" not in ena[0]:
        sys.exit(f"{ena_tsv} predates the secondary_sample_accession field; "
                 f"refetch with: sources.py ena {BIOPROJECT}")

    out, unjoined = [], []
    for r in ena:
        d = s1.get(r["secondary_sample_accession"])
        if d is None:
            unjoined.append(r["run_accession"])
            continue
        country = d["Country"].strip()
        loc = d["Location"].strip()
        place = f"{country}: {loc}" if loc and loc != "-" else country
        out.append(dict(Run=r["run_accession"], IsolateID=d["Isolate ID"], Location=place,
                        CollectionDate=iso_date(d["Date collected"])))

    if len(out) != N_EXPECTED:
        sys.exit(f"REFUSED: expected {N_EXPECTED} joined runs, got {len(out)}. "
                 f"Unjoined runs: {unjoined[:8]}")
    fine = sum(1 for r in out if ":" in r["Location"])
    if fine != N_WITH_LOCALITY:
        sys.exit(f"REFUSED: expected {N_WITH_LOCALITY} runs with a locality, got {fine}. "
                 f"Table S1 changed; recheck before trusting the values.")
    undated = [r["Run"] for r in out if not r["CollectionDate"]]
    if undated:
        sys.exit(f"REFUSED: {len(undated)} rows have an unparsed date: {undated[:5]}")

    print(f"{len(out)} runs joined; {fine} gain a locality finer than a country",
          file=sys.stderr)
    print("  " + "  ".join(f"{k}:{v}" for k, v in collections.Counter(
        r["Location"].split(":")[0] for r in out).most_common()), file=sys.stderr)
    finish(out, OUT, key_col="Run", bioprojects=[BIOPROJECT], sort_key=lambda r: r["Run"])


if __name__ == "__main__":
    main()
