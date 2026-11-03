#!/usr/bin/env python3
"""Per-isolate collection locality for PRJEB84390, from Tsushima et al. 2025 Table S13.

    Tsushima A, et al. 2025. "Chromosome-level genome assemblies and comparative genomics
    of Puccinia triticina reveal ..." BMC Genomics.  doi:10.1186/s12864-025-12230-4
    Additional file 1 (xlsx), sheet "Table S13".

PRJEB84390 is the largest country-only block in the cohort: 100 wheat RNA-seq runs that
entered as `field` with nothing finer than a country, because ENA's `location` column is
EMPTY for all 171 runs in the project. The archive genuinely has nothing. Table S13 is the
per-isolate table behind the paper's "100 RNA-seq datasets ... across 11 countries", and it
carries a Location column the deposit does not.

## The key

Table S13's `Accession` is the ENA SECONDARY sample accession (`ERS22906149`), not the
`SAMEA...` form, and not a run. `sources.py ena` gained `secondary_sample_accession` on
2026-10-12 so this join needs no private fetch; refetch any project cached before that date.

## What it yields

40 of the 100 gain a locality. The other 60 give `-` in Table S13 itself, so the paper
confirms the country is all that was recorded; those keep the country they already have
rather than being written back as an empty string, which would read as lost coverage.

## Two things not fixed here

The Polish place names are mangled in the supplied xlsx, apparently by a font substitution:
`Jaroseawiec` for Jarosławiec and `Guvcz` for Gurcz, both of which are ALREADY in the cohort
with those spellings from PRJEB15280 and PRJEB39201, the same European surveillance network.
`Chrzelive`, `Pagow`, `Gova-Jciow` and `Janvszeno` are presumably the same corruption and are
not guessed at. Spelling repair belongs in `seeds/missing_localities.tsv`, not here: a
resolver emits what its source says, so the audit chain back to the supplement stays intact.

`80200 Marchelepot` keeps its postcode; `Marchelepot` is already placed in the cohort and the
geocoder's own cleanup handles the prefix.

## Join denominator

The project mixes two data types: 71 Puccinia DNA isolate libraries and the 100 plant
RNA-seq runs. Only the latter can carry a plant collection site, so 100/171 = 58% is the
CEILING, not a shortfall, and min_rate is set below the shared bar for that reason. The
exact count is asserted instead.

    python 01_literature/resolvers/tsu25.py
"""

from __future__ import annotations

import collections, csv, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import ENA, STUDIES
from _resolver import finish

BIOPROJECT = "PRJEB84390"
STUDY = STUDIES / "doi_10.1186_s12864-025-12230-4"
XLSX = STUDY / "12864_2025_12230_MOESM1_ESM.xlsx"
SHEET = "Table S13"
OUT = STUDY / "tsu25_runs.csv"

N_EXPECTED = 100         # every Table S13 row joins; the other 71 runs are pathogen DNA
N_WITH_LOCALITY = 40     # the rest give "-" in the supplement itself

# Table S13 writes the country as the paper does, not as the archive does.
COUNTRY_AS_DEPOSITED = {"Türkiye": "Turkey"}


def main() -> None:
    import openpyxl
    if not XLSX.exists():
        sys.exit(f"missing {XLSX}; Additional file 1 of doi:10.1186/s12864-025-12230-4")

    rows = list(openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
                [SHEET].iter_rows(values_only=True))
    hdr = None
    s13 = {}
    for r in rows:
        v = ["" if c is None else str(c).strip() for c in r]
        if hdr is None:
            # The sheet carries a title line above the header, so the header is found by
            # content rather than by position; a fixed row index broke on Table S1.
            if "Location" in v and "Accession" in v:
                hdr = v
            continue
        if not any(v):
            continue
        d = dict(zip(hdr, v))
        if d.get("Accession"):
            s13[d["Accession"]] = d
    if hdr is None:
        sys.exit(f"REFUSED: no header row with both Location and Accession in {SHEET}")

    ena_tsv = ENA / f"{BIOPROJECT}.tsv"
    if not ena_tsv.exists():
        sys.exit(f"missing {ena_tsv}; run: sources.py ena {BIOPROJECT}")
    ena = list(csv.DictReader(ena_tsv.open(), delimiter="\t"))
    if "secondary_sample_accession" not in (ena[0] if ena else {}):
        sys.exit(f"{ena_tsv} predates the secondary_sample_accession field; "
                 f"refetch with: sources.py ena {BIOPROJECT}")
    run_of = {r["secondary_sample_accession"]: r for r in ena}

    out, unmatched = [], []
    for ers, d in s13.items():
        r = run_of.get(ers)
        if r is None:
            unmatched.append(ers)
            continue
        country = COUNTRY_AS_DEPOSITED.get(d.get("Country", ""), d.get("Country", "")).strip()
        loc = (d.get("Location") or "").strip()
        place = f"{country}: {loc}" if loc and loc != "-" else country
        out.append(dict(Run=r["run_accession"], IsolateID=d.get("Isolate name", ""),
                        Location=place, CollectionDate=(d.get("Collection year") or "").strip()))

    if len(out) != N_EXPECTED:
        sys.exit(f"REFUSED: expected {N_EXPECTED} joined runs, got {len(out)}. "
                 f"Unjoined Table S13 accessions: {unmatched[:8]}")
    fine = sum(1 for r in out if ":" in r["Location"])
    if fine != N_WITH_LOCALITY:
        sys.exit(f"REFUSED: expected {N_WITH_LOCALITY} runs with a locality, got {fine}. "
                 f"Table S13 changed; recheck before trusting the values.")

    print(f"{len(out)} runs joined; {fine} gain a locality finer than a country, "
          f"{len(out) - fine} give '-' in the supplement", file=sys.stderr)
    print("  " + "  ".join(f"{k}:{v}" for k, v in collections.Counter(
        r["Location"].split(":")[0] for r in out).most_common()), file=sys.stderr)
    # 100 of 171: the balance is pathogen DNA, see the module docstring.
    finish(out, OUT, key_col="Run", bioprojects=[BIOPROJECT], min_rate=0.55,
           sort_key=lambda r: r["Run"])


if __name__ == "__main__":
    main()
