#!/usr/bin/env python3
"""Per-isolate collection locality for PRJEB65589, from Lewis et al. 2024 Dataset S2.

    Lewis CM, Morier-Gxoyiya C, Hubbard A, Nellist CF, Bebber DP, Saunders DGO. 2024.
    "Resurgence of wheat stem rust infections in western Europe: causes and how to curtail
    them." New Phytologist.  doi:10.1111/nph.19864
    Dataset S2 (separate .xlsx) supplied by leon 2026-10-09.

All 53 runs entered the cohort as `field` on Notes S1, but with no locality finer than a
country: the ENA record carries `country` and nothing else, so 53 samples sat on the map as
country centroids. Dataset S2 is the per-isolate table the article points at, and it has the
collection site.

## The key

ENA `sample_title` is the isolate code - `20-0354`, `21-0118`, `22-002-2` - and Dataset S2's
`Original isolate ID` is the same code. 52 of the 53 join; `22-0031` is absent from the
dataset and is left alone rather than guessed at.

`Latitude and longitude` is present as a column but is `-` for every row, so the place name is
all there is. That is enough: the shared geocoder resolves "Banbury, Oxfordshire" and
"Northeim, Lower Saxony" without help.

Seven of the 52 give `-` for Location. Those keep the country they already have rather than
being written back as an empty string, which would look like a regression in coverage.

Dataset S2 also carries Host: Wheat 50, Rye 2. It is emitted as `HostSpecies` in a CSV,
because `host` is not a settable curation field - per-sample host is read by
host_summary.curated_species(), which takes exactly this shape. The two rye samples are the
f. sp. secalis pair that WHOLE_PROJECT_HOST had folded into wheat.

    python 01a_Literature/resolvers/lewis24.py
"""

from __future__ import annotations

import collections, csv, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import ENA, STUDIES
from _resolver import finish

BIOPROJECT = "PRJEB65589"
XLSX = STUDIES / "doi_10.1111_nph.19864" / "nph19864-sup-0002-datasets2.xlsx"
OUT = STUDIES / "doi_10.1111_nph.19864" / "lewis24_runs.csv"
N_EXPECTED = 52          # 53 runs, one isolate absent from Dataset S2


def main() -> None:
    import openpyxl
    if not XLSX.exists():
        sys.exit(f"missing {XLSX}")
    rows = [r for r in openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
            ["Dataset S2"].iter_rows(values_only=True)]
    hdr = [str(c).strip() if c else "" for c in rows[1]]
    s2 = {}
    for r in rows[2:]:
        d = dict(zip(hdr, r))
        k = str(d.get("Original isolate ID") or "").strip()
        if k:
            s2[k] = d

    ena_tsv = ENA / f"{BIOPROJECT}.tsv"
    if not ena_tsv.exists():
        sys.exit(f"missing {ena_tsv}; run: sources.py ena {BIOPROJECT}")

    out, unmatched = [], []
    for r in csv.DictReader(ena_tsv.open(), delimiter="\t"):
        title = (r.get("sample_title") or "").strip()
        d = s2.get(title)
        if d is None:
            unmatched.append(title)
            continue
        loc = str(d.get("Location") or "").strip()
        country = (r.get("country") or "").strip()
        # "-" means Dataset S2 has no site for this isolate. Keep the country rather than
        # writing an empty Location back over it.
        place = f"{country}: {loc}" if loc and loc != "-" else country
        host = str(d.get("Host") or "").strip()
        out.append(dict(Run=r["run_accession"], IsolateID=title, Location=place,
                        HostSpecies={"Wheat": "Triticum aestivum", "Rye": "Secale cereale"}.get(host, host),
                        CollectionDate=str(d.get("Date collected") or "").strip()))

    if len(out) != N_EXPECTED:
        sys.exit(f"REFUSED: expected {N_EXPECTED} joined runs, got {len(out)}. "
                 f"Unmatched sample_titles: {unmatched}")
    if unmatched:
        print(f"note: {len(unmatched)} run(s) absent from Dataset S2 and left alone: "
              f"{unmatched}", file=sys.stderr)
    fine = sum(1 for r in out if ":" in r["Location"])
    print(f"{len(out)} runs joined; {fine} gain a locality finer than a country",
          file=sys.stderr)
    print("  hosts: " + "  ".join(f"{k}:{v}" for k, v in
          collections.Counter(r["HostSpecies"] for r in out).most_common()), file=sys.stderr)
    finish(out, OUT, key_col="Run", bioprojects=[BIOPROJECT], sort_key=lambda r: r["Run"])


if __name__ == "__main__":
    main()
