#!/usr/bin/env python3
"""Per-isolate collection locality and date for PRJEB36444, from Corredor-Moreno et al. 2021.

    Corredor-Moreno P, Minter F, Davey PE, et al. 2021. "The branched-chain amino acid
    aminotransferase TaBCAT1 modulates amino acid metabolism and positively regulates wheat
    rust susceptibility." The Plant Cell 33:1728.  doi:10.1093/plcell/koab049
    Supplemental Data Set 1 (tpc.00530.2020-s01.xlsx), pulled from Europe PMC PMC8254495.

All 26 runs sat on country centroids: ENA carries `country` and nothing finer. Supplemental
Data Set 1 is the Pst sample table and gives the town for every one of them.

## The key

`Accession number` is the ENA SECONDARY sample accession (ERS). The sheet spans three
studies - Hubbard 2015 (SRX, 33 rows), Bueno-Sancho 2017 (ERS, 97 rows) and this one (26
rows marked `This study`) - so the `Data reference` column must be filtered on, or the join
reaches for 154 accessions of which only 26 are in this deposit.

## The sheet has a key error, and two runs are left alone because of it

`ERS4270296` appears on TWO rows, 15.0040 Cambridge and 15.0130 Whattisham, and
`ERS4270297` appears on none even though ENA has that sample. The rows are otherwise in
strict ascending accession order and 15.0130 sits exactly where 297 would, so the intended
mapping is almost certainly 15.0040 -> 296 and 15.0130 -> 297. That is an inference about a
typo, not something the sheet says, and both towns are 50 km apart in the UK, so BOTH rows
are dropped rather than guessed. 24 of 26 runs gain a locality; the two keep their country.

## Corroboration

`14.0083 Bohl-Iggelheim` in the Bueno-Sancho block independently confirms the reading of
PRJEB15280's `Germany: Bohl Bauch 10` as Böhl-Iggelheim plus a field code, which was made
from the project's sampling footprint alone earlier the same day. `Louville-la-Chenard`
appears here too, as it does in Lewis 2022 and PRJEB39201.

    python 01_literature/resolvers/corredor21.py
"""

from __future__ import annotations

import collections, csv, re, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import ENA, STUDIES
from _resolver import finish

BIOPROJECT = "PRJEB36444"
STUDY = STUDIES / "doi_10.1093_plcell_koab049"
XLSX = STUDY / "tpc.00530.2020-s01.xlsx"
SHEET = "Supplemental Data Set 1"
OUT = STUDY / "corredor21_runs.csv"

N_EXPECTED = 24          # 26 runs, less the two the sheet's duplicated accession makes ambiguous
AMBIGUOUS = "ERS4270296"


def iso_date(v: str) -> str:
    """Three date forms in one column: 2015-06-01 00:00:00, 17/06/2016, 2014."""
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
        sys.exit(f"missing {XLSX}; Supplemental Data Set 1 of doi:10.1093/plcell/koab049")

    rows = [["" if c is None else str(c).strip() for c in r]
            for r in openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
            [SHEET].iter_rows(values_only=True)]
    hdr = next((r for r in rows if "Location" in r and "Accession number" in r), None)
    if hdr is None:
        sys.exit(f"REFUSED: no header row with both Location and Accession number in {SHEET}")

    mine = [d for d in (dict(zip(hdr, r)) for r in rows[rows.index(hdr) + 1:])
            if "This study" in (d.get("Data reference") or "")]
    dupes = [a for a, n in collections.Counter(
        d["Accession number"] for d in mine).items() if n > 1]
    if dupes != [AMBIGUOUS]:
        sys.exit(f"REFUSED: expected exactly {AMBIGUOUS} to be duplicated in the sheet, "
                 f"found {dupes}. The key error changed; re-read before trusting this.")
    s1 = {d["Accession number"]: d for d in mine if d["Accession number"] != AMBIGUOUS}

    ena_tsv = ENA / f"{BIOPROJECT}.tsv"
    if not ena_tsv.exists():
        sys.exit(f"missing {ena_tsv}; run: sources.py ena {BIOPROJECT}")
    ena = list(csv.DictReader(ena_tsv.open(), delimiter="\t"))
    if not ena or "secondary_sample_accession" not in ena[0]:
        sys.exit(f"{ena_tsv} predates the secondary_sample_accession field; "
                 f"refetch with: sources.py ena {BIOPROJECT}")

    out, skipped = [], []
    for r in ena:
        d = s1.get(r["secondary_sample_accession"])
        if d is None:
            skipped.append(r["secondary_sample_accession"])
            continue
        loc, country = d["Location"].strip(), d["Country"].strip()
        # One row contradicts itself: Sample ID 15.0383 has Country 'UK' and Location
        # 'Saultain, NW France'. Saultain is a commune in Nord and is already in the cohort as
        # "France: Saultain" from PRJEB39201, so the Location text is right and the Country
        # column is wrong. Left to itself the row is unplaceable, because the geocoder derives
        # the expected country from the string's prefix and then rejects a French answer.
        if loc.endswith(", NW France"):
            if country != "UK":
                sys.exit(f"REFUSED: the Saultain row's Country is now {country!r}, not 'UK'. "
                         f"The sheet changed; re-read before keeping this correction.")
            loc, country = loc[:-len(", NW France")], "France"
        if loc in ("", "-", "unknown"):
            sys.exit(f"REFUSED: {r['run_accession']} has no Location in the sheet; this "
                     f"resolver assumes every 'This study' row carries one")
        out.append(dict(Run=r["run_accession"], SampleID=d["Sample ID"],
                        Location=f"{country}: {loc}",
                        CollectionDate=iso_date(d["Collection date"])))

    if len(out) != N_EXPECTED:
        sys.exit(f"REFUSED: expected {N_EXPECTED} joined runs, got {len(out)}. "
                 f"Unjoined samples: {skipped}")
    undated = [r["Run"] for r in out if not r["CollectionDate"]]
    if undated:
        sys.exit(f"REFUSED: {len(undated)} rows have an unparsed date: {undated[:5]}")

    print(f"{len(out)} runs joined, all with a locality; {len(skipped)} left alone because the "
          f"sheet reuses {AMBIGUOUS} and never writes ERS4270297", file=sys.stderr)
    print("  " + "  ".join(f"{k}:{v}" for k, v in collections.Counter(
        r["Location"].split(":")[0] for r in out).most_common()), file=sys.stderr)
    # 24 of 26 by design, see the module docstring.
    finish(out, OUT, key_col="Run", bioprojects=[BIOPROJECT], min_rate=0.85,
           sort_key=lambda r: r["Run"])


if __name__ == "__main__":
    main()
