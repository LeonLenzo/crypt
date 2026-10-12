#!/usr/bin/env python3
"""Per-isolate collection town for PRJEB36485, from Visser et al. 2021 Table S2.

    Visser B, Meyer M, Park RF, Gilligan CA, Burgin LE, Hort MC, Hodson DP, Pretorius ZA.
    2021. "Microsatellite analysis and urediniospore dispersal simulations support the
    movement of Puccinia graminis f. sp. tritici from southern Africa to Australia."
    Plant Pathology 70:1875.  doi:10.1111/ppa.13468
    Supplemental Table S2 (ppa13468-sup-0004-tables2.xlsx), supplied by leon 2026-10-12.

Wiley returns 403 to every automated route and the article has no PMC record, so this is
the one supplement in the 2026-10-12 pass that had to be fetched by hand.

## The key

`Isolate ID` is the ENA `sample_title` verbatim - `14/SADL1`, `15/SAZP7`, `14/K4`, `SA1` -
and all 49 of the project's titles are present in Table S2. No accession mapping needed,
which is a first for this group's supplements.

## Only the field arm is emitted

The deposit runs two arms and `visser21-isolates-setting` already separates them: 39 runs
named year/site are field samples, 10 named isolates (SA1-SA4, KE74217, KE89069, ET03b/10,
ET08/10, ER179b/11, ER181a/11) are wheat seedlings inoculated with purified isolates, dated
1974 to 2011. Table S2 gives a `Country of isolation` for those 10 as well, but writing it
as their `location` would record where the ISOLATE was found as if it were where the
sequenced plant grew. That is the trap that put PRJNA1231252 in the field cohort, so the
10 are skipped.

## Two rows are skipped because Ada21 holds a finer value

Table S2 says Bredasdorp for `14/SADL1` and `14/SADL2`; the Ada21 supplement says Ouplaas, a
settlement about 45 km away in the same Overberg district. Ada21 is also the only source
with a day-precision date for exactly those samples (2014/08/23, against Table S2's bare
"2014"), so it is carrying the original collection record while Table S2 rounds to the
district town. The two are skipped rather than overwritten; `14/SADL5` is Bredasdorp in both
and is emitted.

## What it yields

37 runs emitted, 24 of them with a town. The 13 Kenyan samples give "-" in Table S2 itself.

    python 01a_Literature/resolvers/visser21.py
"""

from __future__ import annotations

import collections, csv, re, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import ENA, STUDIES
from _resolver import finish

BIOPROJECT = "PRJEB36485"
STUDY = STUDIES / "doi_10.1111_ppa.13468"
XLSX = STUDY / "ppa13468-sup-0004-tables2.xlsx"
SHEET = "Table S2"
OUT = STUDY / "visser21_runs.csv"

FIELD_TITLE = re.compile(r"^\d{2}/")      # 14/SADL1, 15/SAZP7, 14/K4
KEEP_ADA21 = {"14/SADL1", "14/SADL2"}     # Ouplaas, finer than Table S2's Bredasdorp

N_EXPECTED = 37           # 49 runs, less 10 inoculated isolates and the 2 above
N_WITH_TOWN = 24          # the 13 Kenyan rows give "-"


def main() -> None:
    import openpyxl
    if not XLSX.exists():
        sys.exit(f"missing {XLSX}; Supplemental Table S2 of doi:10.1111/ppa.13468. "
                 f"Wiley 403s every automated route, so this one is a manual download.")

    rows = [["" if c is None else str(c).strip() for c in r]
            for r in openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
            [SHEET].iter_rows(values_only=True)]
    hdr = next((r for r in rows if "Town of isolation" in r and "Isolate ID" in r), None)
    if hdr is None:
        sys.exit(f"REFUSED: no header row with both Isolate ID and Town of isolation in {SHEET}")
    s2 = {d["Isolate ID"]: d for d in (dict(zip(hdr, r)) for r in rows[rows.index(hdr) + 1:])
          if d.get("Isolate ID")}

    ena_tsv = ENA / f"{BIOPROJECT}.tsv"
    if not ena_tsv.exists():
        sys.exit(f"missing {ena_tsv}; run: sources.py ena {BIOPROJECT}")
    ena = list(csv.DictReader(ena_tsv.open(), delimiter="\t"))

    missing = [r["sample_title"] for r in ena if r["sample_title"] not in s2]
    if missing:
        sys.exit(f"REFUSED: {len(missing)} ENA sample_titles absent from Table S2: {missing[:8]}")

    out, inoculated, deferred = [], [], []
    for r in ena:
        title = r["sample_title"]
        if not FIELD_TITLE.match(title):
            inoculated.append(title)
            continue
        if title in KEEP_ADA21:
            deferred.append(title)
            continue
        d = s2[title]
        town = d["Town of isolation"].strip()
        country = d["Country of isolation"].strip()
        out.append(dict(Run=r["run_accession"], IsolateID=title,
                        Location=f"{country}: {town}" if town and town != "-" else country))

    if sorted(deferred) != sorted(KEEP_ADA21):
        sys.exit(f"REFUSED: expected to defer {sorted(KEEP_ADA21)} to Ada21, saw {deferred}")
    if len(out) != N_EXPECTED:
        sys.exit(f"REFUSED: expected {N_EXPECTED} emitted runs, got {len(out)}. "
                 f"{len(inoculated)} inoculated-isolate titles were skipped: {inoculated}")
    fine = sum(1 for r in out if ":" in r["Location"])
    if fine != N_WITH_TOWN:
        sys.exit(f"REFUSED: expected {N_WITH_TOWN} runs with a town, got {fine}. "
                 f"Table S2 changed; recheck before trusting the values.")

    print(f"{len(out)} field-arm runs emitted, {fine} with a town; "
          f"{len(inoculated)} inoculated isolates and {len(deferred)} Ada21-finer rows skipped",
          file=sys.stderr)
    print("  " + "  ".join(f"{k}:{v}" for k, v in collections.Counter(
        r["Location"] for r in out if ":" in r["Location"]).most_common()), file=sys.stderr)
    # 37 of 49 by design, see the module docstring.
    finish(out, OUT, key_col="Run", bioprojects=[BIOPROJECT], min_rate=0.70,
           sort_key=lambda r: r["Run"])


if __name__ == "__main__":
    main()
