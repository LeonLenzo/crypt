#!/usr/bin/env python3
"""Collection year for PRJNA357594, decoded from Che et al. 2018 Supplementary Table 3.

    Che Z, Liu H, Yao F, et al. 2018. "Transcriptome analysis reveals the genetic basis of
    the resistance response to Fusarium graminearum in maize." BMC Genomics 19:131.
    doi:10.1186/s12864-018-4513-4
    Additional file 1 (12864_2018_4513_MOESM1_ESM.xlsx), sheet SupplementaryTable3.

All 16 runs had no collection date. The paper's methods say the RNA-seq material came from
"2004 and 2006 field samples" without saying which sample is which, and the BioSample
attributes carry only line, treatment and days-post-inoculation, which gives eight
combinations for sixteen samples. The year is the missing dimension and the archive's
structured fields do not have it.

## Where the year actually is

Supplementary Table 3 names each sample and footnotes the naming scheme outright:

    "Sample names show the name of inbred followed by type of treatment (B = mock,
     Fg = fungal inoculation), sample collection 1 day or 2 days post inoculation,
     and experimental year 2004 or 2006"

So `B73Fg2d04` is B73, fungal-inoculated, 2 dpi, 2004. The GEO title of every run carries
that same code verbatim - `GSM2429847: B73Fg2d04; Zea mays; RNA-Seq` - so the supplement
supplies the KEY and the archive supplies the JOIN. The titles are cached beside this
resolver in `geo_sample_titles.tsv`; they are not in runinfo or in the BioSample attributes,
which is why a plain curation predicate cannot reach them.

## Year, not month, deliberately

Supplementary Table 2 gives inoculation windows - CO441 Aug 4-6 2004, B73 Aug 24-26 2004,
CO441 July 28-30 2006, B73 Aug 14-16 2006 - and collection is 1 or 2 days after. That would
support MONTH precision for three of the four groups, but CO441's 2006 window closes on
30 July, so 1-2 dpi straddles July and August. Rather than record month for three groups and
year for the fourth, all sixteen get the year, which every group supports outright.

    python 01_literature/resolvers/che18.py
"""

from __future__ import annotations

import collections, csv, re, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import STUDIES
from _resolver import finish

BIOPROJECT = "PRJNA357594"
STUDY = STUDIES / "doi_10.1186_s12864-018-4513-4"
XLSX = STUDY / "12864_2018_4513_MOESM1_ESM.xlsx"
SHEET = "SupplementaryTable3"
TITLES = STUDY / "geo_sample_titles.tsv"
OUT = STUDY / "che18_runs.csv"

N_EXPECTED = 16
CODE = re.compile(r"\b((?:B73|CO441)(?:B|Fg)[12]d(0[46]))\b")
FOOTNOTE = "experimental year 2004 or 2006"


def main() -> None:
    import openpyxl
    for p in (XLSX, TITLES):
        if not p.exists():
            sys.exit(f"missing {p}")

    # The footnote is the authority for reading the code, so refuse if it is not there:
    # without it the trailing 04/06 is just two digits.
    text = " ".join(
        str(c) for r in openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
        [SHEET].iter_rows(values_only=True) for c in r if c is not None)
    if FOOTNOTE not in text:
        sys.exit(f"REFUSED: {SHEET} no longer carries the naming footnote "
                 f"({FOOTNOTE!r}); the 04/06 suffix cannot be read as a year without it")
    in_sheet = {m.group(1) for m in CODE.finditer(text)}

    out, unmatched = [], []
    for r in csv.DictReader(TITLES.open(), delimiter="\t"):
        m = CODE.search(r["geo_title"])
        if not m:
            unmatched.append(r["Run"])
            continue
        code, yy = m.group(1), m.group(2)
        if code not in in_sheet:
            sys.exit(f"REFUSED: {r['Run']}'s code {code!r} is not in {SHEET}; the GEO title "
                     f"and the supplement disagree about which samples exist")
        out.append(dict(Run=r["Run"], SampleCode=code, CollectionDate=f"20{yy}"))

    if len(out) != N_EXPECTED:
        sys.exit(f"REFUSED: expected {N_EXPECTED} runs, got {len(out)}. "
                 f"Titles with no decodable code: {unmatched}")
    by_year = collections.Counter(r["CollectionDate"] for r in out)
    if sorted(by_year) != ["2004", "2006"] or len(set(by_year.values())) != 1:
        sys.exit(f"REFUSED: expected a balanced 2004/2006 split, got {dict(by_year)}")

    print(f"{len(out)} runs dated from the Table S3 naming scheme: "
          + "  ".join(f"{k}:{v}" for k, v in sorted(by_year.items())), file=sys.stderr)
    finish(out, OUT, key_col="Run", bioprojects=[BIOPROJECT], sort_key=lambda r: r["Run"])


if __name__ == "__main__":
    main()
