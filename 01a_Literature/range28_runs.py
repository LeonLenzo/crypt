#!/usr/bin/env python3
"""Resolve collection dates for the Range 28 field sorghum drought time series.

PRJNA1119650 is 100 leaf RNA-seq libraries from the University of Arizona Maricopa
Agricultural Center, deposited by Boyce Thompson Institute / Cornell and reported in:

    Yu L et al. 2026. "Integrated Multi-Omic Analyses Uncover a Regulatory Link Between
    Photosynthesis and Drought Tolerance in Field-Grown Sorghum."
    Plant, Cell & Environment.  doi:10.1111/pce.70649   PMC13436492   open access.

The same BioProject is also the deposit of record for a second paper from the same group
(doi:10.1093/plcell/koag196, PMC13412047), which states the identical accession in its data
availability section. One dataset, two papers; neither is a reuse of the other.

## Why a join is needed at all

Every one of the 100 BioSamples carries `collection_date` = 2023-07-13. That is a single
value on a seven-week weekly time series, and it is in the wrong year:

    "samples were collected weekly (on Thursdays) for 7 weeks (week 1-week 7, from 13 August
     to 24 September 2020)"

and the trial is dated independently in the same methods - "Destructive plant phenotyping was
conducted at harvest, beginning on 4 November 2020" - and again in the microbiome section as
"the same 2020 field trial". So the archive date is three years out and flattens the series,
the same failure the Kashima join was built to undo.

What the archive DOES carry, and carries correctly, is the design code, as the BioSample
attribute `source_material_id`:

    A_WL_TP1_1   =   accession A, water-limited, timepoint 1, replicate 1

and the paper names the timepoints it used - "the time point (week 1, week 3, week 5 and
week 7)" - which matches the four TP values present (TP1, TP3, TP5, TP7) exactly. Weekly
Thursdays counted from 13 August 2020 therefore give TP1 = 13 Aug, TP3 = 27 Aug,
TP5 = 10 Sep, TP7 = 24 Sep 2020. The arithmetic is checked below rather than trusted: all
four land on a Thursday, and TP7 lands on the paper's stated end date.

Because `source_material_id` is a BioSample attribute, `_sra_key` finds it without any
regex, and this table keys on it directly.

## What is NOT taken from the archive

`geo_loc_name` reads `Nigeria` on 18 samples and `USA` on 16, blank on the other 66. Those
are the PLANT INTRODUCTION origins of two of the six accessions, not sampling sites: all 18
Nigeria rows are PI 533871 and all 16 USA rows are PI 656053. Every plant grew in one field
in Arizona. This is the germplasm-origin trap recorded in crypt-archive-is-unreliable, and it
is why location is written as a project-wide curation rule from the paper's coordinates.

`dev_stage` reads `Seedling` on all 100 and is simply wrong: the drought treatment began at
flag leaf appearance, ~47 days after planting, and sampling ran from there to maturity. Not
curated here because the cohort does not use dev_stage, but noted so nobody leans on it.

Run:  python 01a_Literature/range28_runs.py
"""

from __future__ import annotations

import collections, csv, datetime, json, re, sys
from pathlib import Path

from _layout import BIOSAMPLE, RUNS, STUDIES

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ATTRS = BIOSAMPLE / "PRJNA1119650.json"
OUT = STUDIES / "doi_10.1111_pce.70649" / "range28_runs.csv"

BIOPROJECT = "PRJNA1119650"

# "samples were collected weekly (on Thursdays) for 7 weeks (week 1-week 7, from 13 August
# to 24 September 2020)". Week 1 is the anchor; the rest is arithmetic, verified below.
WEEK1 = datetime.date(2020, 8, 13)
LAST_WEEK_STATED = datetime.date(2020, 9, 24)

# "Six selected sorghum accessions (PI 533871: M1, PI 533961: DL/59/1530, PI 656041: 80M,
# PI 656053: N290-B, PI 656076: SC1271 and PI 656096: SC391)" and "accession IDs (A: M1,
# B: DL/59/1530, C: 80M, D: N290-B, E: SC1271 and F: SC391)". The two lists are given in the
# same order, which is how letter maps to PI number; the mapping is cross-checked against the
# BioSample `cultivar` attribute below rather than assumed from the ordering alone.
ACCESSION = {
    "A": ("M1", "PI 533871"), "B": ("DL/59/1530", "PI 533961"), "C": ("80M", "PI 656041"),
    "D": ("N290-B", "PI 656053"), "E": ("SC1271", "PI 656076"), "F": ("SC391", "PI 656096"),
}
TREATMENT = {"WW": "well-watered (control)", "WL": "water-limited (drought)"}

CODE = re.compile(r"^([A-F])_(WW|WL)_TP([1-7])_(\d+)$")


def main() -> None:
    # The date mapping is the entire point of this script, so prove the arithmetic before use.
    weeks = {i: WEEK1 + datetime.timedelta(weeks=i - 1) for i in range(1, 8)}
    off = [f"TP{i}={d}" for i, d in weeks.items() if d.weekday() != 3]
    if off:
        sys.exit(f"REFUSED: the paper says sampling was on Thursdays, but {off} are not")
    if weeks[7] != LAST_WEEK_STATED:
        sys.exit(f"REFUSED: week 7 computes to {weeks[7]}, paper says {LAST_WEEK_STATED}")

    attrs = json.load(ATTRS.open())
    if not attrs.get("complete") or attrs.get("sampled") != attrs.get("total"):
        sys.exit(f"REFUSED: BioSample attribute cache is incomplete: "
                 f"{attrs.get('sampled')}/{attrs.get('total')}")
    a = attrs["attrs"]

    rows, seen_cv = [], collections.defaultdict(set)
    for samn, v in sorted(a.items()):
        code = (v.get("source_material_id") or "").strip()
        m = CODE.match(code)
        if not m:
            sys.exit(f"REFUSED: {samn} has source_material_id {code!r}, which does not match "
                     f"the design grammar; the date mapping would silently miss it")
        geno, trt, tp, rep = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        if tp not in weeks:
            sys.exit(f"REFUSED: {samn} is TP{tp}, outside the 7 stated weeks")
        name, pi = ACCESSION[geno]
        seen_cv[geno].add((v.get("cultivar") or "").strip())
        rows.append(dict(
            SourceMaterialID=code, BioSample=samn,
            CollectionDate=weeks[tp].isoformat(), Timepoint=f"TP{tp}", Week=tp,
            GenotypeCode=geno, Accession=name, PlantIntroduction=pi,
            Treatment=trt, TreatmentDesc=TREATMENT[trt], Replicate=rep,
            ArchiveDate=(v.get("collection_date") or ""),
            ArchiveGeoLoc=(v.get("geo_loc_name") or "")))

    # The letter-to-PI mapping comes from two lists printed in the same order. If the
    # BioSample `cultivar` attribute disagrees with it, the ordering assumption is wrong.
    for geno, cvs in sorted(seen_cv.items()):
        expect = ACCESSION[geno][1]
        if cvs != {expect}:
            sys.exit(f"REFUSED: genotype {geno} should be {expect} per the paper's accession "
                     f"list, but BioSample cultivar says {sorted(cvs)}")

    if len({r["SourceMaterialID"] for r in rows}) != len(rows):
        dup = [k for k, n in collections.Counter(r["SourceMaterialID"] for r in rows).items() if n > 1]
        sys.exit(f"REFUSED: source_material_id is the join key and is not unique: {dup}")

    days = sorted({r["CollectionDate"] for r in rows})
    if len(days) < 4:
        sys.exit(f"REFUSED: only {len(days)} distinct date(s); the series did not expand")

    rows.sort(key=lambda r: (r["GenotypeCode"], r["Treatment"], r["Week"], r["Replicate"]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT.relative_to(ROOT)}: {len(rows)} samples", file=sys.stderr)
    print(f"  dates: {', '.join(days)}", file=sys.stderr)
    print(f"  per date: {dict(collections.Counter(r['CollectionDate'] for r in rows))}", file=sys.stderr)
    print(f"  treatment: {dict(collections.Counter(r['TreatmentDesc'] for r in rows))}", file=sys.stderr)
    print(f"  accession: {dict(collections.Counter(r['Accession'] for r in rows))}", file=sys.stderr)
    print(f"  archive said: {dict(collections.Counter(r['ArchiveDate'] for r in rows))}", file=sys.stderr)

    # How many of our runs will the join actually reach?
    ours = [r for r in csv.DictReader(RUNS.open(), delimiter="\t")
            if r["BioProject"] == BIOPROJECT]
    keys = {r["SourceMaterialID"] for r in rows}
    hit = sum(1 for r in ours
              if (a.get(r["BioSample"], {}).get("source_material_id") or "") in keys)
    print(f"  joins {hit}/{len(ours)} of our runs on source_material_id", file=sys.stderr)


if __name__ == "__main__":
    main()
