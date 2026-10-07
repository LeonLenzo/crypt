#!/usr/bin/env python3
"""Resolve the EPICON field-droughted sorghum time series from its GEO records.

PRJNA527782 is 395 runs with a completely empty archive record: no tissue, no location, no
date. It is a GEO submission, and GEO holds everything SRA does not.

    Varoquaux N et al. 2019. "Transcriptomic analysis of field-droughted sorghum from
    seedling to maturity reveals biotic and metabolic responses." PNAS 116:27124-27132.
    doi:10.1073/pnas.1907500116   PMC6936495   open access.   GEO: GSE128441

Weekly sampling of two genotypes (RTx430, BTx642) under control watering and two drought
regimes, across the 2016 growing season. The GEO records give, per sample:

    source_name          Leaves (197) or Roots (198)
    characteristics      genotype, watering regime, age (days), plot
    description[1]       a 10-character code whose first six digits are MMDDYY
    molecule             polyA RNA, on all 395

Those dates resolve to 16 weekly timepoints from 2016-06-15 to 2016-09-28, which is the
"tightly-resolved time series" the paper describes and which nothing in SRA records.

## Why half of it is out of scope

198 of the 395 are ROOT. The cohort is aerial tissue only, because soil community makes a
fungal hit in root material as likely a saprotroph or endophyte as a pathogen. Roots are
written with their real tissue value and left out of the field cohort by the tissue rule
rather than by deletion.

Relevant to how this study is read: the authors DID report a fungal finding, "evidence of a
disruption in the plant's symbiosis with arbuscular mycorrhizal fungi". That is a root
symbiont, so it bears on the root half and not the leaf half, but it means this is not a
study that reported nothing about fungi.

## The site needs a cross-reference

The 2019 paper's main text never names the field site; it is in the SI Appendix, whose
Europe PMC zip is truncated and unreadable. The site is instead taken from the group's 2026
companion paper, which names it and states that year 1 followed the same protocol:

    "Sorghum plants were planted, grown, and treated as described by Varoquaux et al. (2019)
     with some differences in years 2 and 3. These included, for year 2, plants being sown in
     pre-wetted plots at Kearny Agricultural Research and Extension Center (Parlier, CA, USA)
     on June 6th, 2017"

Deliberately NOT taken from the author affiliations, which also list Kearney: an affiliation
is where an author works, not where a plant grew, and that inference has already been wrong
once in this project.

Run:  python 01a_Literature/resolvers/epicon.py
"""

from __future__ import annotations

import collections, csv, re, sys, urllib.request
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import RUNS, STUDIES

ROOT = Path(__file__).resolve().parents[2]
STUDY = STUDIES / "doi_10.1073_pnas.1907500116"
CACHE = STUDY / "GSE128441_samples.txt"
OUT = STUDY / "epicon_runs.csv"

GEO = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
       "?acc=GSE128441&targ=gsm&form=text&view=brief")

SITE = "USA: California, Parlier (UC Kearney Agricultural Research and Extension Center)"
SETTING = "field"
TISSUE = {"Leaves": "leaf", "Roots": "root"}


def fetch_geo() -> str:
    if CACHE.exists():
        return CACHE.read_text(encoding="utf-8", errors="replace")
    with urllib.request.urlopen(GEO, timeout=180) as fh:
        txt = fh.read().decode("utf-8", errors="replace")
    STUDY.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(txt)
    print(f"cached {CACHE.relative_to(ROOT)}", file=sys.stderr)
    return txt


def parse(txt: str) -> list[dict]:
    out = []
    for blk in txt.split("^SAMPLE = ")[1:]:
        gsm = blk.split("\n", 1)[0].strip()
        d = {"GSM": gsm, "desc": []}
        for line in blk.splitlines():
            m = re.match(r"!Sample_(\w+) = (.*)$", line)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if k == "characteristics_ch1":
                kk, _, vv = v.partition(":")
                d[kk.strip()] = vv.strip()
            elif k == "description":
                d["desc"].append(v)
            elif k in ("source_name_ch1", "molecule_ch1", "title"):
                d[k] = v
        out.append(d)
    return out


def main() -> None:
    recs = parse(fetch_geo())
    print(f"GEO: {len(recs)} samples", file=sys.stderr)

    # Every library must be polyA mRNA; the cohort's whole detection premise assumes it.
    mol = collections.Counter(r.get("molecule_ch1", "?") for r in recs)
    if set(mol) != {"polyA RNA"}:
        sys.exit(f"REFUSED: expected every sample to be 'polyA RNA', got {dict(mol)}")

    rows, undated = [], 0
    for r in recs:
        src = r.get("source_name_ch1", "")
        tis = TISSUE.get(src, "")
        if not tis:
            sys.exit(f"REFUSED: unknown source_name {src!r} on {r['GSM']}; "
                     f"the tissue mapping is incomplete and roots must not be mislabelled")
        code = r["desc"][1] if len(r["desc"]) > 1 else ""
        m = re.match(r"^(\d{2})(\d{2})(\d{2})", code)
        if m:
            mm, dd, yy = m.groups()
            date = f"20{yy}-{mm}-{dd}"
        else:
            date, undated = "", undated + 1
        rows.append(dict(
            SampleName=r["GSM"], Tissue=tis, CollectionDate=date, Location=SITE,
            Setting=SETTING, SamplingSelection="unselected",
            Genotype=r.get("genotype", ""), WateringRegime=r.get("watering regime", ""),
            AgeDays=r.get("age (days)", ""), Plot=r.get("plot", ""), SourceCode=code))

    # The dates are the point of this study; a silent failure to parse them would flatten the
    # series exactly as the Kashima join once did.
    if undated:
        sys.exit(f"REFUSED: {undated} samples have no parseable MMDDYY code")
    days = sorted({r["CollectionDate"] for r in rows})
    if len(days) < 10:
        sys.exit(f"REFUSED: only {len(days)} distinct dates; the weekly series did not parse")

    rows.sort(key=lambda r: r["SampleName"])
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    tis = collections.Counter(r["Tissue"] for r in rows)
    print(f"\nwrote {OUT.relative_to(ROOT)}: {len(rows)} samples", file=sys.stderr)
    print(f"  tissue: {dict(tis)}", file=sys.stderr)
    print(f"  {len(days)} weekly dates, {days[0]} to {days[-1]}", file=sys.stderr)
    print(f"  genotypes: {dict(collections.Counter(r['Genotype'] for r in rows))}", file=sys.stderr)
    print(f"  regimes: {dict(collections.Counter(r['WateringRegime'] for r in rows))}", file=sys.stderr)

    # How many of ours will actually join?
    ours = [r for r in csv.DictReader(RUNS.open(), delimiter="\t")
            if r["BioProject"] == "PRJNA527782"]
    keys = {r["SampleName"] for r in rows}
    hit = sum(1 for r in ours if r["SampleName"] in keys)
    print(f"  joins {hit}/{len(ours)} of our runs on SampleName", file=sys.stderr)


if __name__ == "__main__":
    main()
