#!/usr/bin/env python3
"""Screen the manuscripts on host, by searching title and abstract for the target crops.

    python release/01_cohort/host_screen.py

Reads  raw/raw_papers.tsv
Writes curation/host_screen_papers.tsv — one row per paper, with the terms behind the call.

A paper is on target if its title or abstract names wheat, maize, rice or sorghum. There is
no middle category: it names a target crop or it does not.

An intermediate tier was tried and dropped. Papers using only generic language ("cereal
crops", "phyllosphere") looked like they needed rescuing, because Colvin & Chopra reads that
way and is linked to PRJEB67964 and PRJEB83049, 1,570 analysed runs. It is a REANALYSIS: the
primary sources are Tor23d and Man24b, which name maize and sorghum outright. Checked against
every cohort accession, the generic tier retained exactly one that a named paper did not, and
that one is PRJNA721299, foxtail millet, contributing no analysed runs.

`text_searched` records whether an abstract was available. 14 papers are judged on title
alone, which is a weaker exclusion and is marked rather than hidden; none of them is a source
for a cohort accession.

VALIDATED, 2026-11-02, against the 67 papers then known to supply a cohort accession: 64 were
retained. The three it dropped were all accounted for -- two reanalyses whose accessions come
from named primary sources, and the foxtail-millet deposit contributing no analysed runs -- so
39 of 41 cohort accessions were reachable through a named primary-source paper. Reinstate that
check as a test once this pipeline produces its own cohort; it cannot run against a cohort it
has not built yet.
"""

import collections
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPERS = HERE / "raw" / "raw_papers.tsv"
OUT = HERE / "curation" / "host_screen_papers.tsv"

# Common name, genus, and the abbreviated binomial that follows a first mention.
CROPS = {
    "wheat": [r"wheat", r"triticum\b", r"\bt\. ?(aestivum|turgidum|durum|spelta)\b"],
    "maize": [r"maize", r"\bcorn\b", r"zea mays", r"\bz\. ?mays\b"],
    "rice": [r"\brice\b", r"oryza\b", r"\bo\. ?sativa\b", r"paddy"],
    "sorghum": [r"sorghum", r"\bs\. ?bicolor\b", r"broomcorn", r"\bmilo\b"],
}

# Phrases carrying a crop word that do not mean the crop. Removed before matching, and
# counted, so the screen can be audited instead of trusted.
DECOYS = [r"buckwheat", r"rice university", r"cornell", r"corn salad", r"wheatgrass",
          r"rice ?bean"]

DECOY_RE = re.compile("|".join(DECOYS), re.I)
CROP_RE = {c: re.compile("|".join(p), re.I) for c, p in CROPS.items()}


def normalise(text):
    """Lowercase, mend words broken across a line, collapse whitespace."""
    t = (text or "").lower().replace("­", "")
    t = re.sub(r"(\w)-\s+(\w)", r"\1\2", t)
    return re.sub(r"\s+", " ", t)


def screen(title, abstract):
    text = normalise(f"{title} {abstract}")
    decoys = sorted({m.group(0) for m in DECOY_RE.finditer(text)})
    text = DECOY_RE.sub(" ", text)
    crops, terms = [], []
    for crop, rx in CROP_RE.items():
        found = sorted({m.group(0).strip() for m in rx.finditer(text)})
        if found:
            crops.append(crop)
            terms.extend(found)
    return crops, sorted(set(terms)), decoys


def main():
    rows = []
    with open(PAPERS, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            crops, terms, decoys = screen(r.get("title", ""), r.get("abstract", ""))
            rows.append({
                "doi": r.get("doi", ""),
                "on_target": int(bool(crops)),
                "crops": ";".join(crops),
                "terms": ";".join(terms),
                "text_searched": "title+abstract" if (r.get("abstract") or "").strip()
                                 else "title only",
                "decoys_removed": ";".join(decoys)})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    on = sum(r["on_target"] for r in rows)
    print(f"-> {OUT.name}  {len(rows)} papers")
    print(f"   on target  {on:>5}")
    print(f"   off target {len(rows) - on:>5}")
    title_only = sum(1 for r in rows if r["text_searched"] == "title only")
    print(f"   judged on title alone, no abstract retrieved  {title_only}")
    crops = collections.Counter(c for r in rows for c in r["crops"].split(";") if c)
    print(f"   crops named: {dict(crops.most_common())}")
    decoys = collections.Counter(d for r in rows for d in r["decoys_removed"].split(";") if d)
    if decoys:
        print(f"   decoys removed: {dict(decoys.most_common())}")


if __name__ == "__main__":
    main()
