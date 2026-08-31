#!/usr/bin/env python3
"""
metadata/figures/prep_lit_resolution.py — prep data for lit_resolution_alluvial.R

Reads meta_search.py + meta_text.py outputs, classifies every BioProject by
(a) how its DOI was found (or wasn't) and (b) how its full text was retrieved
(or wasn't), and writes one row per BioProject to a TSV for ggalluvial.

Classification matches lit_resolution_sankey.py (Plotly version, retired in
favour of this + lit_resolution_alluvial.R — same underlying logic, just
feeding an R plotting script instead of building the figure directly):
  - Serper only ever runs when NCBI found nothing at all (pmid_source=="none")
  - "DOI only" is an NCBI-XML outcome (bare <Publication> DOI, no PMID),
    never touches Serper/CrossRef discovery — CrossRef only enriches it
  - every non-"none" pmid_source has a doi (0 exceptions) except 14
    migration-era cache entries that carry a doi despite pmid_source=="none"
    — folded into "NCBI (DOI only)" here, same as before

Output: metadata/output/figures/sankey/lit_resolution_data.tsv
  columns: BioProject, search, text, outcome

Run from crypt/: python metadata/figures/prep_lit_resolution.py
"""

import csv
import json
from pathlib import Path

BIOPROJECTS_PATH = Path("metadata/output/meta_search/data/bioprojects.json")
TEXT_CACHE_PATH  = Path("metadata/output/meta_text/data/text_cache.jsonl")
OUT_DIR  = Path("metadata/output/figures/sankey")
OUT_TSV  = OUT_DIR / "lit_resolution_data.tsv"
OUT_DIR.mkdir(parents=True, exist_ok=True)

bps = json.loads(BIOPROJECTS_PATH.read_text())

text_cache: dict[str, dict] = {}
with open(TEXT_CACHE_PATH, encoding="utf-8") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        doi, rest = line.split("\t", 1)
        text_cache[doi] = json.loads(rest)   # last line wins (append-only cache)

NCBI_SRCS = {"bioproject_xml", "pmc_search", "pubmed_search"}


def classify_search(v: dict) -> str:
    doi = v.get("doi")
    if not doi:
        return "No DOI found"
    src = v.get("pmid_source") or "none"
    if src in NCBI_SRCS:
        return "NCBI (PMID)"
    if src.startswith("serper"):
        return "Web search"
    return "NCBI (DOI only)"   # doi_only→*, bioproject_xml_doi_only, or "none" w/ a doi


def classify_text(v: dict) -> str:
    doi = v.get("doi")
    if not doi:
        return "No DOI"
    if not v.get("full_text"):
        return "No OA copy found"
    pdf_url = text_cache.get(doi, {}).get("pdf_url", "")
    if pdf_url.startswith("epmc:"):
        return "PMC OA"
    if pdf_url.startswith("manual:"):
        return "Manual PDF"
    return "Unpaywall PDF"


rows = []
for bp, v in bps.items():
    search = classify_search(v)
    text   = classify_text(v)
    outcome = "Full text retrieved" if v.get("full_text") else "No full text available"
    rows.append({"BioProject": bp, "search": search, "text": text, "outcome": outcome})

with open(OUT_TSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["BioProject", "search", "text", "outcome"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"Written {OUT_TSV}  ({len(rows):,} BioProjects)")
from collections import Counter
print("\nsearch:")
for k, n in Counter(r["search"] for r in rows).most_common():
    print(f"  {k:<16} {n:>5,}")
print("\ntext:")
for k, n in Counter(r["text"] for r in rows).most_common():
    print(f"  {k:<18} {n:>5,}")
print("\noutcome:")
for k, n in Counter(r["outcome"] for r in rows).most_common():
    print(f"  {k:<24} {n:>5,}")
