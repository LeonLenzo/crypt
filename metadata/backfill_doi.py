#!/usr/bin/env python3
"""
backfill_doi.py — populate primary_doi for BPs that already have a PMID.

For each cache entry that has a primary_pmid but no primary_doi, fetches
the DOI from PubMed efetch ArticleIdList. Batch-fetches in groups of 100
to minimise API calls. Updates lit_cache.json and re-exports literature.json.

Usage:
    python metadata/backfill_doi.py [--dry-run]
"""

import argparse
import json
import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys_path = str(Path(__file__).resolve().parent.parent)
import sys; sys.path.insert(0, sys_path)
from _util import load_json, save_json

CACHE_PATH = Path("metadata/output/fetch_lit/data/lit_cache.json")
LIT_OUT    = Path("metadata/output/fetch_lit/data/literature.json")

API_KEY = os.environ.get("NCBI_API_KEY", "")
RATE    = 9.0 if API_KEY else 2.5
ENTREZ  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"User-Agent": "crypt/backfill_doi (leon.lenzo@curtin.edu.au)"}

_LIT_FIELDS = [
    "primary_pmid", "primary_doi", "primary_pub_date", "primary_publication",
    "abstract", "methods_text", "pmcid", "n_papers_found", "pmid_source",
]

_last_req: float = 0.0


def _rate_wait() -> None:
    global _last_req
    gap  = 1.0 / RATE
    wait = _last_req + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_req = time.monotonic()


def _get(path: str, **params) -> bytes:
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{ENTREZ}/{path}?{urllib.parse.urlencode(params)}"
    _rate_wait()
    return urllib.request.urlopen(
        urllib.request.Request(url, headers=HEADERS), timeout=30
    ).read()


def fetch_dois_for_pmids(pmids: list[str]) -> dict[str, str]:
    """Return {pmid: doi} for a batch of PMIDs via efetch."""
    if not pmids:
        return {}
    try:
        raw  = _get("efetch.fcgi", db="pubmed",
                    id=",".join(pmids), rettype="xml", retmode="xml")
        root = ET.fromstring(raw)
    except Exception as e:
        print(f"  WARNING: efetch failed: {e}", flush=True)
        return {}
    result: dict[str, str] = {}
    for article in root.findall(".//PubmedArticle"):
        pid_el = article.find(".//MedlineCitation/PMID")
        if pid_el is None:
            continue
        pid = pid_el.text or ""
        doi = ""
        for aid in article.findall(".//ArticleIdList/ArticleId"):
            if aid.get("IdType") == "doi":
                doi = (aid.text or "").strip()
                break
        result[pid] = doi
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cache = load_json(CACHE_PATH) if CACHE_PATH.exists() else {}

    # Find BPs with PMID but no DOI
    to_backfill = {
        bp: entry["primary_pmid"]
        for bp, entry in cache.items()
        if entry.get("primary_pmid") and not entry.get("primary_doi")
    }
    print(f"BPs with PMID but no DOI: {len(to_backfill)}")

    if args.dry_run:
        print("(dry-run — no changes written)")
        return

    # Batch fetch in groups of 100
    BATCH = 100
    items  = list(to_backfill.items())
    n_filled = 0
    for start in range(0, len(items), BATCH):
        batch = items[start:start + BATCH]
        pmids = [pmid for _, pmid in batch]
        print(f"  Fetching DOIs for PMIDs {start+1}–{start+len(batch)} …", flush=True)
        doi_map = fetch_dois_for_pmids(pmids)
        for bp, pmid in batch:
            doi = doi_map.get(pmid, "")
            if doi:
                cache[bp]["primary_doi"] = doi
                n_filled += 1
        save_json(cache, CACHE_PATH)

    print(f"\nFilled: {n_filled}/{len(to_backfill)} BPs got a DOI")
    print(f"  (remaining {len(to_backfill) - n_filled} PMIDs have no DOI in PubMed — "
          f"preprints or data notes)")

    # Re-export literature.json
    lit_out = {}
    for bp, entry in cache.items():
        lit_out[bp] = {k: entry.get(k, "") for k in _LIT_FIELDS}
    LIT_OUT.write_text(json.dumps(lit_out, indent=2))
    print(f"Written: {LIT_OUT}")


if __name__ == "__main__":
    main()
