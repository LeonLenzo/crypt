#!/usr/bin/env python3
"""
04_meta.py — fetch BioProject metadata for co-infection BioProjects.

Reads the co-infection TSV produced by 03_crypt.py and fetches BioProject
titles, descriptions, and linked PubMed articles via Entrez.

By default targets high-confidence co-infections only (same_genus_secondary
== "False"), which are the most credible signal and the primary focus for
distinguishing genuine cryptic detections from k-mer bleed artefacts.
Use --all to include same-genus secondaries as well.

Output (output/04_meta/):
  {mode}_bioproject_meta.tsv   one row per BioProject, sorted by run count
  {mode}_meta_cache.json       Entrez response cache (enables resumability)
  {mode}.log
  {mode}_summary.txt

Usage:
  python 04_meta.py --mode mal
  python 04_meta.py --mode hal
  python 04_meta.py --mode mal --all
"""

import argparse
import csv
import json
import os
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from _util import _Tee, http_get, load_json, save_json

# ── Settings ──────────────────────────────────────────────────────────────────

CRYPT_DIR = Path("output/03_crypt")
OUT_DIR   = Path("output/04_meta")

API_KEY = os.environ.get("NCBI_API_KEY", "")
RATE    = 9.0 if API_KEY else 2.5
HEADERS = {"User-Agent": "crypt/04_meta (leon.lenzo@curtin.edu.au)"}
ENTREZ  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

OUTPUT_FIELDS = [
    "BioProject", "n_runs", "title", "description",
    "pmids", "publications", "primaries", "secondaries",
]

# ── Rate limiting ──────────────────────────────────────────────────────────────

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
    return http_get(url, HEADERS)


# ── XML helpers ────────────────────────────────────────────────────────────────

def _xtext(root: ET.Element, xpath: str) -> str:
    el = root.find(xpath)
    return (el.text or "").strip() if el is not None else ""


# ── Entrez fetchers ────────────────────────────────────────────────────────────

def _bp_uid(accession: str) -> str | None:
    """Return Entrez UID for a BioProject accession, or None if not found."""
    raw  = _get("esearch.fcgi", db="bioproject",
                term=f"{accession}[Project Accession]", retmax=1)
    root = ET.fromstring(raw)
    ids  = [el.text for el in root.findall(".//Id") if el.text]
    return ids[0] if ids else None


def fetch_bioproject_meta(accession: str, cache: dict) -> dict:
    """
    Fetch title, description, PMIDs, and publication titles for a BioProject.
    Mutates cache in place; caller is responsible for persisting it.
    """
    if accession in cache:
        return cache[accession]

    result: dict = {"title": "", "description": "", "pmids": [], "publications": []}

    try:
        uid = _bp_uid(accession)
        if uid is None:
            cache[accession] = result
            return result
    except Exception as e:
        print(f"WARNING: esearch failed for {accession}: {e}", flush=True)
        return result

    # BioProject title + description
    try:
        raw  = _get("efetch.fcgi", db="bioproject", id=uid)
        root = ET.fromstring(raw)
        result["title"]       = _xtext(root, ".//ProjectDescr/Title")
        result["description"] = _xtext(root, ".//ProjectDescr/Description")
    except Exception as e:
        print(f"WARNING: efetch failed for {accession}: {e}", flush=True)

    # Linked PubMed articles
    try:
        raw   = _get("elink.fcgi", dbfrom="bioproject", db="pubmed", id=uid)
        root  = ET.fromstring(raw)
        pmids = [el.text for el in root.findall(".//LinkSetDb/Link/Id") if el.text]
        result["pmids"] = pmids
    except Exception as e:
        print(f"WARNING: elink failed for {accession}: {e}", flush=True)

    # PubMed titles
    if result["pmids"]:
        try:
            raw  = _get("efetch.fcgi", db="pubmed",
                        id=",".join(result["pmids"]), rettype="xml", retmode="xml")
            root = ET.fromstring(raw)
            pub_map: dict[str, str] = {}
            for article in root.findall(".//PubmedArticle"):
                pid_el   = article.find(".//MedlineCitation/PMID")
                title_el = article.find(".//MedlineCitation/Article/ArticleTitle")
                if pid_el is not None and title_el is not None:
                    pub_map[pid_el.text or ""] = (title_el.text or "").strip()
            result["publications"] = [
                f"[{pid}] {pub_map.get(pid, '')}" for pid in result["pmids"]
            ]
        except Exception as e:
            print(f"WARNING: pubmed efetch failed for {accession}: {e}", flush=True)

    cache[accession] = result
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["mal", "hal"])
    ap.add_argument("--all", action="store_true",
                    help="include all co-infected runs, not just high-confidence")
    args = ap.parse_args()
    mode = args.mode

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sys.stdout = _Tee(OUT_DIR / f"{mode}.log")

    crypt_tsv  = CRYPT_DIR / f"{mode}_crypt.tsv"
    cache_path = OUT_DIR / f"{mode}_meta_cache.json"
    out_tsv    = OUT_DIR / f"{mode}_bioproject_meta.tsv"

    print(f"Mode: {mode.upper()}")
    print(f"NCBI_API_KEY set: {'yes' if API_KEY else 'no'}", flush=True)

    if not crypt_tsv.exists():
        sys.exit(f"ERROR: {crypt_tsv} not found — run 03_crypt.py --mode {mode} first")

    # Load co-infection table
    with open(crypt_tsv, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"Loaded {len(rows):,} runs from {crypt_tsv.name}")

    # Filter to target set
    if args.all:
        target = [r for r in rows if r.get("co_infection_flag") != "single"]
        label  = "all co-infected"
    else:
        target = [
            r for r in rows
            if r.get("co_infection_flag") != "single"
            and r.get("same_genus_secondary") == "False"
        ]
        label = "high-confidence (same_genus_secondary=False)"

    print(f"Filter: {label} → {len(target):,} runs", flush=True)

    # Aggregate per BioProject
    bp_agg: dict[str, dict] = defaultdict(
        lambda: {"n_runs": 0, "primaries": set(), "secondaries": set()}
    )
    for r in target:
        bp = r.get("BioProject", "").strip()
        if not bp or bp in ("", "NA"):
            continue
        bp_agg[bp]["n_runs"] += 1
        if r.get("primary_pathogen"):
            bp_agg[bp]["primaries"].add(r["primary_pathogen"].strip())
        for sp in r.get("secondary_pathogens", "").split(";"):
            sp = sp.strip()
            if sp:
                bp_agg[bp]["secondaries"].add(sp)

    bioprojects = sorted(bp_agg, key=lambda b: bp_agg[b]["n_runs"], reverse=True)
    print(f"BioProjects to fetch: {len(bioprojects)}", flush=True)

    cache = load_json(cache_path)
    n_cached = sum(1 for bp in bioprojects if bp in cache)
    print(f"Already cached: {n_cached}", flush=True)

    results = []
    for i, bp in enumerate(bioprojects, 1):
        cached = bp in cache
        tag    = " (cached)" if cached else ""
        print(f"  [{i:>3}/{len(bioprojects)}] {bp} "
              f"({bp_agg[bp]['n_runs']} runs){tag}", end=" ... ", flush=True)

        try:
            meta = fetch_bioproject_meta(bp, cache)
            status = "ok"
        except Exception as e:
            print(f"ERROR: {e}", flush=True)
            meta   = {"title": "", "description": "", "pmids": [], "publications": []}
            status = "error"

        results.append({
            "BioProject":   bp,
            "n_runs":       bp_agg[bp]["n_runs"],
            "title":        meta["title"],
            "description":  meta["description"],
            "pmids":        "; ".join(meta["pmids"]),
            "publications": "; ".join(meta["publications"]),
            "primaries":    "; ".join(sorted(bp_agg[bp]["primaries"])),
            "secondaries":  "; ".join(sorted(bp_agg[bp]["secondaries"])),
        })

        if not cached:
            print(status, flush=True)
            save_json(cache, cache_path)
        else:
            print(status, flush=True)

    # Write TSV
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    n_with_pmids = sum(1 for r in results if r["pmids"])
    n_with_title = sum(1 for r in results if r["title"])

    summary = (
        f"04_meta {mode.upper()} summary\n"
        f"Filter:      {label}\n"
        f"BioProjects: {len(results)}\n"
        f"With title:  {n_with_title}\n"
        f"With PMIDs:  {n_with_pmids}\n"
    )
    (OUT_DIR / f"{mode}_summary.txt").write_text(summary)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}")


if __name__ == "__main__":
    main()
