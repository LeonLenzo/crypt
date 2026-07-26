#!/usr/bin/env python3
"""
05_meta.py — fetch BioProject metadata for co-infection BioProjects.

Reads the co-infection TSV produced by 04_crypt.py and fetches BioProject
titles, descriptions, and linked PubMed articles via Entrez.

By default targets high-confidence co-infections only (same_genus_secondary
== "False"), which are the most credible signal and the primary focus for
distinguishing genuine cryptic detections from k-mer bleed artefacts.
Use --all to include same-genus secondaries as well.

Output (output/05_meta/):
  {mode}_bioproject_meta.tsv   one row per BioProject, sorted by run count
  {mode}_meta_cache.json       Entrez response cache (enables resumability)
  {mode}.log
  {mode}_summary.txt

Usage:
  python 05_meta.py --mode mal
  python 05_meta.py --mode hal
  python 05_meta.py --mode mal --hc    # restrict to same_genus_secondary=False
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

CRYPT_DIR = Path("output/04_crypt/data")
OUT_DIR   = Path("output/05_meta")

API_KEY = os.environ.get("NCBI_API_KEY", "")
RATE    = 9.0 if API_KEY else 2.5
HEADERS = {"User-Agent": "crypt/05_meta (leon.lenzo@curtin.edu.au)"}
ENTREZ  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

OUTPUT_FIELDS = [
    "BioProject", "n_runs", "study_design", "title",
    "primary_pmid", "primary_pub_date", "primary_publication",
    "n_papers_found", "primaries", "secondaries",
]

# ── Study design inference ─────────────────────────────────────────────────────

_COINF_KEYWORDS = {
    "co-infection", "coinfection", "co infection", "dual infection",
    "mixed infection", "dual inoculation", "co-inoculation", "coinoculation",
    "double infection", "co-inoculated", "co-infected",
}
_FIELD_KEYWORDS = {
    "field", "survey", "surveillance", "epidemiology", "epidemiological",
    "natural infection", "naturally infected", "wild", "farm", "orchard",
    "commercial", "pathogenomics",
}


def _study_design(title: str, description: str, pub_text: str = "") -> str:
    text = (title + " " + description + " " + pub_text).lower()
    if any(kw in text for kw in _COINF_KEYWORDS):
        return "coinf_experiment"
    if any(kw in text for kw in _FIELD_KEYWORDS):
        return "field_survey"
    return "unclear"

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


_MONTH_MAP = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
              "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
              "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}


def _pub_date(article: ET.Element) -> str:
    """Return YYYY-MM-DD for a PubmedArticle element, or '9999-12-31' if unknown."""
    # Electronic pub date is often earlier than print
    for el in article.findall(".//ArticleDate[@DateType='Electronic']"):
        yr = el.findtext("Year")
        mo = _MONTH_MAP.get(el.findtext("Month") or "", el.findtext("Month") or "01")
        dy = (el.findtext("Day") or "01").zfill(2)
        if yr and yr.isdigit():
            return f"{yr}-{mo.zfill(2)}-{dy}"
    # Journal pub date
    pd = article.find(".//JournalIssue/PubDate")
    if pd is not None:
        yr  = pd.findtext("Year")
        raw = pd.findtext("Month") or "01"
        mo  = _MONTH_MAP.get(raw, raw.zfill(2) if raw.isdigit() else "01")
        dy  = pd.findtext("Day") or "01"
        dy  = dy.zfill(2) if dy.isdigit() else "01"
        if yr and yr.isdigit():
            return f"{yr}-{mo}-{dy}"
        ml = pd.findtext("MedlineDate") or ""
        if ml[:4].isdigit():
            return f"{ml[:4]}-01-01"
    return "9999-12-31"


# ── Entrez fetchers ────────────────────────────────────────────────────────────

def _bp_uid(accession: str) -> str | None:
    """Return Entrez UID for a BioProject accession, or None if not found."""
    raw  = _get("esearch.fcgi", db="bioproject",
                term=f"{accession}[Project Accession]", retmax=1)
    root = ET.fromstring(raw)
    ids  = [el.text for el in root.findall(".//Id") if el.text]
    return ids[0] if ids else None


def _pmc_search(accession: str) -> list[str]:
    """Full-text search PMC for accession string. Returns PMCIDs."""
    raw  = _get("esearch.fcgi", db="pmc", term=f'"{accession}"', retmax=50)
    root = ET.fromstring(raw)
    return [el.text for el in root.findall(".//Id") if el.text]


def _pmcids_to_pmids(pmcids: list[str]) -> list[str]:
    """Get each PMC article's own PMID via esummary (avoids elink which returns references)."""
    if not pmcids:
        return []
    raw  = _get("esummary.fcgi", db="pmc", id=",".join(pmcids), retmode="json")
    data = json.loads(raw)
    pmids = []
    for pmcid in pmcids:
        doc = data.get("result", {}).get(pmcid, {})
        for aid in doc.get("articleids", []):
            if aid.get("idtype") == "pmid" and aid.get("value"):
                pmids.append(str(aid["value"]))
                break
    return pmids


def _pubmed_search(accession: str) -> list[str]:
    """Search PubMed title/abstract for accession string. Returns PMIDs."""
    raw  = _get("esearch.fcgi", db="pubmed", term=f'"{accession}"', retmax=50)
    root = ET.fromstring(raw)
    return [el.text for el in root.findall(".//Id") if el.text]


def fetch_bioproject_meta(accession: str, cache: dict) -> dict:
    """
    Fetch title, description, PMIDs, and publication titles for a BioProject.
    PMIDs are collected from three sources: elink bioproject→pubmed, PMC
    full-text search, and PubMed text search.  Abstracts are pulled and used
    for study_design inference but not written to the TSV.
    Mutates cache in place; caller is responsible for persisting it.
    """
    if accession in cache and cache[accession].get("_v") == 4:
        return cache[accession]

    result: dict = {"title": "", "description": "",
                    "primary_pmid": "", "primary_pub_date": "",
                    "primary_publication": "", "n_papers_found": 0,
                    "primary_pub_text": "", "_v": 4}

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

    # PMIDs from three sources
    all_pmids: set[str] = set()

    # Source 1: elink bioproject→pubmed
    try:
        raw  = _get("elink.fcgi", dbfrom="bioproject", db="pubmed", id=uid)
        root = ET.fromstring(raw)
        all_pmids.update(
            el.text for el in root.findall(".//LinkSetDb/Link/Id") if el.text)
    except Exception as e:
        print(f"WARNING: elink failed for {accession}: {e}", flush=True)

    # Source 2: PMC full-text search for accession string
    try:
        pmcids = _pmc_search(accession)
        all_pmids.update(_pmcids_to_pmids(pmcids))
    except Exception as e:
        print(f"WARNING: PMC search failed for {accession}: {e}", flush=True)

    # Source 3: PubMed text search for accession string
    try:
        all_pmids.update(_pubmed_search(accession))
    except Exception as e:
        print(f"WARNING: PubMed search failed for {accession}: {e}", flush=True)

    result["n_papers_found"] = len(all_pmids)

    # Fetch titles + pub dates; identify the earliest paper (most likely the depositing paper)
    if all_pmids:
        try:
            raw  = _get("efetch.fcgi", db="pubmed",
                        id=",".join(sorted(all_pmids)), rettype="xml", retmode="xml")
            root = ET.fromstring(raw)
            pub_info: dict[str, dict] = {}   # pmid → {title, date, abstract}
            for article in root.findall(".//PubmedArticle"):
                pid_el = article.find(".//MedlineCitation/PMID")
                ttl_el = article.find(".//MedlineCitation/Article/ArticleTitle")
                abs_el = article.find(".//MedlineCitation/Article/Abstract/AbstractText")
                if pid_el is None:
                    continue
                pid = pid_el.text or ""
                pub_info[pid] = {
                    "title":    (ttl_el.text or "").strip() if ttl_el is not None else "",
                    "date":     _pub_date(article),
                    "abstract": (abs_el.text  or "").strip() if abs_el  is not None else "",
                }
            # Primary = earliest publication date
            if pub_info:
                primary_pid = min(pub_info, key=lambda p: pub_info[p]["date"])
                pri = pub_info[primary_pid]
                result["primary_pmid"]        = primary_pid
                result["primary_pub_date"]    = pri["date"]
                result["primary_publication"] = f"[{primary_pid}] {pri['title']}"
                result["primary_pub_text"]    = f"{pri['title']} {pri['abstract']}"
        except Exception as e:
            print(f"WARNING: pubmed efetch failed for {accession}: {e}", flush=True)

    cache[accession] = result
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["mal", "hal"])
    ap.add_argument("--hc", action="store_true",
                    help="restrict to high-confidence runs (same_genus_secondary=False)")
    args = ap.parse_args()
    mode = args.mode

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "logs").mkdir(parents=True, exist_ok=True)
    sys.stdout = _Tee(OUT_DIR / "logs" / f"{mode}.log")

    crypt_tsv  = CRYPT_DIR / f"{mode}_crypt.tsv"
    cache_path = OUT_DIR / "data" / f"{mode}_meta_cache.json"
    out_tsv    = OUT_DIR / "data" / f"{mode}_bioproject_meta.tsv"

    print(f"Mode: {mode.upper()}")
    print(f"NCBI_API_KEY set: {'yes' if API_KEY else 'no'}", flush=True)

    if not crypt_tsv.exists():
        sys.exit(f"ERROR: {crypt_tsv} not found — run 04_crypt.py --mode {mode} first")

    # Load co-infection table
    with open(crypt_tsv, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"Loaded {len(rows):,} runs from {crypt_tsv.name}")

    # Filter to target set
    if args.hc:
        target = [
            r for r in rows
            if r.get("co_infection_flag") != "single"
            and r.get("same_genus_secondary") == "False"
        ]
        label = "high-confidence (same_genus_secondary=False)"
    else:
        target = [r for r in rows if r.get("co_infection_flag") != "single"]
        label  = "all co-infected"

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
            meta   = {"title": "", "primary_pmid": "", "primary_pub_date": "",
                      "primary_publication": "", "n_papers_found": 0,
                      "primary_pub_text": ""}
            status = "error"

        results.append({
            "BioProject":          bp,
            "n_runs":              bp_agg[bp]["n_runs"],
            "study_design":        _study_design(meta["title"], meta.get("description", ""),
                                                 meta.get("primary_pub_text", "")),
            "title":               meta["title"],
            "primary_pmid":        meta.get("primary_pmid", ""),
            "primary_pub_date":    meta.get("primary_pub_date", ""),
            "primary_publication": meta.get("primary_publication", ""),
            "n_papers_found":      meta.get("n_papers_found", 0),
            "primaries":           "; ".join(sorted(bp_agg[bp]["primaries"])),
            "secondaries":         "; ".join(sorted(bp_agg[bp]["secondaries"])),
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

    n_with_pmids = sum(1 for r in results if r["primary_pmid"])
    n_with_title = sum(1 for r in results if r["title"])
    design_counts = {}
    for r in results:
        design_counts[r["study_design"]] = design_counts.get(r["study_design"], 0) + 1

    summary = (
        f"05_meta {mode.upper()} summary\n"
        f"Filter:              {label}\n"
        f"BioProjects:         {len(results)}\n"
        f"With title:          {n_with_title}\n"
        f"With PMIDs:          {n_with_pmids}\n"
        f"Study design (auto):\n"
        f"  coinf_experiment:  {design_counts.get('coinf_experiment', 0)}\n"
        f"  field_survey:      {design_counts.get('field_survey', 0)}\n"
        f"  unclear:           {design_counts.get('unclear', 0)}\n"
        f"  (coinf_experiment → manually verify then exclude from novel count)\n"
    )
    (OUT_DIR / "logs" / f"{mode}_summary.txt").write_text(summary)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}")


if __name__ == "__main__":
    main()
