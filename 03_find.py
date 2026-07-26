#!/usr/bin/env python3
"""
03_find.py — fetch BioProject metadata for co-infection BioProjects.

Reads the unified co-infection TSV produced by 02_filter.py (crypt.tsv) and
fetches BioProject titles, descriptions, and linked PubMed articles via Entrez.

By default targets all co-infected runs (co_infection_flag != 'single').
Use --hc to restrict to high-confidence runs (same_genus_secondary == 'False').

PMIDs are collected from seven sources and the earliest publication is
identified as the likely depositing paper:
  1. BioProject XML <Publication> elements — submitter-supplied, zero extra API calls
  2. elink bioproject→pubmed (NCBI directly linked articles)
  3. PMC full-text search for the BioProject accession string
  4. PubMed text search for the BioProject accession string
  (sources 5–7 only run when 1–4 find nothing)
  5. Europe PMC full-text search — catches preprints and supplementary-method mentions
  6. ENA XML API for PRJEB accessions — submitter-supplied PUBMED_ID fields
  7. OpenAlex full-text search — 250M+ works, different indexing from NCBI/EBI

Output (output/03_meta/):
  data/bioproject_meta.tsv   one row per BioProject, sorted by run count
  data/meta_cache.json       Entrez response cache (enables resumability)
  logs/meta.log
  logs/meta_summary.txt

Usage:
  python 03_meta.py
  python 03_meta.py --hc   # restrict to same_genus_secondary=False
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

CRYPT_TSV = Path("output/02_filter/data/crypt.tsv")
OUT_DIR   = Path("output/03_find")

CACHE_VERSION = 4

API_KEY = os.environ.get("NCBI_API_KEY", "")
RATE    = 9.0 if API_KEY else 2.5
HEADERS = {"User-Agent": "crypt/03_meta (leon.lenzo@curtin.edu.au)"}
ENTREZ  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

EPMC        = "https://www.ebi.ac.uk/europepmc/webservices/rest"
ENA         = "https://www.ebi.ac.uk/ena/browser/api"
OPENALEX    = "https://api.openalex.org"
EXT_HEADERS = {"User-Agent": "crypt/03_find (leon.lenzo@curtin.edu.au)"}

OUTPUT_FIELDS = [
    "BioProject", "modes", "n_runs", "study_design", "title",
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


_last_ext: float = 0.0


def _ext_wait() -> None:
    global _last_ext
    gap  = 0.5   # 2 req/s — respectful of EBI servers
    wait = _last_ext + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_ext = time.monotonic()


def _europepmc_search(accession: str) -> list[str]:
    """Europe PMC full-text search for accession string; returns PMIDs found."""
    query = urllib.parse.quote(f'"{accession}"')
    url   = (f"{EPMC}/search?query={query}"
             f"&format=json&pageSize=10"
             f"&fields=pmid,title,firstPublicationDate,source,abstractText")
    _ext_wait()
    try:
        raw  = http_get(url, EXT_HEADERS)
        data = json.loads(raw)
    except Exception:
        return []
    pmids = []
    for r in data.get("resultList", {}).get("result", []):
        pmid = r.get("pmid", "").strip()
        if pmid:
            pmids.append(pmid)
    return pmids


def _ena_pmids(accession: str) -> list[str]:
    """Fetch ENA XML for a PRJEB accession; return any linked PUBMED IDs."""
    if not accession.startswith("PRJEB"):
        return []
    url = f"{ENA}/xml/{accession}"
    _ext_wait()
    try:
        raw  = http_get(url, EXT_HEADERS)
        root = ET.fromstring(raw)
    except Exception:
        return []
    pmids: set[str] = set()
    for el in root.iter("XREF_LINK"):
        if (el.findtext("DB") or "").upper() == "PUBMED":
            pid = (el.findtext("ID") or "").strip()
            if pid:
                pmids.add(pid)
    for el in root.iter("EXTERNAL_ID"):
        if el.get("namespace", "").upper() == "PUBMED":
            pid = (el.text or "").strip()
            if pid:
                pmids.add(pid)
    return list(pmids)


def _openalex_search(accession: str) -> list[str]:
    """OpenAlex full-text search for accession string; returns PMIDs found."""
    params = urllib.parse.urlencode({
        "filter": f'fulltext.search:"{accession}"',
        "per-page": 5,
        "mailto": "leon.lenzo@curtin.edu.au",
        "select": "ids",
    })
    _ext_wait()
    try:
        raw  = http_get(f"{OPENALEX}/works?{params}", EXT_HEADERS)
        data = json.loads(raw)
    except Exception:
        return []
    pmids = []
    for work in data.get("results", []):
        pmid_url = (work.get("ids") or {}).get("pmid", "")
        if pmid_url:
            pmid = pmid_url.rstrip("/").split("/")[-1]
            if pmid.isdigit():
                pmids.append(pmid)
    return pmids


# ── XML helpers ────────────────────────────────────────────────────────────────

def _xtext(root: ET.Element, xpath: str) -> str:
    el = root.find(xpath)
    return (el.text or "").strip() if el is not None else ""


_MONTH_MAP = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
              "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
              "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}


def _pub_date(article: ET.Element) -> str:
    """Return YYYY-MM-DD for a PubmedArticle element, or '9999-12-31' if unknown."""
    for el in article.findall(".//ArticleDate[@DateType='Electronic']"):
        yr = el.findtext("Year")
        mo = _MONTH_MAP.get(el.findtext("Month") or "", el.findtext("Month") or "01")
        dy = (el.findtext("Day") or "01").zfill(2)
        if yr and yr.isdigit():
            return f"{yr}-{mo.zfill(2)}-{dy}"
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
    raw  = _get("esearch.fcgi", db="bioproject",
                term=f"{accession}[Project Accession]", retmax=1)
    root = ET.fromstring(raw)
    ids  = [el.text for el in root.findall(".//Id") if el.text]
    return ids[0] if ids else None


def _pmc_search(accession: str) -> list[str]:
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
    raw  = _get("esearch.fcgi", db="pubmed", term=f'"{accession}"', retmax=50)
    root = ET.fromstring(raw)
    return [el.text for el in root.findall(".//Id") if el.text]


def fetch_bioproject_meta(accession: str, cache: dict) -> dict:
    """
    Fetch title, description, and earliest publication for a BioProject.
    Mutates cache in place; caller persists it.
    """
    if accession in cache and cache[accession].get("_v") == CACHE_VERSION:
        # Re-try external sources for entries that previously found nothing
        if cache[accession].get("n_papers_found", 0) == 0:
            cached = cache[accession]
            extra: set[str] = set()
            try:
                extra.update(_europepmc_search(accession))
            except Exception:
                pass
            try:
                extra.update(_ena_pmids(accession))
            except Exception:
                pass
            try:
                extra.update(_openalex_search(accession))
            except Exception:
                pass
            if not extra:
                return cached
            # Found something new — fall through to full fetch to get pub metadata
        else:
            return cache[accession]

    result: dict = {
        "title": "", "description": "",
        "primary_pmid": "", "primary_pub_date": "",
        "primary_publication": "", "n_papers_found": 0,
        "primary_pub_text": "", "_v": CACHE_VERSION,
    }

    try:
        uid = _bp_uid(accession)
        if uid is None:
            cache[accession] = result
            return result
    except Exception as e:
        print(f"WARNING: esearch failed for {accession}: {e}", flush=True)
        return result

    bp_xml_pmids: set[str] = set()
    try:
        raw  = _get("efetch.fcgi", db="bioproject", id=uid)
        root = ET.fromstring(raw)
        result["title"]       = _xtext(root, ".//ProjectDescr/Title")
        result["description"] = _xtext(root, ".//ProjectDescr/Description")
        # BioProject XML sometimes carries submitter-supplied PMIDs directly
        for pub in root.findall(".//Publication"):
            pid = (pub.get("id") or "").strip()
            if pid.isdigit():
                bp_xml_pmids.add(pid)
    except Exception as e:
        print(f"WARNING: efetch failed for {accession}: {e}", flush=True)

    all_pmids: set[str] = set(bp_xml_pmids)

    try:
        raw  = _get("elink.fcgi", dbfrom="bioproject", db="pubmed", id=uid)
        root = ET.fromstring(raw)
        all_pmids.update(
            el.text for el in root.findall(".//LinkSetDb/Link/Id") if el.text)
    except Exception as e:
        print(f"WARNING: elink failed for {accession}: {e}", flush=True)

    try:
        pmcids = _pmc_search(accession)
        all_pmids.update(_pmcids_to_pmids(pmcids))
    except Exception as e:
        print(f"WARNING: PMC search failed for {accession}: {e}", flush=True)

    try:
        all_pmids.update(_pubmed_search(accession))
    except Exception as e:
        print(f"WARNING: PubMed search failed for {accession}: {e}", flush=True)

    # External sources (5–7) only run when NCBI+BioProjectXML (1–4) found nothing —
    # avoids redundant external calls on re-runs where PMIDs were already found.
    if not all_pmids:
        try:
            all_pmids.update(_europepmc_search(accession))
        except Exception as e:
            print(f"WARNING: Europe PMC search failed for {accession}: {e}", flush=True)

        try:
            all_pmids.update(_ena_pmids(accession))
        except Exception as e:
            print(f"WARNING: ENA XML fetch failed for {accession}: {e}", flush=True)

        try:
            all_pmids.update(_openalex_search(accession))
        except Exception as e:
            print(f"WARNING: OpenAlex search failed for {accession}: {e}", flush=True)

    result["n_papers_found"] = len(all_pmids)

    if all_pmids:
        try:
            raw  = _get("efetch.fcgi", db="pubmed",
                        id=",".join(sorted(all_pmids)), rettype="xml", retmode="xml")
            root = ET.fromstring(raw)
            pub_info: dict[str, dict] = {}
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
    ap.add_argument("--hc", action="store_true",
                    help="restrict to high-confidence runs (same_genus_secondary=False)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "logs").mkdir(parents=True, exist_ok=True)
    sys.stdout = _Tee(OUT_DIR / "logs" / "meta.log")

    cache_path = OUT_DIR / "data" / "meta_cache.json"
    out_tsv    = OUT_DIR / "data" / "bioproject_meta.tsv"

    print(f"NCBI_API_KEY set: {'yes' if API_KEY else 'no'}", flush=True)

    if not CRYPT_TSV.exists():
        sys.exit(f"ERROR: {CRYPT_TSV} not found — run 02_filter.py first")

    with open(CRYPT_TSV, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"Loaded {len(rows):,} runs from {CRYPT_TSV}", flush=True)

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
        lambda: {"n_runs": 0, "modes": set(), "primaries": set(), "secondaries": set()}
    )
    for r in target:
        bp = r.get("BioProject", "").strip()
        if not bp or bp in ("", "NA"):
            continue
        bp_agg[bp]["n_runs"]  += 1
        bp_agg[bp]["modes"].add(r.get("mode", ""))
        if r.get("primary_pathogen"):
            bp_agg[bp]["primaries"].add(r["primary_pathogen"].strip())
        for sp in r.get("secondary_pathogens", "").split(";"):
            sp = sp.strip()
            if sp:
                bp_agg[bp]["secondaries"].add(sp)

    bioprojects = sorted(bp_agg, key=lambda b: bp_agg[b]["n_runs"], reverse=True)
    print(f"BioProjects to fetch: {len(bioprojects)}", flush=True)

    cache    = load_json(cache_path)
    n_cached = sum(1 for bp in bioprojects
                   if bp in cache and cache[bp].get("_v") == CACHE_VERSION)
    print(f"Already cached (v{CACHE_VERSION}): {n_cached}", flush=True)

    results = []
    for i, bp in enumerate(bioprojects, 1):
        cached = bp in cache and cache[bp].get("_v") == CACHE_VERSION
        tag    = " (cached)" if cached else ""
        print(f"  [{i:>3}/{len(bioprojects)}] {bp} "
              f"({bp_agg[bp]['n_runs']} runs){tag}", end=" ... ", flush=True)

        try:
            meta   = fetch_bioproject_meta(bp, cache)
            status = "ok"
        except Exception as e:
            print(f"ERROR: {e}", flush=True)
            meta   = {"title": "", "primary_pmid": "", "primary_pub_date": "",
                      "primary_publication": "", "n_papers_found": 0,
                      "primary_pub_text": ""}
            status = "error"

        modes_str = "+".join(sorted(bp_agg[bp]["modes"] - {""}))

        results.append({
            "BioProject":          bp,
            "modes":               modes_str,
            "n_runs":              bp_agg[bp]["n_runs"],
            "study_design":        _study_design(meta["title"],
                                                 meta.get("description", ""),
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

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    n_with_pmids = sum(1 for r in results if r["primary_pmid"])
    n_with_title = sum(1 for r in results if r["title"])
    design_counts: dict[str, int] = {}
    for r in results:
        design_counts[r["study_design"]] = design_counts.get(r["study_design"], 0) + 1

    n_both_modes = sum(1 for r in results if "+" in r["modes"])

    summary = (
        f"03_meta summary\n"
        f"Filter:              {label}\n"
        f"BioProjects:         {len(results)}\n"
        f"  from MAL only:     {sum(1 for r in results if r['modes'] == 'mal')}\n"
        f"  from HAL only:     {sum(1 for r in results if r['modes'] == 'hal')}\n"
        f"  from both modes:   {n_both_modes}\n"
        f"With title:          {n_with_title}\n"
        f"With PMIDs:          {n_with_pmids}\n"
        f"Study design (auto):\n"
        f"  coinf_experiment:  {design_counts.get('coinf_experiment', 0)}\n"
        f"  field_survey:      {design_counts.get('field_survey', 0)}\n"
        f"  unclear:           {design_counts.get('unclear', 0)}\n"
        f"  (coinf_experiment → manually verify then exclude from novel count)\n"
    )
    (OUT_DIR / "logs" / "meta_summary.txt").write_text(summary)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}")


if __name__ == "__main__":
    main()
