#!/usr/bin/env python3
"""
03_find.py — fetch BioProject metadata for all BioProjects in crypt.tsv.

Fetches titles, submission dates, study design, and linked PubMed articles for
all 1,797 BioProjects (single + co-infected) to enable field prevalence analysis.
Use --hc to restrict co-infection counts to same_genus_secondary=False runs.

PMIDs are sourced via six strategies (short-circuit, stop at first hit):
  1. BioProject XML <Publication> elements — free, extracted from title fetch
  2. PMC full-text search (40% yield)
  3. Europe PMC (35% yield; catches preprints + supplementary-method mentions)
  4. ENA XML API for PRJEB accessions — submitter-supplied PUBMED_ID
  5. Semantic Scholar (S2_API_KEY required; skipped without key)
  6. PubMed text search (last resort)

Output (output/03_find/):
  data/bioproject_meta.tsv   one row per BioProject
  data/find_cache.json       Entrez + EBI + S2 cache (v4; resumable)
  logs/find.log
  logs/find_summary.txt

Usage:
  python 03_find.py
  python 03_find.py --hc
"""

import argparse
import csv
from datetime import datetime
import json
import os
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from _util import _Tee, http_get, link_latest, load_json, make_log_dir, save_json

# ── Settings ──────────────────────────────────────────────────────────────────

CRYPT_TSV = Path("output/02_filter/data/crypt.tsv")
OUT_DIR   = Path("output/03_find")

CACHE_VERSION = 4

API_KEY    = os.environ.get("NCBI_API_KEY", "")
S2_API_KEY = os.environ.get("S2_API_KEY", "")
RATE       = 9.0 if API_KEY else 2.5
HEADERS    = {"User-Agent": "crypt/03_meta (leon.lenzo@curtin.edu.au)"}
ENTREZ     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

EPMC        = "https://www.ebi.ac.uk/europepmc/webservices/rest"
ENA         = "https://www.ebi.ac.uk/ena/browser/api"
EXT_HEADERS = {"User-Agent": "crypt/03_find (leon.lenzo@curtin.edu.au)"}

OUTPUT_FIELDS = [
    "BioProject", "modes", "n_runs", "n_coinf", "n_single", "coinf_rate",
    "bp_submission_date", "study_design", "pmid_source", "title",
    "primary_pmid", "primary_pub_date", "primary_publication", "abstract",
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
    gap  = 0.5   # 2 req/s for EBI (EuropePMC, ENA)
    wait = _last_ext + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_ext = time.monotonic()


_last_s2: float = 0.0


def _s2_wait() -> None:
    global _last_s2
    # 10 req/s with key, 1 req/s without — stay comfortably under both
    gap  = 1.1  # S2 rate limit: 1 req/s with or without key
    wait = _last_s2 + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_s2 = time.monotonic()


def _europepmc_search(accession: str) -> list[str]:
    """Europe PMC full-text search for accession string; returns PMIDs found."""
    query = urllib.parse.quote(f'"{accession}"')
    url   = (f"{EPMC}/search?query={query}"
             f"&format=json&pageSize=10"
             f"&fields=pmid,title,firstPublicationDate,source,abstractText")
    _ext_wait()
    try:
        raw  = http_get(url, EXT_HEADERS, no_retry_429=True)
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
        raw  = http_get(url, EXT_HEADERS, no_retry_429=True)
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


def _semantic_scholar_search(accession: str) -> list[str]:
    """Semantic Scholar full-text search for accession string; returns PMIDs found."""
    if not S2_API_KEY:
        return []
    params = urllib.parse.urlencode({
        "query":  accession,
        "fields": "externalIds",
        "limit":  5,
    })
    headers = {
        "User-Agent": "crypt/03_find (leon.lenzo@curtin.edu.au)",
        "x-api-key":  S2_API_KEY,
    }
    _s2_wait()
    try:
        raw  = http_get(
            f"https://api.semanticscholar.org/graph/v1/paper/search?{params}",
            headers,
            no_retry_429=True,
        )
        data = json.loads(raw)
    except Exception:
        return []
    pmids = []
    for paper in data.get("data", []):
        pmid = str((paper.get("externalIds") or {}).get("PubMed", "")).strip()
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


# Strategies in order of expected hit rate.
# Each entry: (label, callable(accession, uid) → list[str])
# uid may be None only for strategies that don't need it.
_STRATEGIES = [
    ("BioProject XML",    lambda acc, uid, xml_pmids: list(xml_pmids)),
    ("PMC full-text",     lambda acc, uid, _: _pmcids_to_pmids(_pmc_search(acc))),
    ("Europe PMC",        lambda acc, uid, _: _europepmc_search(acc)),
    ("ENA XML",           lambda acc, uid, _: _ena_pmids(acc)),
    ("Semantic Scholar",  lambda acc, uid, _: _semantic_scholar_search(acc)),
    ("PubMed",            lambda acc, uid, _: _pubmed_search(acc)),
]


def fetch_bioproject_meta(accession: str, cache: dict) -> dict:
    """
    Fetch title, description, and earliest publication for a BioProject.
    Tries strategies in order; stops on the first that returns a PMID.
    Mutates cache in place; caller persists it.
    """
    if accession in cache and cache[accession].get("_v") == CACHE_VERSION:
        if cache[accession].get("n_papers_found", 0) > 0:
            result = dict(cache[accession])
            result.setdefault("pmid_source", "cached")
            return result
        # Found nothing previously — re-try with current strategy order

    result: dict = {
        "title": "", "description": "",
        "bp_submission_date": "",
        "primary_pmid": "", "primary_pub_date": "",
        "primary_publication": "", "abstract": "",
        "n_papers_found": 0,
        "primary_pub_text": "", "pmid_source": "none",
        "_v": CACHE_VERSION,
    }

    try:
        uid = _bp_uid(accession)
        if uid is None:
            cache[accession] = result
            return result
    except Exception as e:
        print(f"WARNING: esearch failed for {accession}: {e}", flush=True)
        return result

    # Always fetch BioProject XML for title + description
    bp_xml_pmids: set[str] = set()
    try:
        raw  = _get("efetch.fcgi", db="bioproject", id=uid)
        root = ET.fromstring(raw)
        result["title"]       = _xtext(root, ".//ProjectDescr/Title")
        result["description"] = _xtext(root, ".//ProjectDescr/Description")
        sub_el = root.find(".//Submission[@submitted]")
        if sub_el is not None:
            result["bp_submission_date"] = sub_el.get("submitted", "")
        for pub in root.findall(".//Publication"):
            pid = (pub.get("id") or "").strip()
            if pid.isdigit():
                bp_xml_pmids.add(pid)
    except Exception as e:
        print(f"WARNING: efetch failed for {accession}: {e}", flush=True)

    # Try each strategy in order; stop at first hit
    all_pmids: list[str] = []
    pmid_source = "none"
    for name, fn in _STRATEGIES:
        try:
            pmids = fn(accession, uid, bp_xml_pmids)
        except Exception as e:
            print(f"WARNING: {name} failed for {accession}: {e}", flush=True)
            continue
        if pmids:
            all_pmids   = pmids
            pmid_source = name
            break

    result["n_papers_found"] = len(all_pmids)
    result["pmid_source"]    = pmid_source

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
                result["abstract"]            = pri["abstract"]
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
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    sys.stdout = _Tee(log_dir / "find.log")
    link_latest(logs_base, log_dir / "find.log")

    cache_path = OUT_DIR / "data" / "find_cache.json"
    out_tsv    = OUT_DIR / "data" / "bioproject_meta.tsv"

    print(f"NCBI_API_KEY set: {'yes' if API_KEY else 'no'}", flush=True)

    if not CRYPT_TSV.exists():
        sys.exit(f"ERROR: {CRYPT_TSV} not found — run 02_filter.py first")

    with open(CRYPT_TSV, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"Loaded {len(rows):,} runs from {CRYPT_TSV}", flush=True)

    if args.hc:
        coinf_rows = [
            r for r in rows
            if r.get("co_infection_flag") != "single"
            and r.get("same_genus_secondary") == "False"
        ]
        label = "high-confidence (same_genus_secondary=False)"
    else:
        coinf_rows = [r for r in rows if r.get("co_infection_flag") != "single"]
        label      = "all co-infected"

    print(f"Co-infection filter: {label} → {len(coinf_rows):,} runs", flush=True)

    # Aggregate ALL runs per BioProject (single + co-infected)
    bp_agg: dict[str, dict] = defaultdict(
        lambda: {"n_coinf": 0, "n_single": 0,
                 "modes": set(), "primaries": set(), "secondaries": set()}
    )
    for r in rows:
        bp = r.get("BioProject", "").strip()
        if not bp or bp in ("", "NA"):
            continue
        is_coinf = r.get("co_infection_flag") != "single"
        # Only count runs that pass the hc filter as coinf when --hc is set
        if args.hc:
            is_coinf = is_coinf and r.get("same_genus_secondary") == "False"
        if is_coinf:
            bp_agg[bp]["n_coinf"] += 1
            bp_agg[bp]["modes"].add(r.get("mode", ""))
            if r.get("primary_pathogen"):
                bp_agg[bp]["primaries"].add(r["primary_pathogen"].strip())
            for sp in r.get("secondary_pathogens", "").split(";"):
                sp = sp.strip()
                if sp:
                    bp_agg[bp]["secondaries"].add(sp)
        else:
            bp_agg[bp]["n_single"] += 1

    bioprojects = sorted(bp_agg,
                         key=lambda b: bp_agg[b]["n_coinf"], reverse=True)
    print(f"BioProjects to fetch: {len(bioprojects)}", flush=True)

    cache    = load_json(cache_path)
    n_cached = sum(1 for bp in bioprojects
                   if bp in cache and cache[bp].get("_v") == CACHE_VERSION)
    print(f"Already cached (v{CACHE_VERSION}): {n_cached}", flush=True)

    results = []
    for i, bp in enumerate(bioprojects, 1):
        n_coinf  = bp_agg[bp]["n_coinf"]
        n_single = bp_agg[bp]["n_single"]
        n_total  = n_coinf + n_single
        ts       = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] [{i:>4}/{len(bioprojects)}] {bp} — "
              f"{n_total} runs, {n_single} infected, {n_coinf} co-infected",
              end=" · ", flush=True)

        try:
            meta = fetch_bioproject_meta(bp, cache)
        except Exception as e:
            print(f"ERROR: {e}", flush=True)
            meta = {"title": "", "primary_pmid": "", "primary_pub_date": "",
                    "primary_publication": "", "n_papers_found": 0,
                    "primary_pub_text": "", "pmid_source": "error"}

        src      = meta.get("pmid_source", "none")
        cached   = src == "cached"
        modes_str = "+".join(sorted(bp_agg[bp]["modes"] - {""}))

        results.append({
            "BioProject":          bp,
            "modes":               modes_str,
            "n_runs":              n_total,
            "n_coinf":             n_coinf,
            "n_single":            n_single,
            "coinf_rate":          round(n_coinf / n_total, 4) if n_total else "",
            "bp_submission_date":  meta.get("bp_submission_date", ""),
            "study_design":        _study_design(meta["title"],
                                                 meta.get("description", ""),
                                                 meta.get("primary_pub_text", "")),
            "pmid_source":         src,
            "title":               meta["title"],
            "primary_pmid":        meta.get("primary_pmid", ""),
            "primary_pub_date":    meta.get("primary_pub_date", ""),
            "primary_publication": meta.get("primary_publication", ""),
            "abstract":            meta.get("abstract", ""),
            "n_papers_found":      meta.get("n_papers_found", 0),
            "primaries":           "; ".join(sorted(bp_agg[bp]["primaries"])),
            "secondaries":         "; ".join(sorted(bp_agg[bp]["secondaries"])),
        })

        print(f"source: {src}", flush=True)
        if not cached:
            save_json(cache, cache_path)

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    n_with_pmids = sum(1 for r in results if r["primary_pmid"])
    n_with_title = sum(1 for r in results if r["title"])
    design_counts: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    for r in results:
        design_counts[r["study_design"]] = design_counts.get(r["study_design"], 0) + 1
        source_counts[r["pmid_source"]]  = source_counts.get(r["pmid_source"],  0) + 1

    n_both_modes = sum(1 for r in results if "+" in r["modes"])

    n_coinf_bp  = sum(1 for r in results if r["n_coinf"] > 0)
    n_single_bp = sum(1 for r in results if r["n_coinf"] == 0)
    source_lines = "\n".join(
        f"  {src:<20} {cnt}"
        for src, cnt in sorted(source_counts.items(), key=lambda x: -x[1])
    )
    summary = (
        f"03_find summary\n"
        f"Co-infection filter: {label}\n"
        f"BioProjects total:   {len(results)}\n"
        f"  with co-infection: {n_coinf_bp}\n"
        f"  single only:       {n_single_bp}\n"
        f"  from MAL only:     {sum(1 for r in results if r['modes'] == 'mal')}\n"
        f"  from HAL only:     {sum(1 for r in results if r['modes'] == 'hal')}\n"
        f"  from both modes:   {n_both_modes}\n"
        f"With title:          {n_with_title}\n"
        f"With PMIDs:          {n_with_pmids}\n"
        f"PMID source breakdown:\n{source_lines}\n"
        f"Study design (auto):\n"
        f"  coinf_experiment:  {design_counts.get('coinf_experiment', 0)}\n"
        f"  field_survey:      {design_counts.get('field_survey', 0)}\n"
        f"  unclear:           {design_counts.get('unclear', 0)}\n"
        f"  (coinf_experiment → manually verify then exclude from novel count)\n"
    )
    summary_path = log_dir / "find_summary.txt"
    summary_path.write_text(summary)
    link_latest(logs_base, summary_path)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}")


if __name__ == "__main__":
    main()
