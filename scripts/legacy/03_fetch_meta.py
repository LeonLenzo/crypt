#!/usr/bin/env python3
"""
03_fetch_meta.py — fetch BioProject metadata for all BioProjects in crypt.tsv.

Fetches titles, submission dates, linked PubMed articles, and PMC full-text
methods sections for all BioProjects. All results are written to find_cache.json
for downstream classification in 04_filter_meta.py (no TSV output here).

PMIDs are sourced via six strategies (short-circuit, stop at first hit):
  1. BioProject XML <Publication> elements — free, extracted from title fetch
  2. PMC full-text search (40% yield)
  3. Europe PMC (35% yield; catches preprints + supplementary-method mentions)
  4. ENA XML API for PRJEB accessions — submitter-supplied PUBMED_ID
  5. Semantic Scholar (S2_API_KEY required; skipped without key)
  6. PubMed text search (last resort)

After finding a PMID, the script checks whether the primary paper is in PMC
and fetches the full-text methods section if available. This is stored in
find_cache.json for use by 04_filter_meta.py study-design classification.

Cache behaviour:
  v5 with PMID    → instant skip (already fully processed)
  v5 no PMID      → skip (use --retry to re-run strategies)
  v4 with PMID    → upgrade in-place: re-fetch abstract (fixes structured-abstract
                    parsing bug), fetch PMCID, fetch PMC methods section (~2-3 calls)
  v4 no PMID      → silently bumped to v5 (skip unless --retry)
  not in cache    → full fetch (all six PMID strategies + PMC methods if found)

Output:
  output/03_fetch_meta/data/find_cache.json   Entrez + EBI + S2 cache (v5; resumable)
  output/03_fetch_meta/logs/fetch.log
  output/03_fetch_meta/logs/fetch_summary.txt

Usage:
  python 03_fetch_meta.py           # upgrade v4→v5 + fetch new BPs
  python 03_fetch_meta.py --retry   # also retry no-PMID entries
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
from pathlib import Path

from _util import _Tee, http_get, link_latest, load_json, make_log_dir, save_json

# ── Settings ──────────────────────────────────────────────────────────────────

CRYPT_TSV = Path("output/02_filter_runs/data/crypt.tsv")
OUT_DIR   = Path("output/03_fetch_meta")

CACHE_VERSION = 5

API_KEY    = os.environ.get("NCBI_API_KEY", "")
S2_API_KEY = os.environ.get("S2_API_KEY", "")
RATE       = 9.0 if API_KEY else 2.5
HEADERS    = {"User-Agent": "crypt/03_fetch_meta (leon.lenzo@curtin.edu.au)"}
ENTREZ     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

EPMC        = "https://www.ebi.ac.uk/europepmc/webservices/rest"
ENA         = "https://www.ebi.ac.uk/ena/browser/api"
EXT_HEADERS = {"User-Agent": "crypt/03_fetch_meta (leon.lenzo@curtin.edu.au)"}

# Maximum characters to store from PMC methods section (keeps cache manageable)
METHODS_MAX_CHARS = 8000

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
    gap  = 1.1  # S2 rate limit: 1 req/s with or without key
    wait = _last_s2 + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_s2 = time.monotonic()


# ── External search functions ──────────────────────────────────────────────────

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
        "User-Agent": "crypt/03_fetch_meta (leon.lenzo@curtin.edu.au)",
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


def _pmid_to_pmcid(pmid: str) -> str | None:
    """Return the PMC article ID for a PubMed article, or None if not in PMC."""
    raw  = _get("esearch.fcgi", db="pmc", term=f"{pmid}[pmid]", retmax=1)
    root = ET.fromstring(raw)
    ids  = [el.text for el in root.findall(".//Id") if el.text]
    return ids[0] if ids else None


def _fetch_pmc_methods(pmcid: str) -> str:
    """Fetch PMC full-text XML and extract the methods section text."""
    raw  = _get("efetch.fcgi", db="pmc", id=pmcid, rettype="full", retmode="xml")
    root = ET.fromstring(raw)
    for sec in root.iter("sec"):
        sec_type   = sec.get("sec-type", "").lower()
        title_el   = sec.find("title")
        title_text = (title_el.text or "").lower() if title_el is not None else ""
        if ("method" in sec_type or "material" in sec_type
                or "method" in title_text or "material" in title_text):
            parts = []
            for el in sec.iter():
                if el.tag in ("p", "title") and el.text:
                    parts.append(el.text.strip())
                if el.tail:
                    tail = el.tail.strip()
                    if tail:
                        parts.append(tail)
            text = " ".join(p for p in parts if p)
            return text[:METHODS_MAX_CHARS]
    return ""


# Strategies in order of expected hit rate.
_STRATEGIES = [
    ("BioProject XML",   lambda acc, uid, xml_pmids: list(xml_pmids)),
    ("PMC full-text",    lambda acc, uid, _: _pmcids_to_pmids(_pmc_search(acc))),
    ("Europe PMC",       lambda acc, uid, _: _europepmc_search(acc)),
    ("ENA XML",          lambda acc, uid, _: _ena_pmids(acc)),
    ("Semantic Scholar", lambda acc, uid, _: _semantic_scholar_search(acc)),
    ("PubMed",           lambda acc, uid, _: _pubmed_search(acc)),
]


def fetch_bioproject_meta(accession: str, cache: dict) -> dict:
    """
    Fetch title, description, earliest publication, and PMC methods text for a BioProject.
    Tries PMID strategies in order; stops on the first that returns a PMID.
    Fetches PMC full-text methods section if the primary paper is in PMC.
    Mutates cache in place; caller persists it.
    """
    # v5 cache hit
    if accession in cache and cache[accession].get("_v") == CACHE_VERSION:
        if cache[accession].get("n_papers_found", 0) > 0:
            result = dict(cache[accession])
            result.setdefault("pmid_source", "cached")
            return result
        # v5 with no PMID — retry strategies below

    # v4 → v5 upgrade: entries with a PMID get PMC methods + abstract repair
    if accession in cache and cache[accession].get("_v") == 4:
        entry = cache[accession]
        if entry.get("primary_pmid"):
            entry.setdefault("pmcid", "")
            entry.setdefault("methods_text", "")

            # Re-fetch PubMed abstract if missing (fixes structured abstract parsing bug)
            if not entry.get("abstract", "").strip():
                try:
                    raw  = _get("efetch.fcgi", db="pubmed",
                                id=entry["primary_pmid"], rettype="xml", retmode="xml")
                    root = ET.fromstring(raw)
                    for article in root.findall(".//PubmedArticle"):
                        abs_els  = article.findall(
                            ".//MedlineCitation/Article/Abstract/AbstractText")
                        abstract = " ".join(
                            "".join(el.itertext()).strip() for el in abs_els
                        ).strip()
                        ttl_el   = article.find(".//MedlineCitation/Article/ArticleTitle")
                        title    = ("".join(ttl_el.itertext()).strip()
                                    if ttl_el is not None else entry.get("title", ""))
                        if abstract:
                            entry["abstract"]         = abstract
                            entry["title"]            = title
                            entry["primary_pub_text"] = f"{title} {abstract}"
                        break
                except Exception as e:
                    print(f"WARNING: PubMed re-fetch failed for {accession}: {e}",
                          flush=True)

            # Fetch PMC full-text methods
            if not entry["pmcid"]:
                pmcid = _pmid_to_pmcid(entry["primary_pmid"])
                entry["pmcid"] = pmcid or ""
                if pmcid:
                    try:
                        entry["methods_text"] = _fetch_pmc_methods(pmcid)
                    except Exception as e:
                        print(f"WARNING: PMC fetch failed for {pmcid}: {e}", flush=True)

            entry["_v"]      = CACHE_VERSION
            cache[accession] = entry
            result = dict(entry)
            result["pmid_source"] = "cached"
            return result
        # v4 with no PMID — fall through to full re-fetch

    result: dict = {
        "title": "", "description": "",
        "bp_submission_date": "",
        "primary_pmid": "", "primary_pub_date": "",
        "primary_publication": "", "abstract": "",
        "n_papers_found": 0,
        "primary_pub_text": "", "pmid_source": "none",
        "pmcid": "", "methods_text": "",
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
                pid_el  = article.find(".//MedlineCitation/PMID")
                ttl_el  = article.find(".//MedlineCitation/Article/ArticleTitle")
                abs_els = article.findall(".//MedlineCitation/Article/Abstract/AbstractText")
                if pid_el is None:
                    continue
                pid = pid_el.text or ""
                # Join all AbstractText elements with itertext() to handle structured
                # abstracts (Label sections) and inline markup (<i>, <sup>, etc.)
                abstract = " ".join(
                    "".join(el.itertext()).strip() for el in abs_els
                ).strip()
                # ArticleTitle may also contain inline markup
                title = "".join(ttl_el.itertext()).strip() if ttl_el is not None else ""
                pub_info[pid] = {
                    "title":    title,
                    "date":     _pub_date(article),
                    "abstract": abstract,
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

    # Fetch PMC full-text methods section for the primary paper
    if result["primary_pmid"]:
        try:
            pmcid = _pmid_to_pmcid(result["primary_pmid"])
            result["pmcid"] = pmcid or ""
            if pmcid:
                result["methods_text"] = _fetch_pmc_methods(pmcid)
        except Exception as e:
            print(f"WARNING: PMC methods fetch failed for {accession}: {e}", flush=True)

    cache[accession] = result
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--retry", action="store_true",
                    help="re-run PMID strategies for BPs that previously returned nothing "
                         "(only useful when a new search strategy has been added)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    sys.stdout = _Tee(log_dir / "fetch.log")
    link_latest(logs_base, log_dir / "fetch.log")

    cache_path = OUT_DIR / "data" / "find_cache.json"

    print(f"NCBI_API_KEY set: {'yes' if API_KEY else 'no'}", flush=True)
    print(f"S2_API_KEY set:   {'yes' if S2_API_KEY else 'no'}", flush=True)
    if args.retry:
        print("--retry: will re-run strategies on all no-PMID entries", flush=True)

    if not CRYPT_TSV.exists():
        sys.exit(f"ERROR: {CRYPT_TSV} not found — run 02_filter_runs.py first")

    with open(CRYPT_TSV, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # Collect unique BioProjects in order of first appearance
    seen: set[str] = set()
    bioprojects: list[str] = []
    for r in rows:
        bp = r.get("BioProject", "").strip()
        if bp and bp not in ("", "NA") and bp not in seen:
            seen.add(bp)
            bioprojects.append(bp)

    cache = load_json(cache_path)

    # Classify BPs into work categories before the main loop
    to_upgrade: list[str] = []   # v4 with PMID → abstract repair + PMC methods
    to_fetch:   list[str] = []   # not in cache → full fetch
    to_retry:   list[str] = []   # no PMID → retry strategies (only with --retry)
    n_skip_pmid    = 0            # v5 with PMID → already done
    n_skip_no_pmid = 0            # no PMID, not retrying → skip

    for bp in bioprojects:
        entry = cache.get(bp, {})
        v        = entry.get("_v", 0)
        has_pmid = bool(entry.get("primary_pmid"))

        if v == CACHE_VERSION and has_pmid:
            n_skip_pmid += 1
        elif v == CACHE_VERSION and not has_pmid:
            if args.retry:
                to_retry.append(bp)
            else:
                n_skip_no_pmid += 1
        elif v == 4 and has_pmid:
            to_upgrade.append(bp)
        elif v == 4 and not has_pmid:
            if args.retry:
                to_retry.append(bp)
            else:
                # Silently bump to v5 so it won't be reconsidered on the next run
                entry["_v"] = CACHE_VERSION
                entry.setdefault("pmcid", "")
                entry.setdefault("methods_text", "")
                cache[bp] = entry
                n_skip_no_pmid += 1
        else:
            to_fetch.append(bp)

    total_work = len(to_upgrade) + len(to_fetch) + len(to_retry)
    print(f"BioProjects: {len(bioprojects)} total", flush=True)
    print(f"  skip — v5 with PMID:    {n_skip_pmid}", flush=True)
    print(f"  skip — no PMID found:   {n_skip_no_pmid}"
          + (" (use --retry to force)" if n_skip_no_pmid else ""), flush=True)
    print(f"  upgrade v4→v5:          {len(to_upgrade)}", flush=True)
    print(f"  fetch (new):            {len(to_fetch)}", flush=True)
    if args.retry:
        print(f"  retry (no-PMID):        {len(to_retry)}", flush=True)
    print(f"  → work items:           {total_work}", flush=True)

    if n_skip_no_pmid and not args.retry:
        # Persist the silent v4→v5 bumps before doing real work
        save_json(cache, cache_path)

    n_upgraded = n_fetched = n_retried = 0
    work_list = (
        [("upgrade", bp) for bp in to_upgrade]
        + [("fetch",   bp) for bp in to_fetch]
        + [("retry",   bp) for bp in to_retry]
    )

    for i, (kind, bp) in enumerate(work_list, 1):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] [{i:>4}/{total_work}] {kind:<8} {bp}", end=" · ", flush=True)

        try:
            meta = fetch_bioproject_meta(bp, cache)
        except Exception as e:
            print(f"ERROR: {e}", flush=True)
            continue

        src = meta.get("pmid_source", "none")
        print(f"source: {src} | pmcid: {meta.get('pmcid') or '-'} | "
              f"methods: {'yes' if meta.get('methods_text') else 'no'}", flush=True)

        save_json(cache, cache_path)
        if kind == "upgrade":
            n_upgraded += 1
        elif kind == "retry":
            n_retried += 1
        else:
            n_fetched += 1

    save_json(cache, cache_path)

    n_with_pmid    = sum(1 for bp in bioprojects if cache.get(bp, {}).get("primary_pmid"))
    n_with_methods = sum(1 for bp in bioprojects if cache.get(bp, {}).get("methods_text"))
    source_counts: dict[str, int] = {}
    for bp in bioprojects:
        src = cache.get(bp, {}).get("pmid_source", "none")
        source_counts[src] = source_counts.get(src, 0) + 1
    source_lines = "\n".join(
        f"  {src:<20} {cnt}"
        for src, cnt in sorted(source_counts.items(), key=lambda x: -x[1])
    )

    summary = (
        f"03_fetch_meta summary\n"
        f"BioProjects total:       {len(bioprojects)}\n"
        f"  skipped (v5+PMID):     {n_skip_pmid}\n"
        f"  skipped (no PMID):     {n_skip_no_pmid}\n"
        f"  upgraded v4→v5:        {n_upgraded}\n"
        f"  fetched (new):         {n_fetched}\n"
        f"  retried (no-PMID):     {n_retried}\n"
        f"With PMID:               {n_with_pmid}\n"
        f"With PMC methods text:   {n_with_methods}\n"
        f"PMID source breakdown:\n{source_lines}\n"
        f"\nRun 04_filter_meta.py to classify study designs and write bioproject_meta.tsv\n"
    )
    summary_path = log_dir / "fetch_summary.txt"
    summary_path.write_text(summary)
    link_latest(logs_base, summary_path)
    print(f"\n{summary}")
    print(f"Cache written: {cache_path}")


if __name__ == "__main__":
    main()
