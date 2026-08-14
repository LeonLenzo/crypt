#!/usr/bin/env python3
"""
03b_fetch_literature.py — fetch PubMed / PMC literature for all BioProjects in runs.tsv.

Reads BioProject accessions from runs.tsv and runs up to six PMID strategies
(short-circuit, stop at first hit) per BioProject:
  1. BioProject XML <Publication> elements
  2. PMC full-text search  (40% yield)
  3. Europe PMC            (35% yield)
  4. ENA XML               (PRJEB accessions only)
  5. Semantic Scholar      (S2_API_KEY required)
  6. PubMed text search    (last resort)

After finding a PMID, fetches the PMC full-text methods section if available.

Cache: output/03b_fetch_literature/data/lit_cache.json (existing v5 entries are
instant skips; use --retry to re-run no-PMID entries).

Output:
  output/03b_fetch_literature/data/literature.json   one entry per BioProject
  output/03b_fetch_literature/logs/

Usage:
  python 03b_fetch_literature.py           # upgrade v4→v5 + fetch new BPs
  python 03b_fetch_literature.py --retry   # also retry no-PMID entries
"""

import argparse
import csv
import re
from datetime import datetime
import json
import os
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, http_get, link_latest, load_json, make_log_dir, save_json

RUNS_TSV      = Path("stat/output/filter_runs/data/runs.tsv")
CACHE_PATH    = Path("metadata/output/fetch_lit/data/lit_cache.json")
OUT_DIR       = Path("metadata/output/fetch_lit")

CACHE_VERSION = 5
METHODS_MAX_CHARS = 8000

API_KEY       = os.environ.get("NCBI_API_KEY", "")
S2_API_KEY    = os.environ.get("S2_API_KEY", "")
SERPER_API_KEY = os.environ.get("SERPER_API_KEY", "")
RATE       = 9.0 if API_KEY else 2.5
HEADERS    = {"User-Agent": "crypt/03b_fetch_literature (leon.lenzo@curtin.edu.au)"}
ENTREZ     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

EPMC        = "https://www.ebi.ac.uk/europepmc/webservices/rest"
ENA         = "https://www.ebi.ac.uk/ena/browser/api"
CROSSREF    = "https://api.crossref.org/works"
EXT_HEADERS = {"User-Agent": "crypt/03b_fetch_literature (leon.lenzo@curtin.edu.au)"}
CROSSREF_HEADERS = {"User-Agent": "crypt/fetch_lit (mailto:leon.lenzo@curtin.edu.au)"}

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
    gap  = 0.5
    wait = _last_ext + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_ext = time.monotonic()


_last_s2: float = 0.0


def _s2_wait() -> None:
    global _last_s2
    gap  = 1.1
    wait = _last_s2 + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_s2 = time.monotonic()


_last_serper: float = 0.0


def _serper_wait() -> None:
    global _last_serper
    gap  = 1.0
    wait = _last_serper + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_serper = time.monotonic()


_last_crossref: float = 0.0


def _crossref_wait() -> None:
    global _last_crossref
    gap  = 1.0
    wait = _last_crossref + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_crossref = time.monotonic()


_DOI_RE  = re.compile(r'10\.\d{4,}/[^\s"\'<>]+')
_PMID_RE = re.compile(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d{6,9})')

_ACADEMIC_DOMAINS = {
    "pubmed.ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov",
    "biorxiv.org", "medrxiv.org", "preprints.org",
    "frontiersin.org", "mdpi.com", "journals.plos.org", "elifesciences.org",
    "peerj.com", "f1000research.com", "nature.com", "science.org",
    "onlinelibrary.wiley.com", "academic.oup.com", "link.springer.com",
    "cell.com", "sciencedirect.com", "tandfonline.com",
    "royalsocietypublishing.org", "cambridge.org", "apsnet.org", "mpmi.org",
    "plantcell.org", "plantphysiol.org", "journals.asm.org", "phytobiomes.org",
    "europepmc.org", "semanticscholar.org", "doi.org",
    "researchsquare.com", "researchgate.net",
}


def _clean_doi(doi: str) -> str:
    """Strip URL artifacts from extracted DOI strings."""
    doi = doi.split("?")[0]
    doi = doi.rstrip(".,;)")
    if "/attachment/" in doi:
        doi = doi.split("/attachment/")[0]
    if "~" in doi:                              # URL slug separator (some CMS systems)
        doi = doi.split("~")[0]
    parts = doi.split("/")
    if len(parts) >= 4 and parts[3].isdigit():
        doi = "/".join(parts[:3])
    for sfx in ("/full.pdf", ".full.pdf", "/full", ".pdf"):
        if doi.endswith(sfx):
            doi = doi[: -len(sfx)]
    return doi.rstrip(".,;)")


def _doi_to_pmid(doi: str) -> str:
    try:
        raw   = _get("esearch.fcgi", db="pubmed", term=f"{doi}[doi]", retmax=1)
        root  = ET.fromstring(raw)
        count = int(root.findtext(".//Count") or "0")
        if count != 1:
            return ""
        ids = [el.text for el in root.findall(".//Id") if el.text]
        return ids[0] if ids else ""
    except Exception:
        return ""


def _europepmc_doi_to_pmid(doi: str) -> str:
    """Look up a DOI in EuropePMC; return PMID when available (broader than PubMed)."""
    query = urllib.parse.quote(f'DOI:"{doi}"')
    url   = f"{EPMC}/search?query={query}&format=json&pageSize=3&fields=pmid"
    _ext_wait()
    try:
        raw  = http_get(url, EXT_HEADERS, no_retry_429=True)
        data = json.loads(raw)
    except Exception:
        return ""
    for r in data.get("resultList", {}).get("result", []):
        pmid = (r.get("pmid") or "").strip()
        if pmid:
            return pmid
    return ""


def _crossref_get(doi: str) -> dict:
    """Fetch title, abstract, and pub_date from CrossRef for any DOI."""
    url = f"{CROSSREF}/{urllib.parse.quote(doi, safe='')}"
    _crossref_wait()
    try:
        raw  = http_get(url, CROSSREF_HEADERS)
        work = json.loads(raw).get("message", {})
    except Exception:
        return {}
    titles   = work.get("title", [])
    title    = titles[0] if titles else ""
    abstract = re.sub(r"<[^>]+>", "", work.get("abstract", "")).strip()
    dp       = ((work.get("published") or work.get("published-print") or {})
                .get("date-parts", [[]])) or [[]]
    parts    = dp[0] if dp else []
    pub_date = ("-".join([str(parts[0])]
                         + [str(x).zfill(2) for x in parts[1:3]])
                if parts else "")
    return {"title": title, "abstract": abstract, "pub_date": pub_date}


# ── External search functions ──────────────────────────────────────────────────

def _europepmc_search(accession: str) -> list[str]:
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
    return [r["pmid"] for r in data.get("resultList", {}).get("result", [])
            if r.get("pmid", "").strip()]


def _ena_pmids(accession: str) -> list[str]:
    if not accession.startswith("PRJEB"):
        return []
    _ext_wait()
    try:
        raw  = http_get(f"{ENA}/xml/{accession}", EXT_HEADERS, no_retry_429=True)
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


def _serper_search(accession: str) -> list[str]:
    if not SERPER_API_KEY:
        return []
    _serper_wait()
    body = json.dumps({"q": f'"{accession}"', "num": 10}).encode()
    req  = urllib.request.Request(
        "https://google.serper.dev/search", data=body,
        headers={"X-API-KEY": SERPER_API_KEY, "Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            data = json.loads(r.read())
    except Exception:
        return []

    pmids: list[str] = []
    for it in data.get("organic", []):
        link    = it.get("link", "")
        snippet = it.get("snippet", "")
        text    = link + " " + snippet

        # Direct PubMed URL hit
        for pmid in _PMID_RE.findall(text):
            if pmid not in pmids:
                pmids.append(pmid)

        # DOI from an academic domain → resolve to PMID
        if not pmids:
            from urllib.parse import urlparse
            domain = urlparse(link).netloc.lstrip("www.")
            if any(domain == d or domain.endswith("." + d) for d in _ACADEMIC_DOMAINS):
                for doi in _DOI_RE.findall(text):
                    doi = doi.rstrip(".,;)")
                    pmid = _doi_to_pmid(doi)
                    if pmid and pmid not in pmids:
                        pmids.append(pmid)
                        break

    return pmids


def _semantic_scholar_search(accession: str) -> list[str]:
    if not S2_API_KEY:
        return []
    params  = urllib.parse.urlencode({"query": accession, "fields": "externalIds", "limit": 5})
    headers = {**EXT_HEADERS, "x-api-key": S2_API_KEY}
    _s2_wait()
    try:
        raw  = http_get(
            f"https://api.semanticscholar.org/graph/v1/paper/search?{params}",
            headers, no_retry_429=True)
        data = json.loads(raw)
    except Exception:
        return []
    return [str((p.get("externalIds") or {}).get("PubMed", "")).strip()
            for p in data.get("data", [])
            if str((p.get("externalIds") or {}).get("PubMed", "")).strip().isdigit()]


# ── XML helpers ────────────────────────────────────────────────────────────────

def _xtext(root: ET.Element, xpath: str) -> str:
    el = root.find(xpath)
    return (el.text or "").strip() if el is not None else ""


_MONTH_MAP = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
              "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
              "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}


def _pub_date(article: ET.Element) -> str:
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
        dy  = (pd.findtext("Day") or "01").zfill(2) if (pd.findtext("Day") or "01").isdigit() else "01"
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
    raw  = _get("esearch.fcgi", db="pmc", term=f"{pmid}[pmid]", retmax=1)
    root = ET.fromstring(raw)
    ids  = [el.text for el in root.findall(".//Id") if el.text]
    return ids[0] if ids else None


def _fetch_pmc_methods(pmcid: str) -> str:
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
            return " ".join(p for p in parts if p)[:METHODS_MAX_CHARS]
    return ""


_STRATEGIES = [
    ("BioProject XML",   lambda acc, uid, xml_pmids: list(xml_pmids)),
    ("PMC full-text",    lambda acc, uid, _: _pmcids_to_pmids(_pmc_search(acc))),
    ("Europe PMC",       lambda acc, uid, _: _europepmc_search(acc)),
    ("ENA XML",          lambda acc, uid, _: _ena_pmids(acc)),
    ("Semantic Scholar", lambda acc, uid, _: _semantic_scholar_search(acc)),
    ("PubMed",           lambda acc, uid, _: _pubmed_search(acc)),
    ("Serper",           lambda acc, uid, _: _serper_search(acc)),
]


# ── Core fetch function ────────────────────────────────────────────────────────

def fetch_bioproject_meta(accession: str, cache: dict) -> dict:
    if accession in cache and cache[accession].get("_v") == CACHE_VERSION:
        if cache[accession].get("n_papers_found", 0) > 0:
            result = dict(cache[accession])
            result.setdefault("pmid_source", "cached")
            return result

    if accession in cache and cache[accession].get("_v") == 4:
        entry = cache[accession]
        if entry.get("primary_pmid"):
            entry.setdefault("pmcid", "")
            entry.setdefault("methods_text", "")
            if not entry.get("abstract", "").strip():
                try:
                    raw  = _get("efetch.fcgi", db="pubmed",
                                id=entry["primary_pmid"], rettype="xml", retmode="xml")
                    root = ET.fromstring(raw)
                    for article in root.findall(".//PubmedArticle"):
                        abs_els  = article.findall(
                            ".//MedlineCitation/Article/Abstract/AbstractText")
                        abstract = " ".join(
                            "".join(el.itertext()).strip() for el in abs_els).strip()
                        ttl_el   = article.find(".//MedlineCitation/Article/ArticleTitle")
                        title    = ("".join(ttl_el.itertext()).strip()
                                    if ttl_el is not None else entry.get("title", ""))
                        if abstract:
                            entry["abstract"]         = abstract
                            entry["title"]            = title
                            entry["primary_pub_text"] = f"{title} {abstract}"
                        break
                except Exception as e:
                    print(f"WARNING: PubMed re-fetch failed for {accession}: {e}", flush=True)
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

    result: dict = {
        "title": "", "description": "", "bp_submission_date": "",
        "primary_pmid": "", "primary_doi": "", "primary_pub_date": "",
        "primary_publication": "", "abstract": "",
        "n_papers_found": 0, "primary_pub_text": "", "pmid_source": "none",
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

    bp_xml_pmids: set[str] = set()
    bp_xml_dois:  set[str] = set()
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
            else:
                m = _DOI_RE.search(pid)
                if m:
                    bp_xml_dois.add(m.group().rstrip(".,;)"))
    except Exception as e:
        print(f"WARNING: efetch failed for {accession}: {e}", flush=True)

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
                pid      = pid_el.text or ""
                abstract = " ".join(
                    "".join(el.itertext()).strip() for el in abs_els).strip()
                title    = "".join(ttl_el.itertext()).strip() if ttl_el is not None else ""
                doi      = ""
                for aid in article.findall(".//ArticleIdList/ArticleId"):
                    if aid.get("IdType") == "doi":
                        doi = (aid.text or "").strip()
                        break
                pub_info[pid] = {"title": title, "date": _pub_date(article),
                                 "abstract": abstract, "doi": doi}
            if pub_info:
                primary_pid = min(pub_info, key=lambda p: pub_info[p]["date"])
                pri = pub_info[primary_pid]
                result["primary_pmid"]        = primary_pid
                result["primary_doi"]         = pri.get("doi", "")
                result["primary_pub_date"]    = pri["date"]
                result["primary_publication"] = f"[{primary_pid}] {pri['title']}"
                result["abstract"]            = pri["abstract"]
                result["primary_pub_text"]    = f"{pri['title']} {pri['abstract']}"
        except Exception as e:
            print(f"WARNING: pubmed efetch failed for {accession}: {e}", flush=True)

    # ── DOI-path fallback: no PMID found via any strategy ──────────────────────
    if not result["primary_pmid"]:
        # Collect candidate DOIs: BioProject XML first, then Serper-style URL scan
        candidate_dois = list(bp_xml_dois)
        for doi in candidate_dois:
            pmid = _doi_to_pmid(doi)
            if not pmid:
                pmid = _europepmc_doi_to_pmid(doi)
            if pmid:
                # Resolved to a PMID — fetch PubMed metadata normally
                try:
                    raw  = _get("efetch.fcgi", db="pubmed",
                                id=pmid, rettype="xml", retmode="xml")
                    root = ET.fromstring(raw)
                    for article in root.findall(".//PubmedArticle"):
                        pid_el  = article.find(".//MedlineCitation/PMID")
                        ttl_el  = article.find(".//MedlineCitation/Article/ArticleTitle")
                        abs_els = article.findall(".//MedlineCitation/Article/Abstract/AbstractText")
                        if pid_el is None:
                            continue
                        abstract = " ".join(
                            "".join(el.itertext()).strip() for el in abs_els).strip()
                        title    = "".join(ttl_el.itertext()).strip() if ttl_el is not None else ""
                        result["primary_pmid"]        = pmid
                        result["primary_doi"]         = doi
                        result["primary_pub_date"]    = _pub_date(article)
                        result["primary_publication"] = f"[{pmid}] {title}"
                        result["abstract"]            = abstract
                        result["primary_pub_text"]    = f"{title} {abstract}"
                        result["n_papers_found"]      = 1
                        result["pmid_source"]         = "BioProject XML DOI"
                        break
                except Exception:
                    pass
                break
            else:
                # Keep the DOI even without a PMID; enrich from CrossRef
                result["primary_doi"] = doi
                cr = _crossref_get(doi)
                if cr.get("title"):
                    result["primary_publication"] = f"[DOI:{doi}] {cr['title']}"
                    result["abstract"]            = cr.get("abstract", "")
                    result["primary_pub_date"]    = cr.get("pub_date", "")
                    result["primary_pub_text"]    = f"{cr['title']} {cr.get('abstract','')}"
                    result["n_papers_found"]      = 1
                    result["pmid_source"]         = "BioProject XML DOI (CrossRef)"

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

_LIT_FIELDS = [
    "primary_pmid", "primary_doi", "primary_pub_date", "primary_publication",
    "abstract", "methods_text", "pmcid", "n_papers_found", "pmid_source",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--retry", action="store_true",
                    help="re-run PMID strategies for BPs with no PMID found")
    args = ap.parse_args()

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    sys.stdout = _Tee(log_dir / "fetch.log")
    link_latest(logs_base, log_dir / "fetch.log")

    print(f"NCBI_API_KEY:    {'yes' if API_KEY else 'no'}", flush=True)
    print(f"S2_API_KEY:      {'yes' if S2_API_KEY else 'no'}", flush=True)
    print(f"SERPER_API_KEY:  {'yes' if SERPER_API_KEY else 'no'}", flush=True)
    print(f"Cache:        {CACHE_PATH}", flush=True)

    if not RUNS_TSV.exists():
        sys.exit(f"ERROR: {RUNS_TSV} not found — run 02_filter_runs.py first")

    with open(RUNS_TSV, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    seen: set[str] = set()
    bioprojects: list[str] = []
    for r in rows:
        bp = r.get("BioProject", "").strip()
        if bp and bp not in ("", "NA") and bp not in seen:
            seen.add(bp)
            bioprojects.append(bp)

    cache = load_json(CACHE_PATH) if CACHE_PATH.exists() else {}

    to_upgrade: list[str] = []
    to_fetch:   list[str] = []
    to_retry:   list[str] = []
    n_skip_pmid = n_skip_no_pmid = 0

    for bp in bioprojects:
        entry    = cache.get(bp, {})
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
                entry["_v"] = CACHE_VERSION
                entry.setdefault("pmcid", "")
                entry.setdefault("methods_text", "")
                cache[bp] = entry
                n_skip_no_pmid += 1
        else:
            to_fetch.append(bp)

    total_work = len(to_upgrade) + len(to_fetch) + len(to_retry)
    print(f"BioProjects: {len(bioprojects)}", flush=True)
    print(f"  skip — v5 with PMID:  {n_skip_pmid}", flush=True)
    print(f"  skip — no PMID found: {n_skip_no_pmid}", flush=True)
    print(f"  upgrade v4→v5:        {len(to_upgrade)}", flush=True)
    print(f"  fetch (new):          {len(to_fetch)}", flush=True)
    if args.retry:
        print(f"  retry (no-PMID):      {len(to_retry)}", flush=True)
    print(f"  → work items:         {total_work}", flush=True)

    if n_skip_no_pmid and not args.retry:
        save_json(cache, CACHE_PATH)

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
        save_json(cache, CACHE_PATH)
        if kind == "upgrade":
            n_upgraded += 1
        elif kind == "retry":
            n_retried += 1
        else:
            n_fetched += 1

    save_json(cache, CACHE_PATH)

    # ── Export literature.json ────────────────────────────────────────────────
    lit_out: dict[str, dict] = {}
    for bp in bioprojects:
        entry = cache.get(bp, {})
        lit_out[bp] = {k: entry.get(k, "") for k in _LIT_FIELDS}

    lit_path = OUT_DIR / "data" / "literature.json"
    lit_path.write_text(json.dumps(lit_out, indent=2))
    print(f"Written: {lit_path}  ({len(lit_out):,} entries)", flush=True)

    # ── Summary ───────────────────────────────────────────────────────────────
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
        f"03b_fetch_literature summary\n"
        f"BioProjects total:       {len(bioprojects)}\n"
        f"  skipped (v5+PMID):     {n_skip_pmid}\n"
        f"  skipped (no PMID):     {n_skip_no_pmid}\n"
        f"  upgraded v4→v5:        {n_upgraded}\n"
        f"  fetched (new):         {n_fetched}\n"
        f"  retried (no-PMID):     {n_retried}\n"
        f"With PMID:               {n_with_pmid}\n"
        f"With PMC methods text:   {n_with_methods}\n"
        f"PMID source breakdown:\n{source_lines}\n"
        f"\nOutput: {lit_path}\n"
        f"Run 04_filter_kw.py → 05_llm_classify.py for study design classification.\n"
    )
    summary_path = log_dir / "fetch_summary.txt"
    summary_path.write_text(summary)
    link_latest(logs_base, summary_path)
    print(f"\n{summary}")


if __name__ == "__main__":
    main()
