#!/usr/bin/env python3
"""
metadata/meta_search.py — resolve BioProject identity (title/DOI/PMID/PMCID),
nothing else. Replaces ncbi_metadata.py + web_metadata.py.

Job: for every BioProject in runs.tsv, answer "what paper is this?" as cheaply
as possible, short-circuiting at the first hit. Never touches PDFs or full
text — that is meta_text.py's job entirely (see meta_text.py --fetch).

Single-pass cascade per BioProject (each stage only runs if the previous
found nothing):
  1. bioproject_xml     — <Publication> elements in BioProject XML (free)
  2. pmc_search          — PMC accession full-text search
  3. pubmed_search        — PubMed text search
  4. serper               — Google search (accession -> title -> keywords),
                            DOI/PMID extraction from links+snippets, page
                            scraping of academic domains for citation meta
                            tags, PubMed title search as last resort
  5. crossref             — DOI-only entries (from any stage above) get
                            title/abstract/pub_date filled from CrossRef

BioSample XML + ENA/DDBJ BioSample attributes are fetched in the same pass
(unrelated to literature resolution, but same NCBI/EBI round-trip).

Outputs:
  metadata/output/meta_search/data/bioprojects.json  — title, description,
    submission_date, pmid, pmcid, doi, pub_date, publication, abstract,
    pmid_source. No full text, no methods_text — that lives in meta_text.py's
    text_cache.jsonl, keyed by doi.
  metadata/output/meta_search/data/biosamples.json   — BioSample XML/ENA attrs

Cache:
  metadata/output/meta_search/data/bp_cache.json      — per-BP, resumable
  metadata/output/meta_search/data/bs_cache.json       — per-BioSample, resumable
  metadata/output/meta_search/data/serper_cache.json   — raw Serper results
                                                          (expensive; kept forever)

Flags:
  --migrate-legacy   one-time: merge metadata/output/ncbi_metadata/data/bp_cache.json
                     + metadata/output/web_metadata/data/bp_cache.json into the new
                     cache (methods_text dropped), so previously-paid-for NCBI/Serper
                     lookups aren't repeated. Run this once before the first normal run.
  --retry            re-attempt BPs with no PMID/DOI found on a previous run
  --bootstrap        import NCBI-sourced entries from legacy fetch_lit lit_cache.json
  --bootstrap-serper import legacy serper_results.tsv into serper_cache.json
  --fetch-ena        fetch ENA/DDBJ BioSample attrs via EBI API only, skip everything else
  --report           print resolution breakdown and exit (no fetching)
  --dry-run          resolve without writing bp_cache (serper_cache still written)
  --limit N          process at most N BPs (0 = all; for testing)

Requires: SERPER_API_KEY env var for stage 4 (missing key just skips that stage
with a warning — stages 1-3 and 5 still run); NCBI_API_KEY strongly recommended.

Run from crypt/:
  python metadata/meta_search.py --migrate-legacy   # once: import old caches
  python metadata/meta_search.py                    # fetch remaining / new BPs
  python metadata/meta_search.py --retry             # re-attempt unresolved BPs
  python metadata/meta_search.py --fetch-ena         # fill in ENA/DDBJ BioSample attrs
  python metadata/meta_search.py --report            # show resolution breakdown
"""

import argparse
import collections
import csv
import html.parser
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, http_get, link_latest, load_json, make_log_dir, save_json

# ── Paths ──────────────────────────────────────────────────────────────────────

RUNS_TSV = Path("stat/output/stat_filter/data/runs.tsv")
OUT_DIR  = Path("metadata/output/meta_search")
BP_CACHE     = OUT_DIR / "data" / "bp_cache.json"
BS_CACHE     = OUT_DIR / "data" / "bs_cache.json"
SERPER_CACHE = OUT_DIR / "data" / "serper_cache.json"

LEGACY_LIT_CACHE    = Path("metadata/output/fetch_lit/data/lit_cache.json")
LEGACY_SERPER_TSV   = Path("metadata/output/serper/serper_results.tsv")
# ncbi_metadata.py / web_metadata.py are retired (superseded by this script,
# 2026-08-25); their outputs were moved to metadata/legacy/output/ — kept here
# only so --migrate-legacy remains re-runnable if bp_cache.json is ever lost.
LEGACY_NCBI_CACHE   = Path("metadata/legacy/output/ncbi_metadata/data/bp_cache.json")
LEGACY_NCBI_BS      = Path("metadata/legacy/output/ncbi_metadata/data/bs_cache.json")
LEGACY_WEB_CACHE    = Path("metadata/legacy/output/web_metadata/data/bp_cache.json")
LEGACY_WEB_SERPER   = Path("metadata/legacy/output/web_metadata/data/serper_cache.json")

CACHE_VERSION = 1
BS_BATCH      = 300

# ── Credentials / rates ────────────────────────────────────────────────────────

NCBI_KEY   = os.environ.get("NCBI_API_KEY", "")
SERPER_KEY = os.environ.get("SERPER_API_KEY", "")
NCBI_RATE  = 9.0 if NCBI_KEY else 2.5

ENTREZ   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest"
CROSSREF = "https://api.crossref.org/works"

_UA = "crypt/meta_search (leon.lenzo@curtin.edu.au)"
_H_NCBI = {"User-Agent": _UA}
_H_EXT  = {"User-Agent": _UA}
_H_CR   = {"User-Agent": f"crypt/meta_search (mailto:leon.lenzo@curtin.edu.au)"}
_H_PAGE = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0",
    "Accept": "text/html,application/xhtml+xml",
    "Accept-Language": "en-US,en;q=0.9",
}

# ENA / EBI BioSamples API — covers SAME* (ENA) and SAMD* (DDBJ) accessions
_ENA_BASE     = "https://www.ebi.ac.uk/biosamples/samples"
_ENA_HEADERS  = {"Accept": "application/json", "User-Agent": _UA}
_ENA_ATTR_MAP = {
    "geographic location (country and/or sea)": "geo_loc_name",
    "collection date":                           "collection_date",
    "organism part":                             "tissue",
    "host scientific name":                      "host",
    "lat_lon":                                   "lat_lon",
    "description":                               "bs_description",
}
_ENA_LOCALITY = "geographic location (region and locality)"
_ENA_PREFIXES = {"SAME", "SAMD"}
_HARMONIZED = {
    "tissue", "geo_loc_name", "collection_date", "lat_lon",
    "host", "isolation_source", "dev_stage",
}

_MONTH_MAP = {
    "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
    "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
    "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12",
}

_DOI_RE    = re.compile(r'10\.\d{4,}/[^\s"\'<>]+')
_PMID_RE   = re.compile(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d{6,9})')
_PREPRT_RE = re.compile(r'(?:biorxiv|medrxiv)\.org/content/[^\s]*?(\d{4}\.\d{2}\.\d{2}\.\d+)')

_ACADEMIC = {
    "pubmed.ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov",
    "biorxiv.org", "medrxiv.org", "preprints.org",
    "frontiersin.org", "mdpi.com", "journals.plos.org", "elifesciences.org",
    "peerj.com", "f1000research.com", "nature.com", "science.org",
    "onlinelibrary.wiley.com", "academic.oup.com", "link.springer.com",
    "cell.com", "sciencedirect.com", "tandfonline.com",
    "royalsocietypublishing.org", "cambridge.org", "apsnet.org", "mpmi.org",
    "plantcell.org", "plantphysiol.org", "journals.asm.org", "phytobiomes.org",
    "europepmc.org", "doi.org", "researchsquare.com",
}
_SCRAPABLE = {
    "mdpi.com", "frontiersin.org", "journals.plos.org", "elifesciences.org",
    "peerj.com", "f1000research.com", "preprints.org", "researchsquare.com",
    "biorxiv.org", "medrxiv.org", "apsjournals.apsnet.org", "phytobiomes.org",
    "apsnet.org", "mpmi.org", "plantcell.org", "plantphysiol.org",
    "journals.asm.org", "academic.oup.com", "europepmc.org", "doi.org",
    "pubs.acs.org",
}
_NO_SCRAPE = {
    "researchgate.net", "omicsdi.org", "seqout.org",
    "agdatacommons.nal.usda.gov", "gold.jgi.doe.gov", "sciencedirect.com",
}

# NCBI-only sources — entries with these pmid_source values are eligible for --bootstrap
# ── Rate limiting ──────────────────────────────────────────────────────────────

_last: dict[str, float] = {
    "ncbi": 0.0, "ext": 0.0, "cr": 0.0, "serper": 0.0, "page": 0.0, "ena": 0.0,
}


def _wait(key: str, rate: float) -> None:
    gap = 1.0 / rate
    w   = _last[key] + gap - time.monotonic()
    if w > 0:
        time.sleep(w)
    _last[key] = time.monotonic()


def _ncbi_get(path: str, **params) -> bytes:
    if NCBI_KEY:
        params["api_key"] = NCBI_KEY
    _wait("ncbi", NCBI_RATE)
    return http_get(f"{ENTREZ}/{path}?{urllib.parse.urlencode(params)}", _H_NCBI)


def _ncbi_post(path: str, **params) -> bytes:
    if NCBI_KEY:
        params["api_key"] = NCBI_KEY
    _wait("ncbi", NCBI_RATE)
    url  = f"{ENTREZ}/{path}"
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=_H_NCBI, method="POST")
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def _ext_get(url: str) -> bytes:
    _wait("ext", 2.0)
    with urllib.request.urlopen(urllib.request.Request(url, headers=_H_EXT), timeout=20) as r:
        return r.read()


def _cr_get(doi: str) -> bytes:
    _wait("cr", 1.0)
    url = f"{CROSSREF}/{urllib.parse.quote(doi, safe='')}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=_H_CR), timeout=20) as r:
        return r.read()


def _page_get(url: str) -> bytes:
    _wait("page", 1.0)
    req = urllib.request.Request(url, headers=_H_PAGE)
    with urllib.request.urlopen(req, timeout=20) as r:
        return r.read(32768)


def _ena_wait() -> None:
    _wait("ena", 2.0)


def _fetch_ena_biosample(accession: str) -> dict:
    _ena_wait()
    try:
        data = json.loads(http_get(f"{_ENA_BASE}/{accession}", _ENA_HEADERS))
    except Exception:
        return {"accession": accession}
    chars = data.get("characteristics", {})
    row: dict = {"accession": accession}
    for ena_key, our_key in _ENA_ATTR_MAP.items():
        items = chars.get(ena_key, [])
        val   = items[0].get("text", "").strip() if items else ""
        if val:
            row[our_key] = val
    loc_items = chars.get(_ENA_LOCALITY, [])
    locality  = loc_items[0].get("text", "").strip() if loc_items else ""
    if locality:
        row["geo_loc_name"] = (
            f"{row['geo_loc_name']}: {locality}" if "geo_loc_name" in row else locality
        )
    return row


# ── BioProject XML ─────────────────────────────────────────────────────────────

def _xtext(root: ET.Element, xpath: str) -> str:
    el = root.find(xpath)
    return (el.text or "").strip() if el is not None else ""


def _bp_uid(accession: str) -> str | None:
    raw  = _ncbi_get("esearch.fcgi", db="bioproject",
                     term=f"{accession}[Project Accession]", retmax=1)
    ids  = [el.text for el in ET.fromstring(raw).findall(".//Id") if el.text]
    return ids[0] if ids else None


def _fetch_bp_xml(uid: str) -> tuple[dict, set[str], set[str]]:
    """Fetch BioProject XML; return (meta_dict, xml_pmids, xml_dois)."""
    raw  = _ncbi_get("efetch.fcgi", db="bioproject", id=uid)
    root = ET.fromstring(raw)
    meta = {
        "title":           _xtext(root, ".//ProjectDescr/Title"),
        "description":     _xtext(root, ".//ProjectDescr/Description"),
        "submission_date": "",
    }
    sub_el = root.find(".//Submission[@submitted]")
    if sub_el is not None:
        meta["submission_date"] = sub_el.get("submitted", "")

    xml_pmids: set[str] = set()
    xml_dois:  set[str] = set()
    for pub in root.findall(".//Publication"):
        pid = (pub.get("id") or "").strip()
        if pid.isdigit():
            xml_pmids.add(pid)
        else:
            m = _DOI_RE.search(pid)
            if m:
                xml_dois.add(m.group().rstrip(".,;)"))

    return meta, xml_pmids, xml_dois


# ── NCBI PMID search strategies ────────────────────────────────────────────────

def _pmc_search(accession: str) -> list[str]:
    raw    = _ncbi_get("esearch.fcgi", db="pmc", term=f'"{accession}"', retmax=50)
    pmcids = [el.text for el in ET.fromstring(raw).findall(".//Id") if el.text]
    if not pmcids:
        return []
    raw  = _ncbi_get("esummary.fcgi", db="pmc", id=",".join(pmcids), retmode="json")
    data = json.loads(raw)
    pmids = []
    for pmcid in pmcids:
        for aid in data.get("result", {}).get(pmcid, {}).get("articleids", []):
            if aid.get("idtype") == "pmid" and aid.get("value"):
                pmids.append(str(aid["value"]))
                break
    return pmids


def _pubmed_search(accession: str) -> list[str]:
    raw = _ncbi_get("esearch.fcgi", db="pubmed", term=f'"{accession}"', retmax=50)
    return [el.text for el in ET.fromstring(raw).findall(".//Id") if el.text]


def _title_to_pmid(title: str) -> str:
    if len(title.split()) < 6:
        return ""
    query = " ".join(title.split()[:12]) + "[ti]"
    try:
        raw  = _ncbi_get("esearch.fcgi", db="pubmed", term=query, retmax=2)
        root = ET.fromstring(raw)
        if int(root.findtext(".//Count") or "0") != 1:
            return ""
        ids = [el.text for el in root.findall(".//Id") if el.text]
        return ids[0] if ids else ""
    except Exception:
        return ""


def _doi_to_pmid_pubmed(doi: str) -> str:
    try:
        raw  = _ncbi_get("esearch.fcgi", db="pubmed", term=f"{doi}[doi]", retmax=1)
        root = ET.fromstring(raw)
        if int(root.findtext(".//Count") or "0") != 1:
            return ""
        ids = [el.text for el in root.findall(".//Id") if el.text]
        return ids[0] if ids else ""
    except Exception:
        return ""


def _doi_to_pmid_europepmc(doi: str) -> str:
    query = urllib.parse.quote(f'DOI:"{doi}"')
    url   = f"{EPMC_URL}/search?query={query}&format=json&pageSize=3&fields=pmid"
    try:
        for r in json.loads(_ext_get(url)).get("resultList", {}).get("result", []):
            pmid = (r.get("pmid") or "").strip()
            if pmid:
                return pmid
    except Exception:
        pass
    return ""


def _pmid_to_pmcid(pmid: str) -> str:
    try:
        raw = _ncbi_get("esearch.fcgi", db="pmc", term=f"{pmid}[pmid]", retmax=1)
        ids = [el.text for el in ET.fromstring(raw).findall(".//Id") if el.text]
        return ids[0] if ids else ""
    except Exception:
        return ""


def _doi_to_pmcid(doi: str) -> str:
    """Search PMC by DOI — catches OA papers not indexed in PubMed."""
    try:
        raw = _ncbi_get("esearch.fcgi", db="pmc", term=f"{doi}[doi]", retmax=1)
        ids = [el.text for el in ET.fromstring(raw).findall(".//Id") if el.text]
        return ids[0] if ids else ""
    except Exception:
        return ""


def _fetch_pubmed_meta(pmids: list[str]) -> dict[str, dict]:
    """Batch-fetch title, abstract, DOI, pub_date for a list of PMIDs."""
    if not pmids:
        return {}
    raw  = _ncbi_get("efetch.fcgi", db="pubmed",
                     id=",".join(sorted(pmids)), rettype="xml", retmode="xml")
    root = ET.fromstring(raw)
    out: dict[str, dict] = {}
    for article in root.findall(".//PubmedArticle"):
        pid_el  = article.find(".//MedlineCitation/PMID")
        ttl_el  = article.find(".//MedlineCitation/Article/ArticleTitle")
        abs_els = article.findall(".//MedlineCitation/Article/Abstract/AbstractText")
        if pid_el is None:
            continue
        pid      = pid_el.text or ""
        title    = "".join(ttl_el.itertext()).strip() if ttl_el is not None else ""
        abstract = " ".join("".join(el.itertext()).strip() for el in abs_els).strip()
        doi      = ""
        for aid in article.findall(".//ArticleIdList/ArticleId"):
            if aid.get("IdType") == "doi":
                doi = (aid.text or "").strip()
                break
        out[pid] = {"title": title, "abstract": abstract, "doi": doi, "date": _pub_date(article)}
    return out


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
        dy  = (pd.findtext("Day") or "01").zfill(2)
        if yr and yr.isdigit():
            return f"{yr}-{mo}-{dy}"
        ml = pd.findtext("MedlineDate") or ""
        if ml[:4].isdigit():
            return f"{ml[:4]}-01-01"
    return "9999-12-31"


def _crossref_get(doi: str) -> dict:
    try:
        work = json.loads(_cr_get(doi)).get("message", {})
    except Exception:
        return {}
    titles   = work.get("title", [])
    abstract = re.sub(r"<[^>]+>", "", work.get("abstract", "")).strip()
    dp       = ((work.get("published") or work.get("published-print") or {})
                .get("date-parts", [[]])) or [[]]
    parts    = dp[0] if dp else []
    pub_date = ("-".join([str(parts[0])] + [str(x).zfill(2) for x in parts[1:3]])
                if parts else "")
    return {"title": titles[0] if titles else "", "abstract": abstract, "pub_date": pub_date}


def _fetch_pubmed_meta_one(pmid: str) -> dict:
    info = _fetch_pubmed_meta([pmid])
    return info.get(pmid, {})


# ── Serper (Google) search ─────────────────────────────────────────────────────

def _serper_query(q: str, num: int = 10) -> list[dict]:
    _wait("serper", 1.0)
    body = json.dumps({"q": q, "num": num}).encode()
    req  = urllib.request.Request(
        "https://google.serper.dev/search", data=body,
        headers={"X-API-KEY": SERPER_KEY, "Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            return json.loads(r.read()).get("organic", [])
    except Exception as e:
        print(f"    serper error: {e}", flush=True)
        return []


def _clean_doi(doi: str) -> str:
    doi = doi.split("?")[0].rstrip(".,;)")
    if "/attachment/" in doi:
        doi = doi.split("/attachment/")[0]
    if "~" in doi:
        doi = doi.split("~")[0]
    parts = doi.split("/")
    if len(parts) >= 4 and parts[3].isdigit():
        doi = "/".join(parts[:3])
    for sfx in ("/full.pdf", ".full.pdf", "/full", ".pdf", "/pdf"):
        if doi.endswith(sfx):
            doi = doi[: -len(sfx)]
    # bioRxiv/medRxiv (10.1101/...) DOIs never include the version suffix
    # (v1, v2, ...) that appears in their URLs — e.g. the canonical DOI for
    # https://www.biorxiv.org/content/10.1101/350785v1 is 10.1101/350785.
    # Scraped text sometimes carries the v-suffix straight into the "DOI".
    if doi.startswith("10.1101/"):
        doi = re.sub(r"v\d+$", "", doi)
    return doi.rstrip(".,;)")


def _domain(url: str) -> str:
    return urllib.parse.urlparse(url).netloc.lstrip("www.")


def _is_academic(url: str) -> bool:
    d = _domain(url)
    return any(d == a or d.endswith("." + a) for a in _ACADEMIC)


def _is_scrapable(url: str) -> bool:
    d = _domain(url)
    if any(d == s or d.endswith("." + s) for s in _NO_SCRAPE):
        return False
    return any(d == s or d.endswith("." + s) for s in _SCRAPABLE)


def _extract_pmids(results: list[dict]) -> list[str]:
    pmids: list[str] = []
    for r in results:
        text = r.get("link", "") + " " + r.get("snippet", "")
        for pmid in _PMID_RE.findall(text):
            if pmid not in pmids:
                pmids.append(pmid)
    return pmids


def _extract_dois(results: list[dict]) -> list[str]:
    dois: list[str] = []
    for r in results:
        link    = r.get("link", "")
        snippet = r.get("snippet", "")
        text    = link + " " + snippet

        m = _PREPRT_RE.search(link)
        if m:
            doi = _clean_doi("10.1101/" + m.group(1))
            if doi and doi not in dois:
                dois.append(doi)
            continue

        if _is_academic(link):
            for doi in _DOI_RE.findall(text):
                doi = _clean_doi(doi)
                if doi and doi not in dois:
                    dois.append(doi)
    return dois


class _MetaParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.metas: dict[str, str] = {}
        self._done = False

    def handle_starttag(self, tag, attrs):
        if self._done or tag != "meta":
            return
        d    = dict(attrs)
        name = (d.get("name") or d.get("property") or "").lower().strip()
        val  = (d.get("content") or "").strip()
        if name and val:
            self.metas[name] = val

    def handle_endtag(self, tag):
        if tag == "head":
            self._done = True


def _scrape_doi(url: str) -> str:
    if not _is_scrapable(url):
        return ""
    try:
        raw = _page_get(url)
    except Exception:
        return ""
    if not raw:
        return ""

    for doi in _DOI_RE.findall(url):
        doi = _clean_doi(doi)
        if doi:
            return doi

    p = _MetaParser()
    try:
        p.feed(raw.decode("utf-8", errors="replace"))
    except Exception:
        pass
    for key in ("citation_doi", "dc.identifier", "prism.doi",
                "bepress_citation_doi", "dc.identifier.uri"):
        val = p.metas.get(key, "")
        m = _DOI_RE.search(val)
        if m:
            return _clean_doi(m.group())

    m = _DOI_RE.search(raw.decode("utf-8", errors="replace"))
    if m:
        return _clean_doi(m.group())
    return ""


def _resolve_from_results(results: list[dict]) -> tuple[str, str, str]:
    """Return (pmid, doi, source) from Serper results, or ('','','') if nothing found."""
    pmids = _extract_pmids(results)
    if pmids:
        return pmids[0], "", "serper_url"

    dois = _extract_dois(results)
    for doi in dois:
        pmid = _doi_to_pmid_pubmed(doi)
        if pmid:
            return pmid, doi, "serper_doi→pubmed"
        pmid = _doi_to_pmid_europepmc(doi)
        if pmid:
            return pmid, doi, "serper_doi→europepmc"
    if dois:
        return "", dois[0], "serper_doi→crossref_only"

    for r in sorted(results, key=lambda x: x.get("position", 99)):
        link = r.get("link", "")
        if not link:
            continue
        doi = _scrape_doi(link)
        if not doi:
            continue
        pmid = _doi_to_pmid_pubmed(doi)
        if pmid:
            return pmid, doi, "serper_scrape→pubmed"
        pmid = _doi_to_pmid_europepmc(doi)
        if pmid:
            return pmid, doi, "serper_scrape→europepmc"
        return "", doi, "serper_scrape→crossref_only"

    r1 = next((r for r in results if r.get("position") == 1), None)
    if r1 and r1.get("title"):
        pmid = _title_to_pmid(r1["title"])
        if pmid:
            return pmid, "", "serper_title→pubmed"

    return "", "", ""


def _serper_cascade(bp: str, entry: dict, serper_cache: dict) -> tuple[dict, str]:
    """Run the Serper cascade for one BP. Returns (updated_entry, source)."""
    title  = entry.get("title", "")
    cached = serper_cache.setdefault(bp, {})

    queries: list[tuple[str, str]] = [("accession", f'"{bp}"')]
    accession_cached_empty = "accession" in cached and not cached["accession"]
    if title and not accession_cached_empty:
        prefix = " ".join(title.split()[:8])
        queries.append(("title", f'"{bp}" {prefix}'))
        kw = " ".join(title.split()[-3:])
        queries.append(("keywords", f'"{bp}" {kw} RNA-Seq plant pathogen'))

    for qtype, q in queries:
        if qtype not in cached:
            results = _serper_query(q)
            cached[qtype] = results
        else:
            results = cached[qtype]
        if not results:
            continue
        pmid, doi, source = _resolve_from_results(results)
        if pmid or doi:
            return _apply_pmid_doi(entry, pmid, doi, source), source

    return entry, "none"


# ── Applying a resolved pmid/doi to an entry ───────────────────────────────────

def _apply_pmid_doi(base: dict, pmid: str, doi: str, source: str) -> dict:
    """Return a copy of base filled with resolved PMID/DOI + bibliographic metadata.
    Never touches full text — pmcid is recorded (meta_text.py uses it directly)
    but the PMC article itself is not fetched here."""
    entry = dict(base)
    entry["pmid_source"] = source

    if pmid:
        meta = _fetch_pubmed_meta_one(pmid)
        entry["pmid"]        = pmid
        entry["doi"]         = doi or meta.get("doi", "") or base.get("doi", "")
        entry["publication"] = f"[{pmid}] {meta.get('title', '')}"
        entry["abstract"]    = meta.get("abstract", "")
        entry["pmcid"]       = _pmid_to_pmcid(pmid) or ""

    if doi and not pmid:
        cr = _crossref_get(doi)
        entry["doi"]         = doi
        entry["pub_date"]    = cr.get("pub_date", "")
        entry["publication"] = f"[DOI:{doi}] {cr.get('title', '')}"
        entry["abstract"]    = cr.get("abstract", "")
        if not entry.get("title") and cr.get("title"):
            entry["title"] = cr["title"]
        entry["pmcid"] = _doi_to_pmcid(doi) or ""

    return entry


def _doi_only_enrichment(entry: dict) -> dict:
    """For entries with a DOI but no PMID (bioproject_xml_doi_only), try
    DOI->PMID first, else CrossRef for title/abstract/pub_date."""
    doi = entry.get("doi", "")
    if not doi:
        return entry
    out  = dict(entry)
    pmid = _doi_to_pmid_pubmed(doi) or _doi_to_pmid_europepmc(doi)
    if pmid:
        meta = _fetch_pubmed_meta_one(pmid)
        out["pmid"]        = pmid
        out["doi"]         = doi or meta.get("doi", "")
        out["publication"] = f"[{pmid}] {meta.get('title', '')}"
        out["abstract"]    = meta.get("abstract", "")
        out["pmid_source"] = "doi_only→pubmed"
        out["pmcid"]       = _pmid_to_pmcid(pmid) or ""
    else:
        cr = _crossref_get(doi)
        if cr:
            out["pub_date"]    = cr.get("pub_date", "")
            out["publication"] = f"[DOI:{doi}] {cr.get('title', '')}"
            out["abstract"]    = cr.get("abstract", "")
            out["pmid_source"] = "doi_only→crossref"
        out["pmcid"] = _doi_to_pmcid(doi) or ""
    return out


# ── Per-BioProject full pipeline ───────────────────────────────────────────────

def _empty_entry() -> dict:
    return {
        "title": "", "description": "", "submission_date": "",
        "pmid": "", "pmcid": "", "doi": "",
        "pub_date": "", "publication": "", "abstract": "",
        "pmid_source": "none",
        "_v": CACHE_VERSION,
    }


def resolve_bioproject(accession: str, cache: dict, use_serper: bool) -> dict:
    """Full single-pass resolution for one BP: NCBI cascade -> Serper -> CrossRef."""
    entry = cache.get(accession, {})
    if entry.get("_v") == CACHE_VERSION and entry.get("pmid"):
        return entry
    if entry.get("_v") == CACHE_VERSION and not entry.get("pmid") and not entry.get("_retry"):
        return entry  # already tried, no PMID — caller decides whether to retry

    result = _empty_entry()

    uid = None
    try:
        uid = _bp_uid(accession)
    except Exception as e:
        print(f"  WARNING: esearch failed for {accession}: {e}", flush=True)

    xml_pmids: set[str] = set()
    xml_dois:  set[str] = set()
    if uid is not None:
        try:
            meta, xml_pmids, xml_dois = _fetch_bp_xml(uid)
            result.update(meta)
        except Exception as e:
            print(f"  WARNING: BP XML fetch failed for {accession}: {e}", flush=True)

    # Stage 1-3: NCBI PMID cascade
    all_pmids: list[str] = []
    pmid_source = "none"
    strategies = [
        ("bioproject_xml", lambda: list(xml_pmids)),
        ("pmc_search",     lambda: _pmc_search(accession)),
        ("pubmed_search",  lambda: _pubmed_search(accession)),
    ]
    for name, fn in strategies:
        try:
            pmids = fn()
        except Exception as e:
            print(f"  WARNING: {name} failed for {accession}: {e}", flush=True)
            continue
        if pmids:
            all_pmids, pmid_source = pmids, name
            break
    result["pmid_source"] = pmid_source

    if all_pmids:
        try:
            pub_info = _fetch_pubmed_meta(all_pmids)
            if pub_info:
                primary = min(pub_info, key=lambda p: pub_info[p]["date"])
                pri = pub_info[primary]
                result["pmid"]        = primary
                result["doi"]         = pri["doi"]
                result["pub_date"]    = pri["date"]
                result["publication"] = f"[{primary}] {pri['title']}"
                result["abstract"]    = pri["abstract"]
                result["pmcid"]       = _pmid_to_pmcid(primary) or ""
        except Exception as e:
            print(f"  WARNING: PubMed meta failed for {accession}: {e}", flush=True)

    if not result["pmid"] and xml_dois:
        result["doi"]         = next(iter(xml_dois))
        result["pmid_source"] = "bioproject_xml_doi_only"

    # Stage 4: Serper web search (only if NCBI found nothing at all)
    if use_serper and result["pmid_source"] == "none":
        try:
            serper_cache = _SERPER_CACHE_HANDLE
            result, source = _serper_cascade(accession, result, serper_cache)
        except Exception as e:
            print(f"  WARNING: serper cascade failed for {accession}: {e}", flush=True)

    # Stage 5: CrossRef enrichment for DOI-only entries
    if result["pmid_source"] == "bioproject_xml_doi_only":
        try:
            result = _doi_only_enrichment(result)
        except Exception as e:
            print(f"  WARNING: DOI-only enrichment failed for {accession}: {e}", flush=True)

    result["_v"] = CACHE_VERSION
    cache[accession] = result
    return result


# Module-level handle so resolve_bioproject can reach the Serper cache without
# threading it through every call site (Serper results are saved after every
# single BP since queries are expensive — see main()).
_SERPER_CACHE_HANDLE: dict = {}


# ── BioSample batch fetch ──────────────────────────────────────────────────────

def _parse_biosample(elem: ET.Element) -> dict:
    acc = elem.get("accession", "")
    row: dict = {"accession": acc}
    for attr in elem.findall(".//Attribute"):
        harm = (attr.get("harmonized_name") or attr.get("attribute_name") or "")
        harm = harm.lower().strip().replace(" ", "_")
        if harm in _HARMONIZED and attr.text and attr.text.strip():
            row.setdefault(harm, attr.text.strip())
    return row


def _fetch_biosamples(accessions: list[str], cache: dict) -> dict:
    needed = [a for a in accessions
              if a not in cache or set(cache[a].keys()) <= {"accession"}]
    print(f"BioSamples: {len(accessions):,} total, {len(needed):,} to fetch", flush=True)
    results = dict(cache)

    for i in range(0, len(needed), BS_BATCH):
        batch = needed[i: i + BS_BATCH]
        try:
            raw  = _ncbi_post("efetch.fcgi", db="biosample",
                              id=",".join(batch), rettype="xml", retmode="xml")
            root = ET.fromstring(raw)
            for bs in root.findall(".//BioSample"):
                parsed = _parse_biosample(bs)
                acc = parsed.get("accession")
                if acc:
                    results[acc] = parsed
        except Exception as exc:
            print(f"  batch {i // BS_BATCH + 1} error: {exc}", flush=True)

        done = min(i + BS_BATCH, len(needed))
        if done % 3000 == 0 or done == len(needed):
            print(f"  {done:,} / {len(needed):,}", flush=True)
            save_json(results, BS_CACHE)

    for acc in needed:
        results.setdefault(acc, {"accession": acc})
    save_json(results, BS_CACHE)
    return results


def _fetch_ebi_biosamples(accessions: list[str], cache: dict) -> dict:
    """Fetch ENA/DDBJ BioSample attributes via EBI BioSamples API for SAME*/SAMD* accessions."""
    ebi_accessions = [a for a in accessions if a[:4] in _ENA_PREFIXES]
    needed = [a for a in ebi_accessions if set(cache.get(a, {}).keys()) <= {"accession"}]
    if not needed:
        return cache
    print(f"EBI BioSamples: {len(ebi_accessions):,} ENA/DDBJ total, "
          f"{len(needed):,} to fetch", flush=True)
    results = dict(cache)
    for i, acc in enumerate(needed, 1):
        results[acc] = _fetch_ena_biosample(acc)
        if i % 500 == 0 or i == len(needed):
            print(f"  {i:,} / {len(needed):,}", flush=True)
            save_json(results, BS_CACHE)
    save_json(results, BS_CACHE)
    return results


# ── Migration from legacy ncbi_metadata.py / web_metadata.py caches ───────────

def migrate_legacy() -> None:
    """One-time: merge the two legacy bp_cache.json files into the new cache,
    dropping methods_text (meta_text.py owns full text now). Also copies
    bs_cache.json and serper_cache.json forward so nothing is re-fetched."""
    if not LEGACY_NCBI_CACHE.exists():
        print(f"  no legacy NCBI cache found at {LEGACY_NCBI_CACHE} — nothing to migrate", flush=True)
        return

    ncbi = load_json(LEGACY_NCBI_CACHE)
    web  = load_json(LEGACY_WEB_CACHE) if LEGACY_WEB_CACHE.exists() else {}
    print(f"  legacy NCBI cache:  {len(ncbi):,} BPs", flush=True)
    print(f"  legacy web cache:   {len(web):,} BPs (overlay)", flush=True)

    merged: dict = {}
    for bp, entry in ncbi.items():
        e = dict(entry)
        e.pop("methods_text", None)
        merged[bp] = e
    for bp, entry in web.items():
        e = dict(entry)
        e.pop("methods_text", None)
        merged[bp] = e  # web-resolved entries win — same precedence as old pipeline

    save_json(merged, BP_CACHE)
    print(f"  written: {BP_CACHE}  ({len(merged):,} BPs, methods_text stripped)", flush=True)

    if LEGACY_NCBI_BS.exists() and not BS_CACHE.exists():
        bs = load_json(LEGACY_NCBI_BS)
        save_json(bs, BS_CACHE)
        print(f"  copied BioSample cache: {len(bs):,} entries", flush=True)

    serper_src = LEGACY_WEB_SERPER if LEGACY_WEB_SERPER.exists() else None
    if serper_src and not SERPER_CACHE.exists():
        sc = load_json(serper_src)
        save_json(sc, SERPER_CACHE)
        print(f"  copied Serper cache: {len(sc):,} BPs", flush=True)


def bootstrap_from_legacy(cache: dict) -> int:
    """Import NCBI-sourced entries from the old fetch_lit lit_cache.json."""
    if not LEGACY_LIT_CACHE.exists():
        print(f"  legacy cache not found: {LEGACY_LIT_CACHE}", flush=True)
        return 0
    old = load_json(LEGACY_LIT_CACHE)
    n = 0
    for bp, entry in old.items():
        if bp in cache and cache[bp].get("_v") == CACHE_VERSION:
            continue
        src = entry.get("pmid_source", "none")
        if any(src.lower().startswith(s.lower()) for s in ("Serper", "manual")):
            continue
        new = _empty_entry()
        new["title"]           = entry.get("title", "")
        new["description"]     = entry.get("description", "")
        new["submission_date"] = entry.get("bp_submission_date", "")
        new["pmid"]            = entry.get("primary_pmid", "")
        new["pmcid"]           = entry.get("pmcid", "")
        new["doi"]             = entry.get("primary_doi", "")
        new["pub_date"]        = entry.get("primary_pub_date", "")
        new["publication"]     = entry.get("primary_publication", "")
        new["abstract"]        = entry.get("abstract", "")
        src_map = {
            "BioProject XML": "bioproject_xml",
            "BioProject XML DOI": "bioproject_xml_doi_only",
            "BioProject XML DOI (CrossRef)": "bioproject_xml_doi_only",
            "PMC full-text": "pmc_search",
            "PubMed": "pubmed_search",
            "cached": "bioproject_xml",
            "none": "none",
        }
        new["pmid_source"] = src_map.get(src, src)
        cache[bp] = new
        n += 1
    return n


def bootstrap_serper(serper_cache: dict) -> tuple[int, int]:
    """Seed serper_cache from the legacy serper_results.tsv."""
    if not LEGACY_SERPER_TSV.exists():
        print(f"  not found: {LEGACY_SERPER_TSV}", flush=True)
        return 0, 0
    by_bp: dict[str, list[dict]] = {}
    with open(LEGACY_SERPER_TSV, encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            by_bp.setdefault(row["bioproject"], []).append(row)
    n_seeded = n_skipped = 0
    for bp, rows in by_bp.items():
        if "accession" in serper_cache.get(bp, {}):
            n_skipped += 1
            continue
        results = [
            {"title": r.get("title", ""), "link": r.get("link", ""),
             "snippet": r.get("snippet", ""), "position": int(r.get("rank") or 99)}
            for r in rows if r.get("rank", "0") not in ("0", "")
        ]
        serper_cache.setdefault(bp, {})["accession"] = results
        n_seeded += 1
    return n_seeded, n_skipped


# ── Report ─────────────────────────────────────────────────────────────────────

def print_report(cache: dict, bioprojects: list[str]) -> None:
    total    = len(bioprojects)
    has_pmid = sum(1 for bp in bioprojects if cache.get(bp, {}).get("pmid"))
    has_doi  = sum(1 for bp in bioprojects if cache.get(bp, {}).get("doi"))
    has_pmc  = sum(1 for bp in bioprojects if cache.get(bp, {}).get("pmcid"))
    has_abs  = sum(1 for bp in bioprojects if cache.get(bp, {}).get("abstract"))

    sources: collections.Counter = collections.Counter(
        cache.get(bp, {}).get("pmid_source", "not_fetched") for bp in bioprojects
    )

    print(f"\n{'─'*55}")
    print(f"meta_search — resolution report")
    print(f"{'─'*55}")
    print(f"BioProjects in runs.tsv:  {total:>6,}")
    print(f"  has PMID:               {has_pmid:>6,}  ({100*has_pmid/max(total,1):.1f}%)")
    print(f"  has DOI:                {has_doi:>6,}  ({100*has_doi/max(total,1):.1f}%)")
    print(f"  has PMCID:              {has_pmc:>6,}  ({100*has_pmc/max(total,1):.1f}%)")
    print(f"  has abstract:           {has_abs:>6,}  ({100*has_abs/max(total,1):.1f}%)")
    print(f"\npmid_source breakdown:")
    for src, n in sorted(sources.items(), key=lambda x: -x[1]):
        print(f"  {src:<35} {n:>5}  ({100*n/max(total,1):.1f}%)")
    print(f"{'─'*55}\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    global _SERPER_CACHE_HANDLE

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--migrate-legacy", action="store_true",
                    help="merge legacy ncbi_metadata/web_metadata bp_cache.json files, "
                         "dropping methods_text, and exit")
    ap.add_argument("--retry",     action="store_true",
                    help="re-attempt BPs with no PMID/DOI found on a previous run")
    ap.add_argument("--bootstrap", action="store_true",
                    help="import NCBI-sourced entries from legacy lit_cache.json")
    ap.add_argument("--bootstrap-serper", action="store_true",
                    help="import legacy serper_results.tsv into serper_cache.json")
    ap.add_argument("--fetch-ena", action="store_true",
                    help="fetch ENA/DDBJ BioSample attrs via EBI API only; skip everything else")
    ap.add_argument("--report",  action="store_true",
                    help="print resolution breakdown and exit")
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve without writing bp_cache (serper_cache still written)")
    ap.add_argument("--limit", type=int, default=0,
                    help="process at most N BPs (0 = all)")
    args = ap.parse_args()

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"

    if not args.report:
        log_dir = make_log_dir(logs_base)
        sys.stdout = _Tee(log_dir / "fetch.log")
        link_latest(logs_base, log_dir / "fetch.log")

    if args.migrate_legacy:
        migrate_legacy()
        return

    print(f"NCBI_API_KEY: {'yes' if NCBI_KEY else 'no'}", flush=True)
    use_serper = bool(SERPER_KEY)
    if not use_serper and not args.report and not args.bootstrap_serper:
        print("WARNING: SERPER_API_KEY not set — stage 4 (web search) will be skipped; "
              "NCBI-only resolution still runs.", flush=True)

    if not RUNS_TSV.exists():
        sys.exit(f"ERROR: {RUNS_TSV} not found — run stat/stat_filter.py first")

    bioprojects: list[str] = []
    biosamples:  set[str]  = set()
    bp_seen: set[str] = set()
    with open(RUNS_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            bp = row.get("BioProject", "").strip()
            if bp and bp not in bp_seen:
                bp_seen.add(bp)
                bioprojects.append(bp)
            bs = row.get("BioSample", "").strip()
            if bs:
                biosamples.add(bs)
    print(f"BioProjects: {len(bioprojects):,}  BioSamples: {len(biosamples):,}", flush=True)

    bp_cache = load_json(BP_CACHE) if BP_CACHE.exists() else {}

    if args.bootstrap:
        n_boot = bootstrap_from_legacy(bp_cache)
        print(f"Bootstrapped {n_boot:,} entries from legacy lit_cache", flush=True)
        save_json(bp_cache, BP_CACHE)

    serper_cache = load_json(SERPER_CACHE) if SERPER_CACHE.exists() else {}
    _SERPER_CACHE_HANDLE = serper_cache

    if args.bootstrap_serper:
        n_seeded, n_skipped = bootstrap_serper(serper_cache)
        save_json(serper_cache, SERPER_CACHE)
        print(f"Bootstrap complete: seeded {n_seeded:,}  skipped {n_skipped:,} (already cached)")
        return

    if args.report:
        print_report(bp_cache, bioprojects)
        return

    if args.fetch_ena:
        bs_cache = load_json(BS_CACHE) if BS_CACHE.exists() else {}
        bs_out   = _fetch_ebi_biosamples(sorted(biosamples), bs_cache)
        bs_clean = {k: v for k, v in bs_out.items() if k in biosamples}
        bs_path  = OUT_DIR / "data" / "biosamples.json"
        bs_path.write_text(json.dumps(bs_clean, indent=2))
        print(f"Written: {bs_path}  ({len(bs_clean):,} entries)", flush=True)
        n = len(bs_clean)
        print("BioSample attribute coverage:")
        for field in sorted(_HARMONIZED | {"bs_description"}):
            count = sum(1 for v in bs_clean.values() if v.get(field))
            print(f"  {field:<20} {count:>6,}  ({100*count/max(n,1):.1f}%)", flush=True)
        return

    # ── BioProject resolution ─────────────────────────────────────────────────
    to_fetch: list[str] = []
    to_retry: list[str] = []
    n_skip = 0
    for bp in bioprojects:
        entry = bp_cache.get(bp, {})
        if entry.get("_v") == CACHE_VERSION and entry.get("pmid"):
            n_skip += 1
        elif entry.get("_v") == CACHE_VERSION and not entry.get("pmid"):
            if args.retry:
                to_retry.append(bp)
            else:
                n_skip += 1
        else:
            to_fetch.append(bp)

    if args.limit:
        to_fetch = to_fetch[: args.limit]
        to_retry = to_retry[: max(0, args.limit - len(to_fetch))]

    total_work = len(to_fetch) + len(to_retry)
    print(f"  skip (cached): {n_skip:,}", flush=True)
    print(f"  fetch (new):   {len(to_fetch):,}", flush=True)
    if args.retry:
        print(f"  retry:         {len(to_retry):,}", flush=True)
    print(f"  → work items:  {total_work:,}", flush=True)

    work_list = [("fetch", bp) for bp in to_fetch] + [("retry", bp) for bp in to_retry]

    for i, (kind, bp) in enumerate(work_list, 1):
        if kind == "retry":
            bp_cache.get(bp, {})["_retry"] = True
        print(f"  [{i:>4}/{total_work}] {kind:<6} {bp}", end=" · ", flush=True)
        try:
            meta = resolve_bioproject(bp, bp_cache, use_serper)
        except Exception as e:
            print(f"ERROR: {e}", flush=True)
            continue
        print(
            f"source: {meta.get('pmid_source','?'):<24} "
            f"pmcid: {meta.get('pmcid') or '-':<12} "
            f"doi: {(meta.get('doi') or '-')[:30]}",
            flush=True,
        )
        if i % 20 == 0:
            save_json(bp_cache, BP_CACHE)
            save_json(serper_cache, SERPER_CACHE)

    if not args.dry_run:
        save_json(bp_cache, BP_CACHE)
    save_json(serper_cache, SERPER_CACHE)

    # ── Emit bioprojects.json ─────────────────────────────────────────────────
    if not args.dry_run:
        bp_out = {bp: bp_cache.get(bp, _empty_entry()) for bp in bioprojects}
        bp_path = OUT_DIR / "data" / "bioprojects.json"
        bp_path.write_text(json.dumps(bp_out, indent=2))
        print(f"Written: {bp_path}  ({len(bp_out):,} entries)", flush=True)

    # ── BioSamples ────────────────────────────────────────────────────────────
    bs_cache = load_json(BS_CACHE) if BS_CACHE.exists() else {}
    bs_out   = _fetch_biosamples(sorted(biosamples), bs_cache)
    bs_out   = _fetch_ebi_biosamples(sorted(biosamples), bs_out)
    bs_clean = {k: v for k, v in bs_out.items() if k in biosamples}
    bs_path  = OUT_DIR / "data" / "biosamples.json"
    bs_path.write_text(json.dumps(bs_clean, indent=2))
    print(f"Written: {bs_path}  ({len(bs_clean):,} entries)", flush=True)

    print_report(bp_cache, bioprojects)

    n = len(bs_clean)
    print("BioSample attribute coverage:")
    for field in sorted(_HARMONIZED | {"bs_description"}):
        count = sum(1 for v in bs_clean.values() if v.get(field))
        print(f"  {field:<20} {count:>6,}  ({100*count/max(n,1):.1f}%)", flush=True)


if __name__ == "__main__":
    main()
