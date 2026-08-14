#!/usr/bin/env python3
"""
serper_scrape.py — scrape journal pages from saved Serper links for DOIs.

For each unresolved BioProject (no PMID, no DOI), fetches the top-ranked
academic links from serper_results.tsv and extracts DOIs from HTML <meta>
tags (Highwire Press citation_doi / DC.Identifier standard). Falls back to
PubMed title search using the Serper result title.

DOI is stored as the primary identifier; PMID resolved when possible via
PubMed → EuropePMC. Metadata enriched via CrossRef.

No Serper API calls — reads the already-saved TSV.

Usage:
    python metadata/serper_scrape.py [--dry-run] [--limit N]
"""

import argparse
import csv
import html.parser
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

SERPER_TSV  = Path("metadata/output/serper/serper_results.tsv")
CACHE_PATH  = Path("metadata/output/fetch_lit/data/lit_cache.json")
LIT_OUT     = Path("metadata/output/fetch_lit/data/literature.json")

API_KEY  = os.environ.get("NCBI_API_KEY", "")
NCBI_RATE = 9.0 if API_KEY else 2.5
PAGE_RATE = 1.0   # req/s for journal page fetches
ENTREZ    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC      = "https://www.ebi.ac.uk/europepmc/webservices/rest"
CROSSREF  = "https://api.crossref.org/works"

HEADERS_NCBI = {"User-Agent": "crypt/serper_scrape (leon.lenzo@curtin.edu.au)"}
HEADERS_EXT  = {"User-Agent": "crypt/serper_scrape (leon.lenzo@curtin.edu.au)"}
HEADERS_PAGE = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0",
    "Accept": "text/html,application/xhtml+xml",
    "Accept-Language": "en-US,en;q=0.9",
}
CROSSREF_HEADERS = {"User-Agent": "crypt/serper_scrape (mailto:leon.lenzo@curtin.edu.au)"}

CACHE_VERSION = 5
_LIT_FIELDS = [
    "primary_pmid", "primary_doi", "primary_pub_date", "primary_publication",
    "abstract", "methods_text", "pmcid", "n_papers_found", "pmid_source",
]

# Domains where Highwire/DC meta tags are expected and page fetching is polite
_SCRAPABLE = {
    "mdpi.com", "frontiersin.org", "journals.plos.org", "elifesciences.org",
    "peerj.com", "f1000research.com", "preprints.org", "researchsquare.com",
    "biorxiv.org", "medrxiv.org",
    "pubs.acs.org",                         # DOI in URL + meta tags
    "apsjournals.apsnet.org",               # APS plant pathology journals
    "phytobiomes.org", "apsnet.org", "mpmi.org",
    "plantcell.org", "plantphysiol.org", "journals.asm.org",
    "academic.oup.com",                     # OUP (sometimes works)
    "europepmc.org",
    "doi.org",                              # follow redirect to landing page
}

# Domains to never attempt (bot-protected or no useful meta)
_SKIP = {
    "researchgate.net", "omicsdi.org", "seqout.org",
    "agdatacommons.nal.usda.gov", "thbif.onep.go.th", "gold.jgi.doe.gov",
    "ddbj.nig.ac.jp", "cucurbitgenomics.org", "triticeaeexpdb.cn",
    "heatomics.sdau.edu.cn", "barleyexp.com", "sciencedirect.com",
}

_DOI_RE  = re.compile(r'10\.\d{4,}/[^\s"\'<>]+')
_PMID_RE = re.compile(r'(?:pubmed\.ncbi\.nlm\.nih\.gov/|pmid[=:/\s]+)(\d{6,9})', re.I)

_last_ncbi: float = 0.0
_last_page: float = 0.0
_last_ext:  float = 0.0
_last_cr:   float = 0.0


def _ncbi_wait():
    global _last_ncbi
    gap = 1.0 / NCBI_RATE
    wait = _last_ncbi + gap - time.monotonic()
    if wait > 0: time.sleep(wait)
    _last_ncbi = time.monotonic()

def _page_wait():
    global _last_page
    gap = 1.0 / PAGE_RATE
    wait = _last_page + gap - time.monotonic()
    if wait > 0: time.sleep(wait)
    _last_page = time.monotonic()

def _ext_wait():
    global _last_ext
    gap = 0.5
    wait = _last_ext + gap - time.monotonic()
    if wait > 0: time.sleep(wait)
    _last_ext = time.monotonic()

def _cr_wait():
    global _last_cr
    gap = 1.0
    wait = _last_cr + gap - time.monotonic()
    if wait > 0: time.sleep(wait)
    _last_cr = time.monotonic()


def _clean_doi(doi: str) -> str:
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


def _ncbi_get(path: str, **params) -> bytes:
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{ENTREZ}/{path}?{urllib.parse.urlencode(params)}"
    _ncbi_wait()
    return urllib.request.urlopen(
        urllib.request.Request(url, headers=HEADERS_NCBI), timeout=30
    ).read()


def _doi_to_pmid(doi: str) -> str:
    try:
        raw   = _ncbi_get("esearch.fcgi", db="pubmed", term=f"{doi}[doi]", retmax=1)
        root  = ET.fromstring(raw)
        count = int(root.findtext(".//Count") or "0")
        if count != 1:
            return ""
        ids = [el.text for el in root.findall(".//Id") if el.text]
        return ids[0] if ids else ""
    except Exception:
        return ""


def _europepmc_doi_to_pmid(doi: str) -> str:
    query = urllib.parse.quote(f'DOI:"{doi}"')
    url   = f"{EPMC}/search?query={query}&format=json&pageSize=3&fields=pmid"
    _ext_wait()
    try:
        raw  = urllib.request.urlopen(
            urllib.request.Request(url, headers=HEADERS_EXT), timeout=20
        ).read()
        for r in json.loads(raw).get("resultList", {}).get("result", []):
            pmid = (r.get("pmid") or "").strip()
            if pmid:
                return pmid
    except Exception:
        pass
    return ""


def _crossref_get(doi: str) -> dict:
    url = f"{CROSSREF}/{urllib.parse.quote(doi, safe='')}"
    _cr_wait()
    try:
        raw  = urllib.request.urlopen(
            urllib.request.Request(url, headers=CROSSREF_HEADERS), timeout=20
        ).read()
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


def _fetch_pubmed_meta(pmid: str) -> dict:
    try:
        raw  = _ncbi_get("efetch.fcgi", db="pubmed", id=pmid, rettype="xml", retmode="xml")
        root = ET.fromstring(raw)
        for article in root.findall(".//PubmedArticle"):
            pid_el  = article.find(".//MedlineCitation/PMID")
            ttl_el  = article.find(".//MedlineCitation/Article/ArticleTitle")
            abs_els = article.findall(".//MedlineCitation/Article/Abstract/AbstractText")
            if pid_el is None:
                continue
            title    = "".join(ttl_el.itertext()).strip() if ttl_el is not None else ""
            abstract = " ".join("".join(el.itertext()).strip() for el in abs_els).strip()
            doi      = ""
            for aid in article.findall(".//ArticleIdList/ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = (aid.text or "").strip()
                    break
            return {"title": title, "abstract": abstract, "doi": doi}
    except Exception:
        pass
    return {}


def _title_to_pmid(title: str) -> str:
    if len(title.split()) < 6:
        return ""
    query = " ".join(title.split()[:12]) + "[ti]"
    try:
        raw   = _ncbi_get("esearch.fcgi", db="pubmed", term=query, retmax=2)
        root  = ET.fromstring(raw)
        count = int(root.findtext(".//Count") or "0")
        ids   = [el.text for el in root.findall(".//Id") if el.text]
        return ids[0] if count == 1 else ""
    except Exception:
        return ""


# ── HTML meta-tag parser ──────────────────────────────────────────────────────

class _MetaParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.metas: dict[str, str] = {}
        self._past_head = False

    def handle_starttag(self, tag, attrs):
        if self._past_head:
            return
        if tag == "body":
            self._past_head = True
            return
        if tag == "meta":
            d = dict(attrs)
            name    = (d.get("name") or d.get("property") or "").lower().strip()
            content = (d.get("content") or "").strip()
            if name and content:
                self.metas[name] = content

    def handle_endtag(self, tag):
        if tag == "head":
            self._past_head = True


def _parse_metas(raw: bytes) -> dict[str, str]:
    p = _MetaParser()
    try:
        p.feed(raw.decode("utf-8", errors="replace"))
    except Exception:
        pass
    return p.metas


def _domain(url: str) -> str:
    return urllib.parse.urlparse(url).netloc.lstrip("www.")


def _is_scrapable(url: str) -> bool:
    d = _domain(url)
    if any(d == s or d.endswith("." + s) for s in _SKIP):
        return False
    return any(d == s or d.endswith("." + s) for s in _SCRAPABLE)


def _fetch_page(url: str) -> bytes:
    _page_wait()
    req = urllib.request.Request(url, headers=HEADERS_PAGE)
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return r.read(32768)
    except Exception:
        return b""


def _doi_from_link(link: str, snippet: str) -> str:
    """Extract a DOI from a single Serper result. Returns cleaned DOI or ''."""
    text = link + " " + snippet

    # DOI directly in URL (Wiley, Springer, APS, ACS embed DOI in path)
    for doi in _DOI_RE.findall(link):
        doi = _clean_doi(doi)
        if doi:
            return doi

    if not _is_scrapable(link):
        return ""

    raw = _fetch_page(link)
    if not raw:
        return ""

    metas = _parse_metas(raw)

    # Citation meta tags — check multiple name conventions
    for key in ("citation_doi", "dc.identifier", "prism.doi",
                "bepress_citation_doi", "dc.identifier.uri"):
        val = metas.get(key, "")
        m = _DOI_RE.search(val)
        if m:
            return _clean_doi(m.group())

    # DOI anywhere in raw HTML body (last resort)
    m = _DOI_RE.search(raw[:32768].decode("utf-8", errors="replace"))
    if m:
        return _clean_doi(m.group())

    return ""


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    if not SERPER_TSV.exists():
        sys.exit(f"ERROR: {SERPER_TSV} not found")

    by_bp: dict[str, list[dict]] = {}
    with open(SERPER_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            by_bp.setdefault(row["bioproject"], []).append(row)

    cache = load_json(CACHE_PATH) if CACHE_PATH.exists() else {}

    to_scrape = [
        bp for bp in by_bp
        if not cache.get(bp, {}).get("primary_pmid")
        and not cache.get(bp, {}).get("primary_doi")
        and any(r.get("rank") not in ("0", "") for r in by_bp[bp])
    ]
    if args.limit:
        to_scrape = to_scrape[:args.limit]

    print(f"Unresolved BPs to scrape: {len(to_scrape)}"
          + (f"  (capped at {args.limit})" if args.limit else ""))

    n_pmid = n_doi_only = 0

    for i, bp in enumerate(to_scrape, 1):
        rows = sorted(by_bp[bp], key=lambda r: int(r.get("rank") or 0))
        print(f"[{i:>4}/{len(to_scrape)}] {bp}", end=" · ", flush=True)

        doi = pmid = source = ""

        # Strategy 1: scrape links in rank order for a DOI
        for r in rows:
            if r.get("rank") in ("0", ""):
                continue
            link    = r.get("link", "")
            snippet = r.get("snippet", "")
            if not link:
                continue
            doi = _doi_from_link(link, snippet)
            if doi:
                source = "scrape"
                break

        # Strategy 2: PubMed title search on rank-1 title
        if not doi:
            r1 = next((r for r in rows if r.get("rank") == "1"), None)
            if r1 and r1.get("title"):
                pmid = _title_to_pmid(r1["title"])
                if pmid:
                    source = "title_search"

        if not doi and not pmid:
            print("nothing found")
            continue

        # Resolve DOI → PMID if possible
        if doi and not pmid:
            pmid = _doi_to_pmid(doi)
            if not pmid:
                pmid = _europepmc_doi_to_pmid(doi)

        print(f"doi={doi or '-'}  pmid={pmid or '-'}  via {source}", end=" · ", flush=True)

        if args.dry_run:
            print("(dry-run)")
            if pmid:
                n_pmid += 1
            else:
                n_doi_only += 1
            continue

        entry = cache.get(bp, {})

        if pmid:
            meta = _fetch_pubmed_meta(pmid)
            entry.update({
                "primary_pmid":        pmid,
                "primary_doi":         doi or meta.get("doi", ""),
                "primary_pub_date":    "",
                "primary_publication": f"[{pmid}] {meta.get('title', '')}",
                "abstract":            meta.get("abstract", ""),
                "n_papers_found":      1,
                "pmid_source":         f"Serper+{source}",
                "pmcid":               entry.get("pmcid", ""),
                "methods_text":        entry.get("methods_text", ""),
                "_v":                  CACHE_VERSION,
            })
            print(f"saved → [{pmid}] {meta.get('title','')[:50]}")
            n_pmid += 1
        else:
            cr = _crossref_get(doi)
            entry.update({
                "primary_pmid":        "",
                "primary_doi":         doi,
                "primary_pub_date":    cr.get("pub_date", ""),
                "primary_publication": f"[DOI:{doi}] {cr.get('title', '')}",
                "abstract":            cr.get("abstract", ""),
                "n_papers_found":      1,
                "pmid_source":         f"Serper+{source}",
                "pmcid":               entry.get("pmcid", ""),
                "methods_text":        entry.get("methods_text", ""),
                "_v":                  CACHE_VERSION,
            })
            print(f"saved DOI → {doi[:50]}  {cr.get('title','')[:35]}")
            n_doi_only += 1

        cache[bp] = entry
        save_json(cache, CACHE_PATH)

    total = n_pmid + n_doi_only
    print(f"\nResolved: {total}/{len(to_scrape)} BPs  "
          f"(PMID: {n_pmid}, DOI-only: {n_doi_only})")

    if not args.dry_run and total > 0:
        lit_cache = load_json(CACHE_PATH)
        lit_out = {bp: {k: entry.get(k, "") for k in _LIT_FIELDS}
                   for bp, entry in lit_cache.items()}
        LIT_OUT.write_text(json.dumps(lit_out, indent=2))
        print(f"Written: {LIT_OUT}")


if __name__ == "__main__":
    main()
