#!/usr/bin/env python3
"""
serper_resolve.py — resolve DOIs and PMIDs from saved serper_results.tsv.

Reads metadata/output/serper/serper_results.tsv, extracts DOIs from links
and snippets, resolves to PMIDs via PubMed → EuropePMC, and enriches
DOI-only results via CrossRef. Updates lit_cache.json + literature.json.

DOI is stored as the primary identifier; PMID is secondary when available.

Usage:
    python metadata/serper_resolve.py [--dry-run]
"""

import argparse
import csv
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
RATE     = 9.0 if API_KEY else 2.5
ENTREZ   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC     = "https://www.ebi.ac.uk/europepmc/webservices/rest"
CROSSREF = "https://api.crossref.org/works"
HEADERS       = {"User-Agent": "crypt/serper_resolve (leon.lenzo@curtin.edu.au)"}
EXT_HEADERS   = {"User-Agent": "crypt/serper_resolve (leon.lenzo@curtin.edu.au)"}
CROSSREF_HEADERS = {"User-Agent": "crypt/serper_resolve (mailto:leon.lenzo@curtin.edu.au)"}

_DOI_RE     = re.compile(r'10\.\d{4,}/[^\s"\'<>]+')
_PMID_RE    = re.compile(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d{6,9})')
_BIORXIV_RE = re.compile(r'biorxiv\.org/content/[^\s]*?(\d{4}\.\d{2}\.\d{2}\.\d+)')
_MEDRXIV_RE = re.compile(r'medrxiv\.org/content/[^\s]*?(\d{4}\.\d{2}\.\d{2}\.\d+)')

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

CACHE_VERSION = 5
_LIT_FIELDS = [
    "primary_pmid", "primary_doi", "primary_pub_date", "primary_publication",
    "abstract", "methods_text", "pmcid", "n_papers_found", "pmid_source",
]

_last_req: float = 0.0
_last_ext: float = 0.0
_last_cr:  float = 0.0


def _rate_wait() -> None:
    global _last_req
    gap  = 1.0 / RATE
    wait = _last_req + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_req = time.monotonic()


def _ext_wait() -> None:
    global _last_ext
    gap  = 0.5
    wait = _last_ext + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_ext = time.monotonic()


def _cr_wait() -> None:
    global _last_cr
    gap  = 1.0
    wait = _last_cr + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_cr = time.monotonic()


def _get(path: str, **params) -> bytes:
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{ENTREZ}/{path}?{urllib.parse.urlencode(params)}"
    _rate_wait()
    return urllib.request.urlopen(
        urllib.request.Request(url, headers=HEADERS), timeout=30
    ).read()


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
    """EuropePMC DOI search — broader coverage than PubMed esearch."""
    query = urllib.parse.quote(f'DOI:"{doi}"')
    url   = f"{EPMC}/search?query={query}&format=json&pageSize=3&fields=pmid"
    _ext_wait()
    try:
        raw  = urllib.request.urlopen(
            urllib.request.Request(url, headers=EXT_HEADERS), timeout=20
        ).read()
        data = json.loads(raw)
    except Exception:
        return ""
    for r in data.get("resultList", {}).get("result", []):
        pmid = (r.get("pmid") or "").strip()
        if pmid:
            return pmid
    return ""


def _crossref_get(doi: str) -> dict:
    """Fetch title, abstract, pub_date from CrossRef (universal DOI coverage)."""
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
        raw  = _get("efetch.fcgi", db="pubmed", id=pmid, rettype="xml", retmode="xml")
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


def _clean_doi(doi: str) -> str:
    """Strip URL artifacts (query strings, attachment paths, PDF suffixes) from a DOI."""
    doi = doi.split("?")[0]                        # remove ?query=... parameters
    doi = doi.rstrip(".,;)")
    if "/attachment/" in doi:                       # Elsevier/Cell supplementary paths
        doi = doi.split("/attachment/")[0]
    if "~" in doi:                                  # URL slug separator (some CMS systems)
        doi = doi.split("~")[0]
    # OUP PDF paths: 10.1093/journal/id/NNN/filename.pdf — 4th segment is numeric
    parts = doi.split("/")
    if len(parts) >= 4 and parts[3].isdigit():
        doi = "/".join(parts[:3])
    for sfx in ("/full.pdf", ".full.pdf", "/full", ".pdf"):
        if doi.endswith(sfx):
            doi = doi[: -len(sfx)]
    return doi.rstrip(".,;)")


def extract_dois_from_rows(rows: list[dict]) -> list[str]:
    """Extract candidate DOIs from Serper result rows (link + snippet)."""
    dois: list[str] = []
    for r in rows:
        if r.get("rank") == "0":
            continue
        link    = r.get("link", "")
        snippet = r.get("snippet", "")
        text    = link + " " + snippet

        # bioRxiv / medRxiv: construct 10.1101/{id} DOI from URL
        for pat in [_BIORXIV_RE, _MEDRXIV_RE]:
            m = pat.search(link)
            if m:
                doi = _clean_doi("10.1101/" + m.group(1))
                if doi not in dois:
                    dois.append(doi)

        # Any DOI from academic domains
        domain = urllib.parse.urlparse(link).netloc.lstrip("www.")
        is_academic = any(domain == d or domain.endswith("." + d)
                          for d in _ACADEMIC_DOMAINS)
        if is_academic:
            for doi in _DOI_RE.findall(text):
                doi = _clean_doi(doi)
                if doi not in dois:
                    dois.append(doi)

    return dois


def extract_pmids_from_rows(rows: list[dict]) -> list[str]:
    """Fast-path: extract PMIDs from direct PubMed URLs in Serper results."""
    pmids: list[str] = []
    for r in rows:
        if r.get("rank") == "0":
            continue
        text = r.get("link", "") + " " + r.get("snippet", "")
        for pmid in _PMID_RE.findall(text):
            if pmid not in pmids:
                pmids.append(pmid)
    return pmids


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be updated without writing cache")
    args = ap.parse_args()

    if not SERPER_TSV.exists():
        sys.exit(f"ERROR: {SERPER_TSV} not found — run serper_dump.py first")

    # Group rows by BP
    by_bp: dict[str, list[dict]] = {}
    with open(SERPER_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            bp = row["bioproject"]
            by_bp.setdefault(bp, []).append(row)

    cache = load_json(CACHE_PATH) if CACHE_PATH.exists() else {}

    # Process BPs with no PMID AND no DOI (fully unresolved)
    to_resolve = [
        bp for bp in by_bp
        if not cache.get(bp, {}).get("primary_pmid")
        and not cache.get(bp, {}).get("primary_doi")
        and any(r.get("rank") != "0" for r in by_bp[bp])
    ]
    print(f"BPs with Serper results, no PMID, no DOI: {len(to_resolve)}")

    n_pmid = n_doi_only = 0
    for i, bp in enumerate(to_resolve, 1):
        rows = by_bp[bp]
        print(f"[{i:>4}/{len(to_resolve)}] {bp}", end=" · ", flush=True)

        # Fast-path: direct PubMed URL in snippet/link
        pmids = extract_pmids_from_rows(rows)
        if pmids:
            pmid = pmids[0]
            doi_found = ""
            source = "Serper+url"
        else:
            # Collect DOIs from URLs/snippets
            dois = extract_dois_from_rows(rows)
            if not dois:
                print("no DOI or PMID extracted")
                continue

            doi_found = ""
            pmid      = ""
            source    = ""

            for doi in dois:
                # Try PubMed first
                pmid = _doi_to_pmid(doi)
                if pmid:
                    doi_found = doi
                    source    = "Serper+doi→PubMed"
                    break
                # Try EuropePMC (broader coverage)
                pmid = _europepmc_doi_to_pmid(doi)
                if pmid:
                    doi_found = doi
                    source    = "Serper+doi→EuropePMC"
                    break
                # No PMID — store DOI directly via CrossRef
                doi_found = doi
                source    = "Serper+doi→CrossRef"
                break  # use first DOI found

        if not pmid and not doi_found:
            print("no DOI or PMID extracted")
            continue

        print(f"doi={doi_found or '-'}  pmid={pmid or '-'}  via {source}", end=" · ", flush=True)

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
                "primary_doi":         doi_found or meta.get("doi", ""),
                "primary_pub_date":    "",
                "primary_publication": f"[{pmid}] {meta.get('title', '')}",
                "abstract":            meta.get("abstract", ""),
                "n_papers_found":      1,
                "pmid_source":         source,
                "pmcid":               entry.get("pmcid", ""),
                "methods_text":        entry.get("methods_text", ""),
                "_v":                  CACHE_VERSION,
            })
            print(f"saved → [{pmid}] {meta.get('title','')[:50]}")
            n_pmid += 1
        else:
            # DOI-only: enrich from CrossRef
            cr = _crossref_get(doi_found)
            entry.update({
                "primary_pmid":        "",
                "primary_doi":         doi_found,
                "primary_pub_date":    cr.get("pub_date", ""),
                "primary_publication": f"[DOI:{doi_found}] {cr.get('title', '')}",
                "abstract":            cr.get("abstract", ""),
                "n_papers_found":      1,
                "pmid_source":         source,
                "pmcid":               entry.get("pmcid", ""),
                "methods_text":        entry.get("methods_text", ""),
                "_v":                  CACHE_VERSION,
            })
            print(f"saved DOI → {doi_found}  {cr.get('title','')[:45]}")
            n_doi_only += 1

        cache[bp] = entry
        save_json(cache, CACHE_PATH)

    total = n_pmid + n_doi_only
    print(f"\nResolved: {total}/{len(to_resolve)} BPs  "
          f"(PMID: {n_pmid}, DOI-only: {n_doi_only})")

    if not args.dry_run:
        lit_cache = load_json(CACHE_PATH)
        lit_out = {}
        for bp, entry in lit_cache.items():
            lit_out[bp] = {k: entry.get(k, "") for k in _LIT_FIELDS}
        LIT_OUT.write_text(json.dumps(lit_out, indent=2))
        print(f"Written: {LIT_OUT}")


if __name__ == "__main__":
    main()
