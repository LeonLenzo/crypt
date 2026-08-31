#!/usr/bin/env python3
"""
metadata/meta_text.py  —  Full manuscript retrieval: PMC OA text, Unpaywall +
PDF extraction, and manual-PDF ingestion for whatever's left.

Job: given a DOI (and pmcid, if meta_search.py already resolved one), get the
WHOLE paper, not just methods. Tries PMC OA full-text XML first when a pmcid
is known (free, complete, no PDF-parsing artifacts), then Unpaywall -> PDF ->
pdfminer extraction. Whatever's still unresolved after --fetch is exhausted
gets written to failed_dois.tsv for manual download; --ingest-manual then
matches hand-downloaded PDFs back to their DOI and ingests them the same way.
Stores results in a JSONL cache (one record per unique DOI) and writes
full_text back into bioprojects.json via --apply.

Modes
-----
  --fetch          Query Unpaywall + download + extract  (resumable; skips cached DOIs)
  --ingest-manual  Match hand-downloaded PDFs (manual_pdfs/) to failed_dois.tsv rows,
                   report coverage, ingest matched text into the cache
  --apply          Write full_text from cache into bioprojects.json
  --report         Coverage summary without writing anything

Options
-------
  --workers N       Parallel download workers, --fetch only (default 8)
  --limit N         Process at most N BPs, --fetch only (useful for testing)
  --max-chars       Maximum characters to store per paper (default 60000 ≈ 10k tokens)
  --write-retrieved --ingest-manual only: write 'x' back into failed_dois.tsv
                    'retrieved' column for newly-matched rows

Run from crypt/:
  python metadata/meta_text.py --fetch
  python metadata/meta_text.py --fetch --limit 10
  python metadata/meta_text.py --ingest-manual
  python metadata/meta_text.py --ingest-manual --write-retrieved
  python metadata/meta_text.py --apply
  python metadata/meta_text.py --report
"""

import argparse
import csv
import difflib
import io
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from _util import load_json, save_json, _Tee

try:
    from pdfminer.high_level import extract_text as pdf_extract
except ImportError:
    sys.exit("pip install pdfminer.six")

BIOPROJECTS    = Path("metadata/output/meta_search/data/bioprojects.json")
OUT_DIR        = Path("metadata/output/meta_text")
CACHE_PATH     = OUT_DIR / "data" / "text_cache.jsonl"
LOG_DIR        = OUT_DIR / "logs"
FAILED_DOIS    = OUT_DIR / "data" / "failed_dois.tsv"
MANUAL_PDF_DIR = OUT_DIR / "data" / "manual_pdfs"

UNPAYWALL_EMAIL = "leon.lenzo@curtin.edu.au"
UA = "meta_text/1.0 (leon.lenzo@curtin.edu.au)"

DEFAULT_WORKERS  = 8
DEFAULT_MAX_CHARS = 60_000   # ~10k tokens; well within gpt-4o-mini 128k context


# ── Cache helpers ─────────────────────────────────────────────────────────────

def load_cache() -> dict[str, dict]:
    """Return {doi: record} from JSONL cache."""
    cache: dict[str, dict] = {}
    if not CACHE_PATH.exists():
        return cache
    with open(CACHE_PATH) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            doi, _, rest = line.partition("\t")
            try:
                cache[doi] = json.loads(rest)
            except json.JSONDecodeError:
                pass
    return cache


def append_cache(doi: str, record: dict) -> None:
    """Append one record to JSONL cache (atomic-safe: appends, never rewrites)."""
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_PATH, "a") as f:
        f.write(f"{doi}\t{json.dumps(record)}\n")


# ── Unpaywall ─────────────────────────────────────────────────────────────────

def unpaywall_lookup(doi: str) -> dict:
    """Return Unpaywall metadata dict for a DOI, or {} on error."""
    url = (f"https://api.unpaywall.org/v2/"
           f"{urllib.parse.quote(doi, safe='/:')}"
           f"?email={UNPAYWALL_EMAIL}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=10) as r:
            return json.loads(r.read())
    except Exception:
        return {}


# ── PDF download + extraction ─────────────────────────────────────────────────

_PDF_HEADERS = {
    "User-Agent": UA,
    "Accept": "application/pdf,*/*",
}
MAX_PDF_BYTES = 25 * 1024 * 1024   # 25 MB hard cap

# Hosts known to serve HTML landing pages even with pdf Accept header
_REDIRECT_HOSTS = {"doi.org", "dx.doi.org"}


def _is_direct_pdf_url(url: str) -> bool:
    """Return True if url is likely a direct PDF (not a doi.org redirect)."""
    host = urllib.parse.urlparse(url).netloc.lower().lstrip("www.")
    return host not in _REDIRECT_HOSTS


def _pdf_candidate_urls(doi: str, up: dict) -> list[str]:
    """
    Return ordered list of direct PDF URLs to try, best first.

    Priority: repository PDF URLs, then publisher PDF URLs.
    doi.org redirect URLs are excluded (they serve HTML landing pages).
    """
    if not up.get("is_oa"):
        return []

    locs = up.get("oa_locations") or []
    if up.get("best_oa_location"):
        locs = [up["best_oa_location"]] + [
            l for l in locs
            if l.get("url_for_pdf") != up["best_oa_location"].get("url_for_pdf")
        ]

    repo_pdfs, pub_pdfs = [], []
    for loc in locs:
        pdf_url = loc.get("url_for_pdf") or ""
        if not pdf_url or not _is_direct_pdf_url(pdf_url):
            continue
        if loc.get("host_type") == "repository":
            repo_pdfs.append(pdf_url)
        else:
            pub_pdfs.append(pdf_url)

    candidates = repo_pdfs + pub_pdfs

    # bioRxiv canonical fallback: 10.1101/YYYY.MM.DD.XXXXXX
    if doi.startswith("10.1101/") and not any("biorxiv" in u for u in candidates):
        candidates.append(f"https://www.biorxiv.org/content/{doi}v1.full.pdf")

    return candidates


def _europepmc_pmcid(doi: str) -> str:
    """Return PMC ID (e.g. 'PMC1234567') for a DOI via Europe PMC, or ''."""
    url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           f"?query=DOI:{urllib.parse.quote(doi)}&format=json&resulttype=lite")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=10) as r:
            data = json.loads(r.read())
        for result in (data.get("resultList") or {}).get("result") or []:
            pmcid = result.get("pmcid") or ""
            if pmcid:
                return pmcid
    except Exception:
        pass
    return ""


def _fetch_epmc_xml_text(pmcid: str, max_chars: int) -> str:
    """Fetch Europe PMC full-text XML and return plain-text body, or ''."""
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=20) as r:
            xml = r.read().decode("utf-8", errors="replace")
        # Strip XML tags; keep whitespace structure
        text = re.sub(r"<[^>]+>", " ", xml)
        text = re.sub(r"\s{2,}", " ", text)
        text = re.sub(r"\n{3,}", "\n\n", text)
        return text[:max_chars].strip()
    except Exception:
        return ""


def download_pdf(url: str, timeout: int = 45) -> bytes | None:
    """Download a PDF URL; return raw bytes or None on failure."""
    try:
        req = urllib.request.Request(url, headers=_PDF_HEADERS)
        with urllib.request.urlopen(req, timeout=timeout) as r:
            ct = r.headers.get("Content-Type", "")
            # Reject HTML landing pages (paywall redirects etc.)
            if "html" in ct.lower() and "pdf" not in ct.lower():
                return None
            # Read in chunks up to max size
            chunks = []
            total  = 0
            while True:
                chunk = r.read(65536)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_PDF_BYTES:
                    return None   # too large
                chunks.append(chunk)
            return b"".join(chunks)
    except Exception:
        return None


def extract_text_full(pdf_bytes: bytes) -> str:
    """Extract plain text from PDF bytes via pdfminer, no truncation."""
    try:
        text = pdf_extract(io.BytesIO(pdf_bytes))
        text = re.sub(r"\n{3,}", "\n\n", text)
        text = re.sub(r" {2,}", " ", text)
        return text.strip()
    except Exception:
        return ""


def extract_text(pdf_bytes: bytes, max_chars: int) -> str:
    """Extract plain text from PDF bytes via pdfminer; truncate to max_chars."""
    return extract_text_full(pdf_bytes)[:max_chars].strip()


# ── Per-DOI worker ────────────────────────────────────────────────────────────

def process_doi(doi: str, max_chars: int, pmcid: str = "") -> dict:
    """
    Full pipeline for one DOI.

    Strategy (in order):
      1. Europe PMC full-text XML — if pmcid is already known (meta_search.py
         resolves this during BioProject identification), this is free, complete,
         and has none of pdfminer's PDF-parsing artifacts. Tried first.
      2. Unpaywall oa_locations — prefer repository PDF URLs over publisher
      3. bioRxiv canonical PDF pattern (if 10.1101/ DOI)
      4. Europe PMC full-text XML again, via a deferred lookup — covers DOIs
         meta_search didn't have a pmcid for yet, only if the PDF path failed
    """
    record: dict = {"doi": doi, "oa_status": "closed",
                    "pdf_url": "", "full_text": "", "error": ""}

    # Step 1: PMC OA full text, if we already know the pmcid (whole manuscript,
    # not just methods — meta_search.py never writes methods_text anymore)
    if pmcid:
        text = _fetch_epmc_xml_text(pmcid, max_chars)
        if text:
            record["oa_status"] = "pmc_oa"
            record["pdf_url"]   = f"epmc:{pmcid}"
            record["full_text"] = text
            return record

    up = unpaywall_lookup(doi)
    record["oa_status"] = up.get("oa_status", "closed") if up else "closed"

    # Step 2 + 3: try direct PDF URLs (repository first, then publisher)
    pdf_urls = _pdf_candidate_urls(doi, up)
    for url in pdf_urls:
        pdf_bytes = download_pdf(url)
        if pdf_bytes is None:
            continue
        text = extract_text(pdf_bytes, max_chars)
        if text:
            record["pdf_url"]   = url
            record["full_text"] = text
            return record

    # Step 4: Europe PMC full-text XML (deferred lookup — only if PDF failed
    # and meta_search hadn't already found a pmcid for us in step 1)
    if not pmcid and (up.get("is_oa") or doi.startswith("10.1101/")):
        found_pmcid = _europepmc_pmcid(doi)
        if found_pmcid:
            text = _fetch_epmc_xml_text(found_pmcid, max_chars)
            if text:
                record["pdf_url"]   = f"epmc:{found_pmcid}"
                record["full_text"] = text
                return record

    if not up.get("is_oa"):
        record["error"] = "no_oa_pdf"
    else:
        record["error"] = "download_failed"
    return record


# ── --fetch ───────────────────────────────────────────────────────────────────

def cmd_fetch(workers: int, limit: int, max_chars: int) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    sys.stdout = _Tee(LOG_DIR / "meta_text.log")

    bps   = load_json(BIOPROJECTS)
    cache = load_cache()

    # Collect unique DOIs not yet cached, carrying the pmcid meta_search.py
    # already resolved (if any) so process_doi can skip straight to PMC OA text
    seen: dict[str, tuple[str, str]] = {}   # doi -> (first bp seen, pmcid)
    for bp, v in bps.items():
        doi = v.get("doi")
        if doi and doi not in seen and doi not in cache:
            seen[doi] = (bp, v.get("pmcid", ""))

    total_doi = len(seen)
    if limit:
        seen = dict(list(seen.items())[:limit])

    already = len(cache)
    print(f"meta_text  --fetch")
    print(f"  BioProjects: {len(bps):,}  unique DOIs with data: "
          f"{sum(1 for v in bps.values() if v.get('doi')):,}")
    print(f"  Cache: {already:,} DOIs already processed")
    print(f"  To fetch: {len(seen):,}"
          + (f"  (limited to {limit})" if limit else ""))
    print(f"  Workers: {workers}   max_chars: {max_chars:,}")
    print()

    n_ok = n_fail = n_closed = 0

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process_doi, doi, max_chars, pmcid): doi
                for doi, (_bp, pmcid) in seen.items()}
        for i, fut in enumerate(as_completed(futs), 1):
            doi    = futs[fut]
            record = fut.result()
            append_cache(doi, record)

            ft = record.get("full_text", "")
            err = record.get("error", "")
            if ft:
                n_ok += 1
            elif record.get("oa_status") in ("closed", ""):
                n_closed += 1
            else:
                n_fail += 1

            if i % 20 == 0 or i == len(seen):
                print(f"  {i:4d}/{len(seen)}  "
                      f"ok={n_ok}  closed={n_closed}  failed={n_fail}")

    print()
    print(f"Done.  full_text retrieved: {n_ok:,}  "
          f"closed: {n_closed:,}  other failures: {n_fail:,}")
    print(f"Cache: {CACHE_PATH}  ({already + len(seen):,} total entries)")

    reloaded_cache = load_cache()  # pick up the entries just appended
    n_new_failed = sync_failed_dois(bps, reloaded_cache)
    if n_new_failed:
        print(f"failed_dois.tsv: {n_new_failed:,} newly-failed DOIs appended for manual pull")
    print()
    print("Run  python metadata/meta_text.py --apply  to write into bioprojects.json")


# ── --ingest-manual ──────────────────────────────────────────────────────────

def load_failed_dois() -> list[dict]:
    with open(FAILED_DOIS, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save_failed_dois(rows: list[dict]) -> None:
    fieldnames = ["retrieved", "oa_status", "doi", "bioproject", "title"]
    with open(FAILED_DOIS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def sync_failed_dois(bps: dict, cache: dict) -> int:
    """Append any DOI that failed automated retrieval and isn't already
    tracked in failed_dois.tsv. Success/failure is read from `cache` (the
    true source of truth — bps' own full_text field is only populated after
    --apply runs, so it lags cache right after a --fetch pass). Without this,
    DOIs that meta_search.py newly resolves but that turn out to be paywalled
    silently fail to surface for manual pull — failed_dois.tsv only reflected
    whatever was true the day it was hand-seeded, not the current cache state.
    Returns the number of new rows appended."""
    rows  = load_failed_dois() if FAILED_DOIS.exists() else []
    known = {r["doi"] for r in rows}

    seen: set[str] = set()
    new_rows: list[dict] = []
    for bp, v in bps.items():
        doi = v.get("doi")
        if not doi or doi in seen or doi in known:
            continue
        seen.add(doi)
        rec = cache.get(doi)
        if not rec or rec.get("full_text"):
            continue  # not attempted yet, or already succeeded — not a failure
        new_rows.append({
            "retrieved": "", "oa_status": rec.get("oa_status", ""),
            "doi": doi, "bioproject": bp, "title": v.get("title", ""),
        })

    if new_rows:
        save_failed_dois(rows + new_rows)
    return len(new_rows)


def _norm(s: str) -> str:
    """Lowercase, strip to alnum+space, collapse whitespace — for fuzzy compare."""
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", s.lower())).strip()


# pdfminer sometimes inserts stray whitespace inside strings because of PDF
# kerning/ligature tables (e.g. a DOI renders as "10. 11 11/mp p .12 464").
# Match DOIs with optional whitespace between every character rather than
# a plain substring check.
def _doi_pattern(doi: str) -> re.Pattern:
    return re.compile(r"\s*".join(re.escape(c) for c in doi), re.IGNORECASE)


HEADER_CHARS = 6000

# A paper's OWN doi is almost always followed shortly by article-metadata
# boilerplate (ScienceDirect/Elsevier "Received ... Accepted ..." lines,
# copyright notices); a doi found only in a reference-list citation isn't.
# Used to disambiguate when a target doi shows up more than once in the body
# (e.g. a paper cites another paper that's also in our failed-doi list).
_OWN_DOI_MARKERS = ("received", "accepted", "revised", "available online",
                     "all rights reserved")


def _looks_like_own_doi(text: str, match_end: int) -> bool:
    tail = text[match_end:match_end + 80].lower()
    return any(m in tail for m in _OWN_DOI_MARKERS)


TITLE_FUZZY_THRESHOLD = 0.5
TITLE_FUZZY_MARGIN = 0.1   # best must beat second-best by this much

# doi -> manual_pdfs filename, for cases the automated matcher can't resolve
# (e.g. a bioRxiv preprint DOI whose only local copy is the later published-
# journal PDF, which never prints the preprint DOI anywhere in its text).
MANUAL_OVERRIDES = {
    # PRJNA1207595 / "Phytophthora capsici multiple host RNA-Seq" — Unpaywall
    # resolved the bioRxiv preprint DOI, but the file on hand is the later
    # Microbial Pathogenesis (2026) publication of the same study.
    "10.1101/2025.05.26.656185": "main_a.pdf",
}


def phase1_exact_matches(pdf_text: dict[Path, str], rows: list[dict]):
    """
    First pass: DOI-string-in-text matching across ALL pdfs before any fuzzy
    matching runs, so a low-confidence fuzzy hit can never block a later
    high-confidence exact hit (order-independent).

    Returns (claimed, dup_pdfs, ambiguous_pdfs, leftover_pdfs):
      claimed         doi -> (pdf, method)
      dup_pdfs        [(pdf, doi, method)] extra pdfs exact-matching a doi
                       that already has a winner (kept for reporting only)
      ambiguous_pdfs  [(pdf, reason)] multiple different target DOIs found
                       in body with none in the header — can't disambiguate
      leftover_pdfs   pdfs with zero exact hits — candidates for phase 2
    """
    patterns = {r["doi"]: _doi_pattern(r["doi"]) for r in rows}

    candidates: dict[str, list[tuple[Path, str]]] = {}
    ambiguous_pdfs = []
    leftover_pdfs = []

    for pdf, text in pdf_text.items():
        head, body = text[:HEADER_CHARS], text
        header_hits = [doi for doi, pat in patterns.items() if pat.search(head)]
        if header_hits:
            best = min(header_hits, key=lambda d: patterns[d].search(head).start())
            candidates.setdefault(best, []).append((pdf, "doi_in_header"))
            continue
        body_matches = {doi: m for doi, pat in patterns.items()
                         if (m := pat.search(body))}
        if len(body_matches) == 1:
            doi = next(iter(body_matches))
            candidates.setdefault(doi, []).append((pdf, "doi_in_body"))
        elif len(body_matches) > 1:
            own = [doi for doi, m in body_matches.items()
                   if _looks_like_own_doi(body, m.end())]
            if len(own) == 1:
                candidates.setdefault(own[0], []).append((pdf, "doi_in_body_context"))
            else:
                ambiguous_pdfs.append((pdf, f"ambiguous_doi_hits({len(body_matches)})"))
        else:
            leftover_pdfs.append(pdf)

    claimed: dict[str, tuple[Path, str]] = {}
    dup_pdfs = []
    for doi, hits in candidates.items():
        hits.sort(key=lambda h: h[0].name)
        claimed[doi] = hits[0]
        for pdf, method in hits[1:]:
            dup_pdfs.append((pdf, doi, method))

    return claimed, dup_pdfs, ambiguous_pdfs, leftover_pdfs


def phase2_fuzzy_matches(leftover_pdfs: list[Path], pdf_text: dict[Path, str],
                          rows: list[dict], claimed: dict[str, tuple[Path, str]]):
    """
    Second pass: fuzzy title match, restricted to pdfs/rows that survived
    phase 1 unclaimed. Requires a score threshold AND a margin over the
    second-best candidate, so a generic bioproject-style title (e.g.
    "Colletotrichum fructicola RNA sequencing") can't confidently steal a
    doi it merely resembles.

    Mutates `claimed` in place. Returns [(pdf, reason)] for pdfs still unmatched.
    """
    unclaimed_rows = [r for r in rows if r["doi"] not in claimed]
    still_unmatched = []

    for pdf in leftover_pdfs:
        text = pdf_text[pdf]
        text_head_norm = _norm(text[:600])
        fname_norm = _norm(pdf.stem)

        scored = []
        for row in unclaimed_rows:
            title_norm = _norm(row["title"])
            if not title_norm:
                continue
            score = max(
                difflib.SequenceMatcher(None, title_norm, text_head_norm).ratio(),
                difflib.SequenceMatcher(None, title_norm, fname_norm).ratio(),
            )
            scored.append((score, row["doi"]))
        scored.sort(reverse=True)

        if not scored:
            still_unmatched.append((pdf, "no_match"))
            continue
        best_score, best_doi = scored[0]
        second_score = scored[1][0] if len(scored) > 1 else 0.0
        if best_score >= TITLE_FUZZY_THRESHOLD and (best_score - second_score) >= TITLE_FUZZY_MARGIN:
            claimed[best_doi] = (pdf, f"title_fuzzy({best_score:.2f})")
            unclaimed_rows = [r for r in unclaimed_rows if r["doi"] != best_doi]
        else:
            still_unmatched.append(
                (pdf, f"no_match(best={best_score:.2f} 2nd={second_score:.2f})"))

    return still_unmatched


def cmd_ingest_manual(max_chars: int, write_retrieved: bool) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    sys.stdout = _Tee(LOG_DIR / "meta_text_manual.log")

    if not FAILED_DOIS.exists():
        raise SystemExit(f"Not found: {FAILED_DOIS}")
    if not MANUAL_PDF_DIR.exists():
        raise SystemExit(f"Not found: {MANUAL_PDF_DIR}")

    rows = load_failed_dois()
    cache = load_cache()
    pdf_files = sorted(MANUAL_PDF_DIR.glob("*.pdf"))

    print(f"meta_text --ingest-manual")
    print(f"  failed_dois.tsv: {len(rows)} rows")
    print(f"  manual_pdfs/:    {len(pdf_files)} PDFs")
    print()

    pdf_text: dict[Path, str] = {}
    unmatched_pdfs: list[tuple[Path, str]] = []

    for pdf in pdf_files:
        try:
            raw = pdf.read_bytes()
        except Exception as e:
            unmatched_pdfs.append((pdf, f"read_error:{e}"))
            continue
        text = extract_text_full(raw)
        if not text:
            unmatched_pdfs.append((pdf, "extract_failed"))
            continue
        pdf_text[pdf] = text

    claimed, dup_pdfs, ambiguous_pdfs, leftover_pdfs = phase1_exact_matches(pdf_text, rows)
    still_unmatched = phase2_fuzzy_matches(leftover_pdfs, pdf_text, rows, claimed)
    unmatched_pdfs += ambiguous_pdfs + still_unmatched

    row_by_doi = {r["doi"]: r for r in rows}
    for doi, fname in MANUAL_OVERRIDES.items():
        pdf = MANUAL_PDF_DIR / fname
        if doi in claimed or pdf not in pdf_text:
            continue
        claimed[doi] = (pdf, "manual_override")
        unmatched_pdfs = [(p, r) for p, r in unmatched_pdfs if p != pdf]

    # claimed values are now (pdf, method) — attach text for the ingest step below
    claimed = {doi: (pdf, method, pdf_text[pdf]) for doi, (pdf, method) in claimed.items()}
    matched_dois = set(claimed)
    marked_x = {r["doi"] for r in rows if r["retrieved"].strip().lower() == "x"}
    missing = [r for r in rows if r["doi"] not in matched_dois]

    print(f"  All matches ({len(claimed)}):")
    for doi, (pdf, method, text) in sorted(claimed.items(), key=lambda kv: kv[1][1]):
        print(f"    {method:22s} {doi:32s} <- {pdf.name}")
    print()

    print(f"  Matched:  {len(matched_dois)}/{len(rows)} failed-DOI rows have a PDF")
    print(f"  Marked 'x' in tsv: {len(marked_x)}/{len(rows)}")
    print(f"  Marked 'x' but NO PDF matched: "
          f"{len(marked_x - matched_dois)}")
    print(f"  PDF matched but NOT marked 'x': "
          f"{len(matched_dois - marked_x)}")
    print()

    if dup_pdfs:
        print(f"  Duplicate PDFs for an already-matched DOI ({len(dup_pdfs)}):")
        for pdf, doi, method in dup_pdfs:
            winner = claimed[doi][0].name if doi in claimed else "?"
            print(f"    {pdf.name}  ({method} -> {doi}, kept {winner})")
        print()

    if unmatched_pdfs:
        print(f"  Unmatched PDFs — no confident DOI/title match ({len(unmatched_pdfs)}):")
        for pdf, reason in unmatched_pdfs:
            print(f"    {pdf.name}  ({reason})")
        print()

    if missing:
        print(f"  failed_dois.tsv rows with NO manual PDF ({len(missing)}):")
        for r in missing:
            tag = r["retrieved"].strip() or "blank"
            print(f"    {r['doi']}  [{tag}]  {r['title'][:70]}")
        print()

    # Ingest matched PDFs into text_cache.jsonl (resumable — skip already cached)
    n_ingested = n_already = 0
    for doi, (pdf, method, text) in claimed.items():
        if doi in cache and cache[doi].get("full_text"):
            n_already += 1
            continue
        row = row_by_doi[doi]
        record = {
            "doi": doi,
            "oa_status": row.get("oa_status", ""),
            "pdf_url": f"manual:{pdf.name}",
            "full_text": text[:max_chars],
            "error": "",
        }
        append_cache(doi, record)
        n_ingested += 1

    print(f"  Ingested into text_cache.jsonl: {n_ingested}")
    print(f"  Already cached (skipped):       {n_already}")
    print()

    if write_retrieved:
        n_marked = 0
        for r in rows:
            if r["doi"] in matched_dois and r["retrieved"].strip().lower() != "x":
                r["retrieved"] = "x"
                n_marked += 1
        save_failed_dois(rows)
        print(f"  Updated failed_dois.tsv: {n_marked} rows newly marked 'x'")
    else:
        n_would = len(matched_dois - marked_x)
        if n_would:
            print(f"  {n_would} rows would be marked 'x' — rerun with "
                  f"--write-retrieved to apply")


# ── --apply ───────────────────────────────────────────────────────────────────

def cmd_apply() -> None:
    cache = load_cache()
    bps   = load_json(BIOPROJECTS)

    if not cache:
        raise SystemExit(f"Cache empty — run --fetch first")

    # Build doi → full_text map (use last non-empty entry per DOI)
    doi_text: dict[str, str] = {}
    for doi, rec in cache.items():
        if rec.get("full_text"):
            doi_text[doi] = rec["full_text"]

    n_updated = n_already = n_no_text = 0
    for bp, v in bps.items():
        doi = v.get("doi")
        if not doi:
            continue
        if v.get("full_text"):
            n_already += 1
            continue
        text = doi_text.get(doi)
        if text:
            v["full_text"] = text
            n_updated += 1
        else:
            n_no_text += 1

    print(f"--apply")
    print(f"  Cache entries with full_text: {len(doi_text):,}")
    print(f"  Updated:        {n_updated:,}")
    print(f"  Already filled: {n_already:,}")
    print(f"  No text in cache: {n_no_text:,}")

    if n_updated:
        save_json(bps, BIOPROJECTS)
        print(f"  Written: {BIOPROJECTS}")
    else:
        print("  Nothing to write.")


# ── --report ──────────────────────────────────────────────────────────────────

def cmd_report() -> None:
    bps   = load_json(BIOPROJECTS)
    cache = load_cache()

    total     = len(bps)
    has_doi   = sum(1 for v in bps.values() if v.get("doi"))
    has_ft    = sum(1 for v in bps.values() if v.get("full_text"))
    no_doi    = total - has_doi

    from collections import Counter
    status_counts = Counter(r.get("oa_status", "?") for r in cache.values())
    err_counts    = Counter(r.get("error", "") for r in cache.values() if r.get("error"))
    cached_ok     = sum(1 for r in cache.values() if r.get("full_text"))

    print(f"meta_text report  ({total:,} BioProjects)")
    print(f"  With DOI:            {has_doi:,} ({100*has_doi/total:.1f}%)")
    print(f"  No DOI:              {no_doi:,} ({100*no_doi/total:.1f}%)")
    print(f"  full_text in BPs:    {has_ft:,} ({100*has_ft/total:.1f}%)")
    print()
    print(f"  Cache: {len(cache):,} DOIs  —  {cached_ok:,} with full_text")
    print(f"  OA status: " + "  ".join(f"{s}={n}" for s, n in status_counts.most_common()))
    if err_counts:
        print(f"  Errors:    " + "  ".join(f"{e}={n}" for e, n in err_counts.most_common()))


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--fetch",  action="store_true",
                      help="query Unpaywall + download PDFs + extract text")
    mode.add_argument("--ingest-manual", action="store_true",
                      help="match manual_pdfs/ to failed_dois.tsv, report "
                           "coverage, ingest matched text into the cache")
    mode.add_argument("--apply",  action="store_true",
                      help="write full_text from cache into bioprojects.json")
    mode.add_argument("--report", action="store_true",
                      help="show coverage stats")
    ap.add_argument("--workers",   type=int, default=DEFAULT_WORKERS)
    ap.add_argument("--limit",     type=int, default=0,
                    help="process at most N BPs (0 = all)")
    ap.add_argument("--max-chars", type=int, default=DEFAULT_MAX_CHARS,
                    dest="max_chars")
    ap.add_argument("--write-retrieved", action="store_true",
                    dest="write_retrieved",
                    help="--ingest-manual only: write 'x' into failed_dois.tsv "
                         "'retrieved' column for newly-matched rows")
    args = ap.parse_args()

    if args.fetch:
        cmd_fetch(args.workers, args.limit, args.max_chars)
    elif args.ingest_manual:
        cmd_ingest_manual(args.max_chars, args.write_retrieved)
    elif args.apply:
        cmd_apply()
    elif args.report:
        cmd_report()


if __name__ == "__main__":
    main()
