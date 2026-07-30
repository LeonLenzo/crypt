#!/usr/bin/env python3
"""
03a_fetch_xml.py — fetch BioProject and BioSample XML attributes.

Reads runs.tsv for BioProject and BioSample accessions.
BioProject data bootstrapped from existing bioprojects.json output (no extra API
calls for already-fetched entries); only new BioProjects are fetched from NCBI.
BioSample attributes fetched in batches from NCBI BioSample efetch.

Output:
  output/03a_fetch_xml/data/bioprojects.json  — title, description, submission_date
  output/03a_fetch_xml/data/biosamples.json   — tissue, geo_loc_name, collection_date,
                                                lat_lon, host, isolation_source, dev_stage

Run from crypt/:
  python 03a_fetch_xml.py
"""

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from _util import _Tee, http_get, link_latest, load_json, make_log_dir, save_json

RUNS_TSV        = Path("output/02_filter_runs/data/runs.tsv")
OUT_DIR         = Path("output/03a_fetch_xml")
BIOSAMPLE_CACHE = OUT_DIR / "data" / "biosample_cache.json"

API_KEY  = os.environ.get("NCBI_API_KEY", "")
RATE     = 9.0 if API_KEY else 2.5
ENTREZ   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS  = {"User-Agent": "crypt/03a_fetch_xml (leon.lenzo@curtin.edu.au)"}

BS_BATCH = 300   # BioSamples per efetch request

_HARMONIZED = {
    "tissue", "geo_loc_name", "collection_date", "lat_lon",
    "host", "isolation_source", "dev_stage",
}

_last_req: float = 0.0


def _wait() -> None:
    global _last_req
    gap  = 1.0 / RATE
    wait = _last_req + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_req = time.monotonic()


def _get(path: str, **params) -> bytes:
    if API_KEY:
        params["api_key"] = API_KEY
    _wait()
    return http_get(f"{ENTREZ}/{path}?{urllib.parse.urlencode(params)}", HEADERS)


def _post(path: str, **params) -> bytes:
    """POST request to Entrez (avoids HTTP 414 for large id lists)."""
    if API_KEY:
        params["api_key"] = API_KEY
    _wait()
    url  = f"{ENTREZ}/{path}"
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=HEADERS, method="POST")
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except Exception as exc:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


# ── BioProject XML fetch ───────────────────────────────────────────────────────

def _fetch_bp_xml(accession: str) -> dict:
    """Fetch title, description, submission_date from NCBI BioProject XML."""
    try:
        raw  = _get("esearch.fcgi", db="bioproject",
                    term=f"{accession}[Project Accession]", retmax=1)
        root = ET.fromstring(raw)
        ids  = [el.text for el in root.findall(".//Id") if el.text]
        if not ids:
            return {}
        uid = ids[0]
        raw  = _get("efetch.fcgi", db="bioproject", id=uid)
        root = ET.fromstring(raw)
        title = (root.findtext(".//ProjectDescr/Title") or "").strip()
        desc  = (root.findtext(".//ProjectDescr/Description") or "").strip()
        sub_el = root.find(".//Submission[@submitted]")
        subdate = sub_el.get("submitted", "") if sub_el is not None else ""
        return {"title": title, "description": desc, "submission_date": subdate}
    except Exception as exc:
        print(f"  WARNING: {accession}: {exc}", flush=True)
        return {}


# ── BioSample XML batch fetch ─────────────────────────────────────────────────

def _parse_biosample(elem: ET.Element) -> dict:
    """Extract harmonized attributes from one BioSample XML element."""
    acc = elem.get("accession", "")
    row: dict = {"accession": acc}
    for attr in elem.findall(".//Attribute"):
        harm = (attr.get("harmonized_name") or attr.get("attribute_name") or "")
        harm = harm.lower().strip().replace(" ", "_")
        if harm in _HARMONIZED and attr.text and attr.text.strip():
            row.setdefault(harm, attr.text.strip())
    return row


def _fetch_biosamples(accessions: list[str], cache: dict) -> dict:
    # stubs: cache entries with only the accession key (previous failed fetches)
    needed = [a for a in accessions
              if a not in cache or set(cache[a].keys()) <= {"accession"}]
    print(f"BioSamples: {len(accessions):,} total, {len(needed):,} to fetch", flush=True)
    results = dict(cache)

    for i in range(0, len(needed), BS_BATCH):
        batch = needed[i: i + BS_BATCH]
        try:
            raw  = _post("efetch.fcgi", db="biosample",
                         id=",".join(batch), rettype="xml", retmode="xml")
            root = ET.fromstring(raw)
            for bs in root.findall(".//BioSample"):
                parsed = _parse_biosample(bs)
                acc    = parsed.get("accession")
                if acc:
                    results[acc] = parsed
        except Exception as exc:
            print(f"  batch {i // BS_BATCH + 1} error: {exc}", flush=True)

        done = min(i + BS_BATCH, len(needed))
        if done % 3000 == 0 or done == len(needed):
            print(f"  {done:,} / {len(needed):,}", flush=True)
            save_json(results, BIOSAMPLE_CACHE)

    for acc in needed:
        results.setdefault(acc, {"accession": acc})

    save_json(results, BIOSAMPLE_CACHE)
    return results


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log       = _Tee(log_dir / "fetch.log")
    link_latest(logs_base, log_dir / "fetch.log")
    sys.stdout = log

    print(f"NCBI_API_KEY: {'yes' if API_KEY else 'no'}", flush=True)

    if not RUNS_TSV.exists():
        sys.exit(f"ERROR: {RUNS_TSV} not found — run 02_filter_runs.py first")

    # ── Collect accessions from runs.tsv ──────────────────────────────────────
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

    print(f"BioProjects: {len(bioprojects):,}  BioSamples: {len(biosamples):,}",
          flush=True)

    # ── BioProjects: bootstrap from own output ────────────────────────────────
    bp_path    = OUT_DIR / "data" / "bioprojects.json"
    bp_existing = load_json(bp_path) if bp_path.exists() else {}
    print(f"BioProject bootstrap from own output: {len(bp_existing):,} entries", flush=True)

    bp_out: dict[str, dict] = {}
    to_fetch: list[str] = []
    for bp in bioprojects:
        existing = bp_existing.get(bp, {})
        if existing.get("title") is not None:
            bp_out[bp] = existing
        else:
            to_fetch.append(bp)

    print(f"BioProjects: {len(bp_out):,} from cache, {len(to_fetch):,} to fetch",
          flush=True)
    for bp in to_fetch:
        result = _fetch_bp_xml(bp)
        bp_out[bp] = result
        label = result.get("title", "")[:60] or "(no title)"
        print(f"  {bp}: {label}", flush=True)

    bp_path.write_text(json.dumps(bp_out, indent=2))
    print(f"Written: {bp_path}  ({len(bp_out):,} entries)", flush=True)

    # ── BioSamples: batch efetch ──────────────────────────────────────────────
    bs_cache = load_json(BIOSAMPLE_CACHE) if BIOSAMPLE_CACHE.exists() else {}
    bs_out   = _fetch_biosamples(sorted(biosamples), bs_cache)

    # filter to only accessions in runs.tsv (NCBI efetch can return extra related entries)
    bs_out_clean = {k: v for k, v in bs_out.items() if k in biosamples}
    bs_path = OUT_DIR / "data" / "biosamples.json"
    bs_path.write_text(json.dumps(bs_out_clean, indent=2))
    print(f"Written: {bs_path}  ({len(bs_out_clean):,} entries)", flush=True)

    # ── Coverage summary ──────────────────────────────────────────────────────
    n = len(bs_out)
    for field in sorted(_HARMONIZED):
        count = sum(1 for v in bs_out.values() if v.get(field))
        print(f"  {field:<20} {count:>6,}  ({100*count/max(n,1):.1f}%)", flush=True)


if __name__ == "__main__":
    main()
