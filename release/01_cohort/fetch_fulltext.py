#!/usr/bin/env python3
"""Fetch open-access full text for the papers that passed the host screen.

    python release/01_cohort/fetch_fulltext.py
    python release/01_cohort/fetch_fulltext.py --all

Reads  curation/host_screen_papers.tsv   on_target = 1 by default, every DOI with --all
Writes raw/fulltext/<doi>.xml            the article, where it is open access
       raw/raw_fulltext.tsv              one row per paper: pmcid, status, bytes

The manifest is the point. Whether full text can be retrieved is itself a screening
criterion, so the per-paper outcome has to be a recorded fact rather than the implicit
contents of a directory. Four statuses:

    retrieved        full text saved
    not_open_access  indexed, but the body is not available to fetch
    no_pmcid         no PubMed Central record, so no full text to ask for
    error            the request failed; re-run and only these are retried

The XML bodies are gitignored, as the older studies/ tree is: the repository holds code,
small tables and manifests, never the publishers' content.
"""

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCREEN = HERE / "curation" / "host_screen_papers.tsv"
BODIES = HERE / "raw" / "fulltext"
OUT = HERE / "raw" / "raw_fulltext.tsv"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
BATCH = 25
PAUSE = 0.6
COLS = ["doi", "pmcid", "status", "bytes", "file"]


def slug(doi):
    return re.sub(r"[^A-Za-z0-9._-]", "_", doi)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ids_for(dois):
    """doi -> (pmcid, is_open_access), one request per 25 DOIs."""
    query = " OR ".join(f'DOI:"{d}"' for d in dois)
    url = EPMC + "/search?" + urllib.parse.urlencode(
        {"query": query, "resultType": "lite", "format": "json", "pageSize": str(BATCH)})
    out = {}
    with urllib.request.urlopen(url, timeout=60) as fh:
        for r in json.load(fh).get("resultList", {}).get("result", []):
            doi = (r.get("doi") or "").strip().lower()
            out[doi] = ((r.get("pmcid") or "").strip(),
                        r.get("isOpenAccess") == "Y")
    return out


def body(pmcid):
    with urllib.request.urlopen(f"{EPMC}/{pmcid}/fullTextXML", timeout=90) as fh:
        return fh.read()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true",
                    help="every paper, not only the ones that passed the host screen")
    args = ap.parse_args()

    rows = read_tsv(SCREEN)
    dois = [r["doi"] for r in rows
            if r["doi"] and (args.all or r["on_target"] == "1")]
    print(f"{len(dois)} papers to fetch full text for")

    # Keep what is already on disk, so a re-run only chases errors and new DOIs.
    prior = {r["doi"]: r for r in read_tsv(OUT)} if OUT.exists() else {}
    BODIES.mkdir(parents=True, exist_ok=True)

    ids = {}
    for i in range(0, len(dois), BATCH):
        chunk = [d for d in dois[i:i + BATCH]
                 if prior.get(d, {}).get("status") not in ("retrieved", "not_open_access",
                                                           "no_pmcid")]
        if chunk:
            try:
                ids.update(ids_for(chunk))
            except Exception as e:
                print(f"  lookup batch {i // BATCH + 1}: {type(e).__name__} {e}")
            time.sleep(PAUSE)
        print(f"  ids {min(i + BATCH, len(dois))}/{len(dois)}")

    out, fetched = [], 0
    for d in dois:
        keep = prior.get(d)
        if keep and keep.get("status") in ("retrieved", "not_open_access", "no_pmcid"):
            out.append({k: keep.get(k, "") for k in COLS})
            continue
        pmcid, oa = ids.get(d, ("", False))
        path, status, size = "", "", 0
        if not pmcid:
            status = "no_pmcid"
        elif not oa:
            status = "not_open_access"
        else:
            try:
                data = body(pmcid)
                path = f"{slug(d)}.xml"
                (BODIES / path).write_bytes(data)
                status, size = "retrieved", len(data)
                fetched += 1
            except Exception:
                status = "error"
            time.sleep(PAUSE)
        out.append({"doi": d, "pmcid": pmcid, "status": status,
                    "bytes": size, "file": path})

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    import collections
    tally = collections.Counter(r["status"] for r in out)
    mb = sum(int(r["bytes"] or 0) for r in out) / 1e6
    print(f"-> {OUT.name}  {len(out)} papers, {fetched} fetched this run, {mb:.1f} MB on disk")
    for s in ("retrieved", "not_open_access", "no_pmcid", "error"):
        print(f"   {s:<16}{tally[s]:>5}")


if __name__ == "__main__":
    main()
