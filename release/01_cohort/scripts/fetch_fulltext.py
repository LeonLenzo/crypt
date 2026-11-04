#!/usr/bin/env python3
"""Fetch the title and the fullest available text for every DOI in data/inputs/pg_dois.txt.

    python release/01_cohort/scripts/fetch_fulltext.py
    python release/01_cohort/scripts/fetch_fulltext.py --retry     # ask again for what is missing

Reads  data/curation/screen_abstract.tsv  the gate: only papers NOT recommended `drop` are fetched
       data/inputs/fulltext_manual/*      anything obtained by hand, checked first
Writes data/raw/fulltext/<doi>.*          the text, by whichever route reached it
       data/raw/fulltext.tsv          one row per paper: title, source, note, bytes, path

`source` records the route, best first:

    manual    obtained by hand and dropped into data/inputs/fulltext_manual/. No script can
              fetch these, which is what makes them an input rather than a retrieval
    pmc       EuropePMC, as JATS XML with sections marked up
    biorxiv   the preprint server's own full-text page. EuropePMC indexes preprints as PPR
              records with an abstract but returns 404 for fullTextXML, so the PMC route
              misses every one of them
    none      no body retrievable; the abstract is kept in the `abstract` column but NOT
              written to fulltext/, which holds full text only. `note` says why

FULL TEXT IS FETCHED ONLY FOR SURVIVORS. screen_abstract.py runs first and marks reviews,
marker assays and book chapters `drop`; those are skipped here, so no reviewer hand-downloads
the PDF of a paper the abstract already ruled out. A `source=none` survivor is one to fetch by
hand: drop its PDF in data/inputs/fulltext_manual/ and convert with
research/01_literature/retrieve/pdf_to_text.py.

The file extension carries the route (.xml, .biorxiv.html, else manual), so
provenance survives the manifest being lost or its columns changing, which has happened twice.

NOTHING IS REFETCHED. A file already in data/raw/fulltext/ settles that DOI before any request is
made, checked by looking on disk rather than by trusting this manifest: the manifest's own
columns have changed under it, and reading them to decide what exists once sent 147 already
downloaded articles back to EuropePMC.
"""

import argparse
import csv
import json
import re
import shutil
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
DATA = MODULE / "data"
SCREEN = DATA / "curation" / "screen_abstract.tsv"
MANUAL = DATA / "inputs" / "fulltext_manual"
BODIES = DATA / "raw" / "fulltext"
OUT = DATA / "raw" / "fulltext.tsv"
MAILTO = "leon.lenzo@curtin.edu.au"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
CROSSREF = "https://api.crossref.org/works/"
BIORXIV = "https://www.biorxiv.org/content"
PREPRINT = re.compile(r"^10\.1101/")
UA = {"User-Agent": f"Mozilla/5.0 (compatible; crypt/1.0; +mailto:{MAILTO})"}
BATCH = 25
PAUSE = 0.6
COLS = ["doi", "title", "pmcid", "source", "note", "bytes", "path", "abstract"]


def slug(doi):
    return re.sub(r"[^A-Za-z0-9._-]", "_", doi.strip())


def clean(text):
    """Strip the markup EuropePMC and Crossref wrap text in, keeping section labels."""
    t = re.sub(r"<(?:jats:)?title[^>]*>(.*?)</(?:jats:)?title>", r"\1. ", text or "",
               flags=re.S | re.I)
    t = re.sub(r"<h4[^>]*>(.*?)</h4>", r"\1. ", t, flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    for a, b in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"), ("&quot;", '"'),
                 ("&#x2019;", "'"), ("&nbsp;", " ")):
        t = t.replace(a, b)
    return " ".join(t.split())


def source_of(path):
    n = path.name.lower()
    if n.endswith(".xml"):
        return "pmc"
    if n.endswith(".biorxiv.html"):
        return "biorxiv"
    return "manual"


def get(url, headers=None, timeout=60):
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers or {}),
                                timeout=timeout) as fh:
        return fh.read()


def from_biorxiv(doi):
    """The preprint's own full-text page. Versions are tried because the DOI carries none."""
    for v in (1, 2, 3):
        try:
            data = get(f"{BIORXIV}/{doi}v{v}.full", UA, 60)
            if len(data) > 20000:        # an error page or a stub is far smaller
                return data
        except Exception:
            continue
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--retry", action="store_true",
                    help="ask again for papers previously recorded as unavailable")
    args = ap.parse_args()

    if not SCREEN.exists():
        raise SystemExit("run screen_abstract.py first; this fetches only its survivors")
    with open(SCREEN, newline="", encoding="utf-8") as fh:
        screened = list(csv.DictReader(fh, delimiter="\t"))
    dropped = {r["doi"] for r in screened if r["recommend"] == "drop"}
    dois = [r["doi"] for r in screened if r["recommend"] != "drop"]
    print(f"{len(screened)} screened, {len(dropped)} dropped at the abstract, "
          f"{len(dois)} to fetch")
    BODIES.mkdir(parents=True, exist_ok=True)
    prior = {}
    if OUT.exists():
        with open(OUT, newline="", encoding="utf-8") as fh:
            prior = {r["doi"]: r for r in csv.DictReader(fh, delimiter="\t")}
    abs_tbl = {}
    ab = DATA / "raw" / "abstracts.tsv"
    if ab.exists():
        with open(ab, newline="", encoding="utf-8") as fh:
            abs_tbl = {r["doi"]: r for r in csv.DictReader(fh, delimiter="\t")}
    manual = {f.stem: f for f in MANUAL.iterdir() if f.is_file()} if MANUAL.is_dir() else {}
    ondisk = {}
    for f in BODIES.iterdir():
        if f.is_file() and f.stat().st_size:
            ondisk.setdefault(re.sub(r"\.(biorxiv\.html|xml|txt|html)$", "",
                                     f.name, flags=re.I), f)
    print(f"{len(dois)} DOIs; {len(manual)} hand-obtained waiting; "
          f"{sum(1 for d in dois if slug(d) in ondisk)} already have text")

    rows, fetched = [], 0
    for d in dois:
        a = abs_tbl.get(d, {})
        title, abstract = a.get("title", ""), a.get("abstract", "")
        pmcid, oa = a.get("pmcid", ""), a.get("open_access") == "1"

        # Hand-obtained text wins: it exists because nothing could fetch it.
        if slug(d) in manual:
            dst = BODIES / manual[slug(d)].name
            if not dst.exists():
                shutil.move(str(manual[slug(d)]), str(dst))
            ondisk[slug(d)] = dst

        f = ondisk.get(slug(d))
        if f:
            rows.append({"doi": d, "title": title, "pmcid": pmcid, "source": source_of(f),
                         "note": "", "bytes": f.stat().st_size,
                         "path": str(f.relative_to(MODULE)), "abstract": abstract})
            continue

        path, source, note, size = "", "none", "", 0
        if pmcid and oa:
            try:
                data = get(f"{EPMC}/{pmcid}/fullTextXML", timeout=90)
                p = BODIES / f"{slug(d)}.xml"
                p.write_bytes(data)
                source, size, path, fetched = "pmc", len(data), str(p.relative_to(MODULE)), fetched + 1
            except Exception:
                note = "full text request failed"
            time.sleep(PAUSE)
        elif pmcid:
            note = "not open access"
        else:
            note = "no PMC record"

        if source == "none" and PREPRINT.match(d):
            data = from_biorxiv(d)
            if data:
                p = BODIES / f"{slug(d)}.biorxiv.html"
                p.write_bytes(data)
                source, size, note, path = "biorxiv", len(data), "", str(p.relative_to(MODULE))
                fetched += 1
            time.sleep(PAUSE)

        # No body. The abstract survives in its own column, but nothing is written to
        # fulltext/, which holds full text only. `note` is the signal for a manual fetch.
        if source == "none" and not note:
            note = "no PMC record"

        rows.append({"doi": d, "title": title, "pmcid": pmcid, "source": source,
                     "note": note, "bytes": size, "path": path, "abstract": abstract})

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    import collections
    src = collections.Counter(r["source"] for r in rows)
    mb = sum(int(r["bytes"] or 0) for r in rows) / 1e6
    print(f"\n-> {OUT.name}  {len(rows)} papers, {fetched} fetched this run, {mb:.1f} MB")
    for s in ("manual", "pmc", "biorxiv", "none"):
        if src[s]:
            print(f"   {s:<9}{src[s]:>4}")
    body = sum(src[s] for s in ("manual", "pmc", "biorxiv"))
    print(f"   full body for {body}/{len(rows)}; {src['none']} still need a manual PDF")
    print(f"   titles: {sum(1 for r in rows if r['title'])}/{len(rows)}")


if __name__ == "__main__":
    main()
