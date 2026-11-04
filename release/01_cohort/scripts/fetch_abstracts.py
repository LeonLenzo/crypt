#!/usr/bin/env python3
"""Fetch the title and abstract for every DOI in data/inputs/pg_dois.txt.

    python release/01_cohort/scripts/fetch_abstracts.py
    python release/01_cohort/scripts/fetch_abstracts.py --refresh   # re-request cached DOIs

Reads  data/inputs/pg_dois.txt      one DOI per line, the corpus
Writes data/raw/abstracts.tsv       doi, title, abstract, pmcid, open_access

This is the CHEAP pass, and it comes first on purpose. screen_abstract.py reads this to
decide a paper's study type, and only papers it does not drop go on to the expensive
fetch_fulltext.py step, which for paywalled papers means a hand-download. There is no point
retrieving a book chapter's body to discover from its title that it is a review.

EuropePMC answers 25 DOIs a request; Crossref fills the gaps (it held titles and abstracts
EuropePMC did not). Both at 2 requests a second. Re-running asks only for DOIs not already
cached, unless --refresh.
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
MODULE = HERE.parent
DATA = MODULE / "data"
DOIS = DATA / "inputs" / "pg_dois.txt"
OUT = DATA / "raw" / "abstracts.tsv"
MAILTO = "leon.lenzo@curtin.edu.au"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works/"
UA = {"User-Agent": f"crypt/1.0 (mailto:{MAILTO})"}
BATCH = 25
PAUSE = 0.6
COLS = ["doi", "title", "abstract", "pmcid", "open_access"]


def clean(text):
    t = re.sub(r"<(?:jats:)?title[^>]*>(.*?)</(?:jats:)?title>", r"\1. ", text or "",
               flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    for a, b in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"), ("&quot;", '"'),
                 ("&#x2019;", "'"), ("&nbsp;", " ")):
        t = t.replace(a, b)
    return " ".join(t.split())


def get(url, headers=None):
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers or {}),
                                timeout=60) as fh:
        return json.load(fh)


def from_epmc(dois):
    url = EPMC + "?" + urllib.parse.urlencode(
        {"query": " OR ".join(f'DOI:"{d}"' for d in dois),
         "resultType": "core", "format": "json", "pageSize": str(BATCH)})
    out = {}
    for r in get(url).get("resultList", {}).get("result", []):
        doi = (r.get("doi") or "").strip().lower()
        out[doi] = {"title": clean(r.get("title") or "").rstrip("."),
                    "abstract": clean(r.get("abstractText") or ""),
                    "pmcid": (r.get("pmcid") or "").strip(),
                    "open_access": "1" if r.get("isOpenAccess") == "Y" else ""}
    return out


def from_crossref(doi):
    try:
        m = get(CROSSREF + urllib.parse.quote(doi), UA)["message"]
    except Exception:
        return "", ""
    return clean((m.get("title") or [""])[0]).rstrip("."), clean(m.get("abstract") or "")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="re-request DOIs already cached")
    args = ap.parse_args()

    dois = [l.strip().lower() for l in DOIS.read_text().splitlines() if l.strip()]
    cache = {}
    if OUT.exists():
        with open(OUT, newline="", encoding="utf-8") as fh:
            cache = {r["doi"]: r for r in csv.DictReader(fh, delimiter="\t")}
    todo = [d for d in dois if args.refresh or d not in cache or not cache[d].get("title")]
    print(f"{len(dois)} DOIs, {len(dois) - len(todo)} cached, {len(todo)} to request")

    got = {}
    for i in range(0, len(todo), BATCH):
        try:
            got.update(from_epmc(todo[i:i + BATCH]))
        except Exception as e:
            print(f"  epmc batch {i // BATCH + 1}: {type(e).__name__} {e}")
        print(f"  epmc {min(i + BATCH, len(todo))}/{len(todo)}")
        time.sleep(PAUSE)

    rows = []
    for d in dois:
        r = got.get(d) or cache.get(d) or {"doi": d, "title": "", "abstract": "",
                                           "pmcid": "", "open_access": ""}
        r["doi"] = d
        if not r.get("title") or not r.get("abstract"):
            ct, ca = from_crossref(d)
            r["title"] = r.get("title") or ct
            r["abstract"] = r.get("abstract") or ca
            time.sleep(PAUSE)
        rows.append({k: r.get(k, "") for k in COLS})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"-> {OUT.name}  {len(rows)} papers, "
          f"{sum(1 for r in rows if r['title'])} titles, "
          f"{sum(1 for r in rows if r['abstract'])} abstracts")


if __name__ == "__main__":
    main()
