#!/usr/bin/env python3
"""Fetch the title and abstract for every DOI in inputs/dois.txt.

    python release/01_cohort/fetch_title_abstract.py

Reads  inputs/dois.txt     one DOI per line, the literature search's result
Writes raw/raw_papers.tsv  doi, title, abstract

EuropePMC answers in batches of 25. Crossref fills the gaps, which is worth the second
request: it supplied 26 abstracts and 37 titles EuropePMC did not hold. Both are polite at
2 requests a second.
"""

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOIS = HERE / "inputs" / "dois.txt"
OUT = HERE / "raw" / "raw_papers.tsv"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works/"
MAILTO = "leon.lenzo@curtin.edu.au"
BATCH = 25
PAUSE = 0.6


def clean(text):
    """Strip the markup both sources wrap abstracts in, keeping section labels."""
    t = re.sub(r"<(?:jats:)?title[^>]*>(.*?)</(?:jats:)?title>", r"\1. ", text or "",
               flags=re.S | re.I)
    t = re.sub(r"<h4[^>]*>(.*?)</h4>", r"\1. ", t, flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    for a, b in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"), ("&quot;", '"'),
                 ("&#x2019;", "'"), ("&nbsp;", " ")):
        t = t.replace(a, b)
    return " ".join(t.split())


def get(url, headers=None):
    req = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.load(fh)


def from_europepmc(dois):
    query = " OR ".join(f'DOI:"{d}"' for d in dois)
    url = EPMC + "?" + urllib.parse.urlencode(
        {"query": query, "resultType": "core", "format": "json", "pageSize": str(BATCH)})
    found = {}
    for r in get(url).get("resultList", {}).get("result", []):
        doi = (r.get("doi") or "").strip().lower()
        found[doi] = (clean(r.get("title") or "").rstrip("."),
                      clean(r.get("abstractText") or ""))
    return found


def from_crossref(doi):
    try:
        m = get(CROSSREF + urllib.parse.quote(doi),
                {"User-Agent": f"crypt/1.0 (mailto:{MAILTO})"})["message"]
    except Exception:
        return "", ""
    return clean((m.get("title") or [""])[0]).rstrip("."), clean(m.get("abstract") or "")


def main():
    dois = [l.strip().lower() for l in DOIS.read_text().splitlines() if l.strip()]
    papers = {d: ("", "") for d in dois}

    for i in range(0, len(dois), BATCH):
        chunk = dois[i:i + BATCH]
        try:
            papers.update(from_europepmc(chunk))
        except Exception as e:
            print(f"  europepmc batch {i // BATCH + 1}: {type(e).__name__} {e}")
        print(f"  europepmc {min(i + BATCH, len(dois))}/{len(dois)}")
        time.sleep(PAUSE)

    gaps = [d for d, (t, a) in papers.items() if not t or not a]
    print(f"  {len(gaps)} missing a title or abstract; asking crossref")
    for d in gaps:
        t, a = from_crossref(d)
        papers[d] = (papers[d][0] or t, papers[d][1] or a)
        time.sleep(PAUSE)

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["doi", "title", "abstract"])
        for d in dois:
            t, a = papers[d]
            w.writerow([d, t, a])

    have_t = sum(1 for d in dois if papers[d][0])
    have_a = sum(1 for d in dois if papers[d][1])
    print(f"-> {OUT.name}  {len(dois)} DOIs, {have_t} titles, {have_a} abstracts")


if __name__ == "__main__":
    main()
