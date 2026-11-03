#!/usr/bin/env python3
"""Fetch titles and abstracts for the screened manuscripts from EuropePMC, by DOI.

    python 01_literature/retrieve/fetch_abstracts.py
    python 01_literature/retrieve/fetch_abstracts.py --refresh

Screening 505 papers on title alone is guesswork; the abstract usually settles whether a
study was field-collected and on which host. `studies/` holds full text for only 132 of
them, and nothing in the paper tables carries an abstract.

The title is fetched too, because the one in papers.tsv is frequently not a title: 156 of
505 are whole formatted citations ("T. Adams et al., \u201cRust expression browser,\u201d BMC
Genomics, vol. 22, Mar. 2021, doi: ..."), picked up that way by the reference-manager route.

Results cache to data/bronze/abstracts.json and are reused, so a re-run costs nothing and
only new DOIs are requested. EuropePMC is queried in batches of 25 DOIs, one request at a
time at 2 req/s, which is the limiter EBI applies across EuropePMC and ENA together.
"""

import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
ROOT = MODULE.parent
REPO = ROOT.parent          # research/ and release/ are siblings
PAPERS = REPO / "release" / "01_cohort" / "curation" / "curated_papers.tsv"
CACHE = MODULE / "data" / "bronze" / "abstracts.json"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works/"
MAILTO = "leon.lenzo@curtin.edu.au"
BATCH = 25
PAUSE = 0.6


def clean(text: str) -> str:
    """Abstracts arrive with markup. Keep the section labels, drop the tags.

    EuropePMC returns structured abstracts as `<h4>Background</h4>...`, and the labels are
    worth keeping: "Background" and "Methods" are exactly where a field collection is
    stated. Everything else goes, because this text is rendered into a page.
    """
    t = re.sub(r"<h4[^>]*>(.*?)</h4>", r"\n\1. ", text or "", flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    t = (t.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
          .replace("&quot;", '"').replace("&#x2019;", "'").replace("&nbsp;", " "))
    return re.sub(r"[ \t]+", " ", t).strip()


def from_crossref(doi):
    """Crossref, for what EuropePMC does not hold. Returns (title, abstract).

    Worth the second request: of the 37 papers EuropePMC left without an abstract, the two
    that gate a cohort accession both resolve here. 10.1111/ppa.13468 gains a 1,886-character
    abstract naming wheat, and 10.1080/17429145.2024.2393803 gains a title naming maize.
    Abstracts arrive as JATS, so the tags go the same way as EuropePMC's.
    """
    req = urllib.request.Request(CROSSREF + urllib.parse.quote(doi),
                                 headers={"User-Agent": f"crypt/1.0 (mailto:{MAILTO})"})
    try:
        with urllib.request.urlopen(req, timeout=40) as fh:
            m = json.load(fh)["message"]
    except Exception:
        return "", ""
    title = clean((m.get("title") or [""])[0]).rstrip(".")
    return title, clean(m.get("abstract") or "")


def fetch(dois):
    q = " OR ".join(f'DOI:"{d}"' for d in dois)
    url = API + "?" + urllib.parse.urlencode(
        {"query": q, "resultType": "core", "format": "json", "pageSize": str(BATCH)})
    with urllib.request.urlopen(url, timeout=60) as fh:
        return json.load(fh).get("resultList", {}).get("result", [])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="re-request DOIs already cached")
    args = ap.parse_args()

    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    dois = []
    for r in csv.DictReader(open(PAPERS), delimiter="\t"):
        d = (r.get("doi") or "").strip()
        if d and (args.refresh or d not in cache):
            dois.append(d)
    print(f"{len(cache)} cached, {len(dois)} to request")

    for i in range(0, len(dois), BATCH):
        chunk = dois[i:i + BATCH]
        try:
            res = fetch(chunk)
        except Exception as e:
            print(f"  batch {i // BATCH + 1}: {type(e).__name__} {e}", file=sys.stderr)
            continue
        got = {}
        for r in res:
            doi = (r.get("doi") or "").strip().lower()
            got[doi] = {"title": clean(r.get("title") or "").rstrip("."),
                        "abstract": clean(r.get("abstractText") or ""),
                        "journal": (r.get("journalTitle") or "").strip(),
                        "year": (r.get("pubYear") or "").strip(),
                        "oa": r.get("isOpenAccess") == "Y"}
        # Record the misses too, or every re-run asks for them again forever.
        for d in chunk:
            cache[d] = got.get(d.lower(),
                               {"title": "", "abstract": "", "journal": "", "year": "",
                                "oa": False})
        print(f"  {min(i + BATCH, len(dois)):>4}/{len(dois)}  "
              f"{sum(1 for d in chunk if cache[d]['abstract'])} of {len(chunk)} with an abstract")
        time.sleep(PAUSE)

    # Second pass: anything still missing a title or an abstract goes to Crossref.
    gaps = [d for d, v in cache.items() if not v.get("abstract") or not v.get("title")]
    if gaps:
        print(f"{len(gaps)} without a title or abstract; asking Crossref")
        fixed_t = fixed_a = 0
        for n, d in enumerate(gaps, 1):
            title, abstract = from_crossref(d)
            if title and not cache[d].get("title"):
                cache[d]["title"] = title
                fixed_t += 1
            if abstract and not cache[d].get("abstract"):
                cache[d]["abstract"] = abstract
                fixed_a += 1
            if n % 25 == 0:
                print(f"  {n}/{len(gaps)}")
            time.sleep(PAUSE)
        print(f"  recovered {fixed_t} titles and {fixed_a} abstracts")

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    have = sum(1 for v in cache.values() if v["abstract"])
    print(f"-> {CACHE.relative_to(ROOT)}: {len(cache)} DOIs, {have} with an abstract")


if __name__ == "__main__":
    main()
