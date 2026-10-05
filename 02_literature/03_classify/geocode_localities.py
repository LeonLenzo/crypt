#!/usr/bin/env python3
"""geocode_localities.py — turn locality strings into coordinates, via OSM Nominatim.

Only ~7% of the cohort carries lat/lon, but 1,677 samples carry a locality NAME across 468
distinct strings. Geocoding those makes a locality-level map possible, which is the only view
that shows what the supplements recovered: PRJEB39201 has 235 localities and PRJEB31334 75,
all invisible on a country-level map.

Three source conventions have to be told apart, because reversing the wrong one puts the point
in the wrong country:

  "USA: California, Davis"   NCBI: Country: region, locality  -> most specific is LAST
  "Nanjing, China"           LLM free text                    -> most specific is FIRST
  "Ickleton"                 supplement, country held apart   -> bare locality

Nominatim's policy is one request per second with a real User-Agent, so a full pass is ~8
minutes. Results are cached on disk, including misses, so reruns cost nothing and adding a
supplement only geocodes what is new.

Two rejections, both found by eyeballing a test batch rather than assumed:

  WRONG TYPE     "De Vries, Netherlands" is a collector surname that leaked into a Location
                 field; Nominatim happily returns a street address in Rotterdam for it. Results
                 whose class is a road, building or address are rejected, so a name that is not
                 a place fails rather than planting a false point.
  WRONG COUNTRY  a bare locality can resolve into the wrong country entirely. The returned
                 address country must match the country we already hold.

    python 02_literature/03_classify/geocode_localities.py
    python 02_literature/03_classify/geocode_localities.py --limit 20   # try a few first
"""
import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROV = HERE / "data/cohort_provenance.tsv"
CACHE = HERE / "data/geocode_cache.json"
OUT = HERE / "data/localities.tsv"
# accept-language=en matters: without it Nominatim returns the country in the local language
# (Slovenija, Deutschland, Србија) and every country check fails
URL = ("https://nominatim.openstreetmap.org/search?q={}&format=json&limit=1"
       "&addressdetails=1&accept-language=en")
# classes that are not settlements: a match here means the string was not a place name
BAD_CLASS = {"highway", "building", "amenity", "shop", "office", "tourism", "man_made"}
# our country name -> what Nominatim calls it, where they differ
COUNTRY_EQ = {"usa": {"united states"}, "united kingdom": {"united kingdom"},
              "czechia": {"czechia", "czech republic"}, "south korea": {"south korea"},
              "netherlands": {"netherlands"}, "turkey": {"turkey", "türkiye"}}
UA = "crypt-research/1.0 (plant pathogen SRA mining; leon.lenzo@curtin.edu.au)"


def query_for(location, country):
    """Build a 'most specific first' query string from whichever convention was used."""
    loc = location.strip()
    if ":" in loc:                                   # NCBI: Country: region, locality
        head, _, tail = loc.partition(":")
        parts = [p.strip() for p in tail.split(",") if p.strip()]
        parts.reverse()                              # locality is last, geocoders want it first
        parts.append(country or head.strip())
    elif "," in loc:                                 # already specific-first free text
        parts = [p.strip() for p in loc.split(",") if p.strip()]
        if country and country.lower() not in {p.lower() for p in parts}:
            parts.append(country)
    else:                                            # bare locality, country held separately
        parts = [loc] + ([country] if country else [])
    seen, out = set(), []
    for p in parts:
        if p.lower() not in seen:
            seen.add(p.lower())
            out.append(p)
    return ", ".join(out)


def country_ok(expected, got):
    """Does Nominatim's country agree with the one we already hold?"""
    e, g = (expected or "").lower(), (got or "").lower()
    if not e or not g:
        return True
    return g in COUNTRY_EQ.get(e, {e}) or e in g or g in e


def geocode(q, cache, expect_country=""):
    if q in cache:
        return cache[q]
    try:
        req = urllib.request.Request(URL.format(urllib.parse.quote(q)),
                                     headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as h:
            js = json.load(h)
        hit = None
        if js:
            r0 = js[0]
            cls = (r0.get("class") or "").lower()
            got_country = (r0.get("address") or {}).get("country", "")
            if cls in BAD_CLASS:
                print(f"    rejected (class={cls}, not a place): {q[:46]}", file=sys.stderr)
            elif not country_ok(expect_country, got_country):
                print(f"    rejected (resolved to {got_country}, expected "
                      f"{expect_country}): {q[:46]}", file=sys.stderr)
            else:
                hit = {"lat": float(r0["lat"]), "lon": float(r0["lon"]),
                       "display": r0.get("display_name", ""), "class": cls}
    except Exception as e:
        print(f"    ! {q[:50]}: {type(e).__name__}", file=sys.stderr)
        hit = None
    cache[q] = hit            # cache misses too, so a rerun does not retry them
    time.sleep(1.1)           # Nominatim policy: max 1 req/s
    return hit


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, help="only geocode this many new strings")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(PROV), delimiter="\t"))
    want = {}
    for r in rows:
        loc, ctry = r["location"].strip(), r["country"].strip()
        if not loc or not ctry:
            continue
        if loc.lower() == ctry.lower():      # country only, nothing to place
            continue
        want.setdefault((loc, ctry), 0)
        want[(loc, ctry)] += 1

    cache = json.load(open(CACHE)) if CACHE.exists() else {}
    queries = {(l, c): query_for(l, c) for l, c in want}
    todo = [q for q in set(queries.values()) if q not in cache]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(want)} distinct localities covering {sum(want.values())} samples")
    print(f"{len(cache)} cached, {len(todo)} to geocode "
          f"(~{len(todo) * 1.1 / 60:.0f} min at Nominatim's 1 req/s)", flush=True)

    q_country = {}
    for (l, c), q in queries.items():
        q_country.setdefault(q, c)
    for i, q in enumerate(todo, 1):
        geocode(q, cache, q_country.get(q, ""))
        if i % 25 == 0:
            print(f"  {i}/{len(todo)}", flush=True)
            CACHE.parent.mkdir(parents=True, exist_ok=True)
            json.dump(cache, open(CACHE, "w"))
    json.dump(cache, open(CACHE, "w"))

    out, hit, miss = [], 0, 0
    for (loc, ctry), n in sorted(want.items(), key=lambda kv: -kv[1]):
        g = cache.get(queries[(loc, ctry)])
        if g:
            hit += n
            out.append({"location": loc, "country": ctry, "n": n,
                        "query": queries[(loc, ctry)], "lat": g["lat"], "lon": g["lon"],
                        "matched": g["display"][:90]})
        else:
            miss += n
    cols = ["location", "country", "n", "lat", "lon", "query", "matched"]
    with open(OUT, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\n{OUT}: {len(out)} localities placed, covering {hit} samples "
          f"({miss} samples' localities unresolved)")


if __name__ == "__main__":
    main()
