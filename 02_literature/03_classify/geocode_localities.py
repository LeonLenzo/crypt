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
import collections
import csv
import json
import re
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
# Canonical name per country, so both sides of the comparison can be normalised. Before
# 2026-10-08 this was a one-way {ours -> {theirs}} map populated only with the names the
# 02_literature pipeline had already normalised. Deriving the country from a raw location
# string instead (01a_Literature holds "Osgodby, UK", not "Osgodby" + country="United
# Kingdom") immediately produced three false rejections reading "resolved to United Kingdom,
# expected UK", because "uk" is not a substring of "united kingdom" in either direction.
_ALIAS = {
    "usa": "united states", "us": "united states", "u.s.a.": "united states",
    "united states of america": "united states",
    "uk": "united kingdom", "u.k.": "united kingdom", "great britain": "united kingdom",
    "england": "united kingdom", "scotland": "united kingdom", "wales": "united kingdom",
    "the netherlands": "netherlands", "holland": "netherlands",
    "czech republic": "czechia", "turkiye": "turkey", "türkiye": "turkey",
    "south korea": "korea", "republic of korea": "korea",
    "russia": "russian federation", "ivory coast": "cote d'ivoire",
}


def canon(name):
    n = re.sub(r"\s+", " ", (name or "").strip().lower())
    return _ALIAS.get(n, n)


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


def query_ladder(location, country):
    """Progressively less specific queries for one locality, most specific first.

    Added 2026-10-08. A single query placed only 304 of 457 localities, and the misses were
    the LARGEST sites, because they name a research facility Nominatim has never heard of:
    "USA: Nebraska, Lincoln (UNL Havelock Research Farm)" (1,572 samples), "China: Shanghai,
    Baihe Experimental Station" (633), "USA: New York, Cornell Musgrave Research Farm" (740).
    Nine thousand seven hundred samples sat unplaced for want of dropping a farm name.

    The ladder is: the full string; the string with any parenthetical removed; then the string
    with its most specific comma-part dropped, repeatedly. Each rung is strictly less precise
    than the one above, so the rung that answers is recorded as `precision` and a figure can
    say which points are a town and which are only a province.
    """
    import re as _re
    loc = location.strip()
    rungs = [loc]
    bare = _re.sub(r"\s*\([^)]*\)", "", loc).strip().rstrip(",")
    if bare and bare != loc:
        rungs.append(bare)
    # Peel the most specific comma-part off the tail of the NCBI form, or the head of the
    # free-text form, since the conventions put specificity at opposite ends.
    cur = bare or loc
    for _ in range(3):
        if ":" in cur:
            head, _, tail = cur.partition(":")
            parts = [x.strip() for x in tail.split(",") if x.strip()]
            if len(parts) <= 1:
                cur = head.strip()
            else:
                cur = head.strip() + ": " + ", ".join(parts[:-1])
        elif "," in cur:
            parts = [x.strip() for x in cur.split(",") if x.strip()]
            cur = ", ".join(parts[1:])
        else:
            break
        if cur and cur not in rungs:
            rungs.append(cur)
    names = ["locality", "locality-nofacility", "region", "subcountry", "country"]
    return [(query_for(r, country), names[min(i, len(names) - 1)])
            for i, r in enumerate(rungs)]


def country_ok(expected, got):
    """Does Nominatim's country agree with the one we already hold?"""
    e, g = canon(expected), canon(got)
    if not e or not g:
        return True
    return e == g or e in g or g in e


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
    ap.add_argument("--input", default=str(PROV),
                    help="TSV with a `location` column; `country` is derived from the NCBI "
                         "`Country: region, locality` prefix when the column is absent. Added "
                         "2026-10-08 so 01a_Literature's cohort can reuse this cache, which "
                         "already holds 460 answers and all of the rejection logic.")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.input), delimiter="\t"))
    want = {}
    for r in rows:
        loc = (r.get("location") or "").strip()
        ctry = (r.get("country") or "").strip()
        if not ctry and ":" in loc:
            ctry = loc.split(":", 1)[0].strip()
        if not ctry and "," in loc:
            ctry = loc.rsplit(",", 1)[1].strip()
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

    ladders = {(l, c): query_ladder(l, c) for l, c in want}
    q_country = {}
    for (l, c), rungs in ladders.items():
        for q, _ in rungs:
            q_country.setdefault(q, c)

    # Walk each locality's ladder and stop at the first rung that answers, so the extra
    # requests are only spent on localities a single query could not place.
    done = 0
    for (l, c), rungs in sorted(ladders.items(), key=lambda kv: -want[kv[0]]):
        for q, _name in rungs:
            if q in cache:
                if cache[q]:
                    break
                continue
            if args.limit and done >= args.limit:
                break
            geocode(q, cache, q_country.get(q, ""))
            done += 1
            if done % 25 == 0:
                print(f"  {done} requests", flush=True)
                CACHE.parent.mkdir(parents=True, exist_ok=True)
                json.dump(cache, open(CACHE, "w"))
            if cache[q]:
                break
    json.dump(cache, open(CACHE, "w"))

    out, hit, miss = [], 0, 0
    prec_n = collections.Counter()
    for (loc, ctry), n in sorted(want.items(), key=lambda kv: -kv[1]):
        g = prec = used = None
        for q, name in ladders[(loc, ctry)]:
            if cache.get(q):
                g, prec, used = cache[q], name, q
                break
        if g:
            hit += n
            prec_n[prec] += n
            out.append({"location": loc, "country": ctry, "n": n, "precision": prec,
                        "query": used, "lat": g["lat"], "lon": g["lon"],
                        "matched": g["display"][:90]})
        else:
            miss += n
    cols = ["location", "country", "n", "lat", "lon", "precision", "query", "matched"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\n{args.out}: {len(out)} localities placed, covering {hit} samples "
          f"({miss} samples' localities unresolved)")
    for k in ("locality", "locality-nofacility", "region", "subcountry", "country"):
        if prec_n[k]:
            print(f"    {k:<22}{prec_n[k]:>7} samples")


if __name__ == "__main__":
    main()
