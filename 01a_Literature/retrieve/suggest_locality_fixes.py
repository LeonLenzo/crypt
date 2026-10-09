#!/usr/bin/env python3
"""suggest_locality_fixes.py — find the real place behind a misspelled locality.

Written 2026-10-09. Leon checked the localities that fell back to a country centroid and found
a typo in every one he looked up: Ada21 writes "Scadden" for Scaddan, and the pattern runs
through the list - Mendelsham for Mendlesham, Edicott for Endicott, Selingenstadt for
Seligenstadt, Medigorria for Mendigorria. Nominatim is an exact matcher, so one transposed
letter sends a real field site to the middle of its country.

## Why not just use a fuzzy geocoder

Tried first, and it is the wrong instrument. Photon indexes the same OSM data with a
typo-tolerant search and does find some - Mendlesham from "Mendelsham", Mendigorria from
"Medigorria", Harvard Forest - but it will not bridge a single INSERTED letter: asked for
"Edicott, Washington" it returns a bookshop in Tuscany and a road in England, while
"Endicott, Washington" returns the village immediately. It is also not country-biased, so a
US query happily answers with Italy. An unguarded fuzzy match is worse than no match.

## What this does instead

Inverts the problem. Rather than asking a gazetteer to guess the spelling, it enumerates every
name within ONE edit of the locality - substitution, deletion, insertion - compiles them into a
single regex, and asks Overpass for place nodes matching it inside the administrative area we
already know. "Edicott" in Washington returns exactly one candidate, Endicott, and nothing
else. The search space is the real gazetteer for that region, so a hit is a place that exists
where the sample says it was.

Single-edit covers what these actually are. A locality needing two edits is reported as having
no candidate rather than being reached for.

## This script PROPOSES

Nothing is applied. Output is a review table; an accepted row is pasted into
geocode_localities.LOCALITY_FIX with its resolved coordinates in the comment. A gazetteer
guess nobody read is how a plausible wrong point gets planted, which is the failure this
module keeps rediscovering.

    python 01a_Literature/retrieve/suggest_locality_fixes.py \\
        --input 01a_Literature/data/gold/cereal_localities.tsv > fixes.tsv
"""

from __future__ import annotations

import argparse, csv, difflib, json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

OVERPASS = "https://overpass-api.de/api/interpreter"
UA = "crypt-research/1.0 (plant pathogen SRA mining; leon.lenzo@curtin.edu.au)"
PLACE_RX = "^(city|town|village|hamlet|locality|suburb|isolated_dwelling|municipality)$"

# OSM's own name for a country, where ours differs.
# US BioSamples write the state as a postal abbreviation; OSM boundaries carry the full name,
# so an area of "MA" matches nothing and the locality falls through as unfindable.
US_STATE = {"AL":"Alabama","AK":"Alaska","AZ":"Arizona","AR":"Arkansas","CA":"California",
 "CO":"Colorado","CT":"Connecticut","DE":"Delaware","FL":"Florida","GA":"Georgia","HI":"Hawaii",
 "ID":"Idaho","IL":"Illinois","IN":"Indiana","IA":"Iowa","KS":"Kansas","KY":"Kentucky",
 "LA":"Louisiana","ME":"Maine","MD":"Maryland","MA":"Massachusetts","MI":"Michigan",
 "MN":"Minnesota","MS":"Mississippi","MO":"Missouri","MT":"Montana","NE":"Nebraska",
 "NV":"Nevada","NH":"New Hampshire","NJ":"New Jersey","NM":"New Mexico","NY":"New York",
 "NC":"North Carolina","ND":"North Dakota","OH":"Ohio","OK":"Oklahoma","OR":"Oregon",
 "PA":"Pennsylvania","RI":"Rhode Island","SC":"South Carolina","SD":"South Dakota",
 "TN":"Tennessee","TX":"Texas","UT":"Utah","VT":"Vermont","VA":"Virginia","WA":"Washington",
 "WV":"West Virginia","WI":"Wisconsin","WY":"Wyoming"}

# Trailing tokens that name a facility rather than the place. "Sinana RS" is a research
# station AT Sinana, and searching the whole string finds nothing.
FACILITY_TAIL = re.compile(r"\s+(RS|ARC|AS|RC|Research\s+(Station|Centre|Center)|Station|"
                           r"Experimental\s+Station|Farm|Field\s+Station|Kebele|Wereda|"
                           r"Woreda|Zone)\s*$", re.I)

AREA = {"usa": "United States", "uk": "United Kingdom", "the netherlands": "Nederland",
        "czechia": "Česko", "turkey": "Türkiye", "germany": "Deutschland",
        "italy": "Italia", "spain": "España", "france": "France", "serbia": "Србија",
        "denmark": "Danmark", "ethiopia": "Ethiopia", "china": "中国"}


def edit1_regex(tok: str) -> str:
    """Every name within one edit of tok, as one anchored regex."""
    ch = list(tok)
    alts = {re.escape(tok)}
    for i in range(len(ch)):
        alts.add(re.escape("".join(ch[:i])) + "." + re.escape("".join(ch[i + 1:])))
        alts.add(re.escape("".join(ch[:i])) + re.escape("".join(ch[i + 1:])))
    for i in range(len(ch) + 1):
        alts.add(re.escape("".join(ch[:i])) + "." + re.escape("".join(ch[i:])))
    return "^(" + "|".join(sorted(alts)) + ")$"


def edit1_match(a_: str, b_: str) -> bool:
    """Is b_ within one edit of a_? Used to assign a batched area result back to a locality."""
    a_, b_ = a_.lower(), b_.lower()
    if a_ == b_:
        return True
    if abs(len(a_) - len(b_)) > 1:
        return False
    if len(a_) == len(b_):
        return sum(x != y for x, y in zip(a_, b_)) == 1
    s_, l_ = (a_, b_) if len(a_) < len(b_) else (b_, a_)
    for i in range(len(l_)):
        if l_[:i] + l_[i + 1:] == s_:
            return True
    return False


def ratio(a_: str, b_: str) -> float:
    return difflib.SequenceMatcher(None, a_.lower(), b_.lower()).ratio()


def overpass(name_rx: str, area: str, tries: int = 3):
    q = (f'[out:json][timeout:45];area["name"="{area}"]["boundary"="administrative"]->.a;'
         f'node["place"~"{PLACE_RX}"]["name"~"{name_rx}",i](area.a);out center 12;')
    for k in range(tries):
        try:
            req = urllib.request.Request(OVERPASS, data=q.encode(), headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as fh:
                return json.load(fh).get("elements", [])
        except Exception as e:
            # Overpass answers 429 and 504 under load; backing off is normal, not an error.
            if k == tries - 1:
                raise
            time.sleep(12 * (k + 1))
    return []


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--precision", default="country")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--sleep", type=float, default=6.0)
    a = ap.parse_args()

    rows = [r for r in csv.DictReader(open(a.input), delimiter="\t")
            if r["precision"] == a.precision]
    rows.sort(key=lambda r: -int(r["n"]))
    if a.limit:
        rows = rows[:a.limit]

    # One Overpass query per AREA, not per locality. Each locality contributes several
    # candidate spellings (as written, facility suffix stripped, each word alone) and each of
    # those expands to an edit-1 regex; querying them separately is up to four requests per
    # locality, 236 for this cohort, against a public API. Unioned per area it is 17.
    jobs = []
    for r in rows:
        loc, ctry = r["location"], r["country"]
        tail = loc.split(":", 1)[1].strip() if ":" in loc else loc
        parts = [p.strip() for p in tail.split(",") if p.strip()]
        name = parts[0]
        area = parts[-1] if len(parts) > 1 else AREA.get(ctry.lower(), ctry)
        area = US_STATE.get(area.upper(), area)
        if len(area) <= 3:
            area = AREA.get(ctry.lower(), ctry)
        cands = [name, FACILITY_TAIL.sub("", name).strip()] + (
            name.split() if len(name.split()) > 1 else [])
        cands = [c for c in dict.fromkeys(cands) if c and len(c) >= 4]
        jobs.append(dict(row=r, name=name, area=area, cands=cands))

    by_area = {}
    for j in jobs:
        by_area.setdefault(j["area"], []).append(j)
    print(f"{len(rows)} localities, {sum(int(r['n']) for r in rows)} samples, "
          f"{len(by_area)} areas", file=sys.stderr)

    found = {}
    for n, (area, js) in enumerate(sorted(by_area.items(), key=lambda kv: -len(kv[1])), 1):
        rx = "^(" + "|".join(sorted({edit1_regex(c)[2:-2] for j in js for c in j["cands"]})) + ")$"
        try:
            els = overpass(rx, area)
        except Exception as e:
            print(f"  {area}: {type(e).__name__}", file=sys.stderr)
            els = []
        found[area] = [(e.get("tags", {}).get("name", ""), e.get("tags", {}).get("place", ""),
                        e["lat"], e["lon"]) for e in els]
        print(f"  [{n}/{len(by_area)}] {area}: {len(found[area])} place(s) within one edit",
              file=sys.stderr)
        time.sleep(a.sleep)

    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["location", "country", "n", "queried", "area", "candidate", "place",
                "lat", "lon", "n_candidates"])
    for j in jobs:
        r, pool = j["row"], found.get(j["area"], [])
        best, used = None, j["name"]
        for cand in j["cands"]:
            hits = [h for h in pool if edit1_match(cand, h[0])]
            if hits:
                best, used = max(hits, key=lambda h: ratio(cand, h[0])), cand
                break
        if best:
            w.writerow([r["location"], r["country"], r["n"], used, j["area"], best[0],
                        best[1], f"{best[2]:.5f}", f"{best[3]:.5f}",
                        sum(1 for h in pool if edit1_match(used, h[0]))])
        else:
            w.writerow([r["location"], r["country"], r["n"], j["name"], j["area"],
                        "", "", "", "", 0])


if __name__ == "__main__":
    main()
