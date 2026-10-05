#!/usr/bin/env python3
"""cohort_provenance.py — where and when the analysed samples came from.

Three sources of location, in descending priority:

  supplementary table    author-curated, per sample, usually the richest: locality, collection
                         date and the author's own field/lab call together. Built by
                         02_literature/02_text/supp_provenance.py. One supplement (Adams et al.
                         2021) covers 1,017 BioSamples across 12 studies.

  geo_loc_name            BioSample attribute, per sample, the only source of WITHIN-project
                          spatial structure (PRJEB39201 has 237 distinct localities)
  llm_geographic_location extracted from the manuscript, per BIOPROJECT, so it fills in
                          projects whose BioSamples carry no attribute at all, but gives one
                          copied value to every sample in a multi-site study

A supplement wins over the BioSample attribute, which wins over the LLM. The attribute beats
the LLM because it is sample-level where the LLM value is a project-level judgment copied down
(`geo_agreement` disagrees on 280 such rows). The supplement beats the attribute because it
carries locality where the attribute often carries only a country, and because it is the only
source for the author's field/lab distinction.

Collection year comes only from `collection_date`. `submission_date` and `pub_date` both
postdate collection, so they are recorded as an upper bound and never as a value.

    python 02_literature/03_classify/cohort_provenance.py
"""
import argparse
import collections
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _paths import CLASSIFIED, ROOT, SAMPLES  # noqa: E402

OUT = Path(__file__).resolve().parent / "data/cohort_provenance.tsv"
SUPP = ROOT / "02_literature/02_text/data/supp_provenance.tsv"
NULL = {"", "missing", "not applicable", "na", "not collected", "unknown", "none", "unclear",
        "not specified", "not reported", "not_stated", "not stated", "not_applicable",
        "not_collected", "restricted access", "not provided"}

# Country strings arrive from three sources with three spelling conventions, plus a few
# sub-national values that leaked into a country field. Without this, UK and United Kingdom plot
# as two places and 61 "countries" are really 54.
COUNTRY_ALIAS = {
    "uk": "United Kingdom", "united kingdom": "United Kingdom",
    "usa": "USA", "united states": "USA", "united states of america": "USA",
    # US states found in the country column
    "california": "USA", "north carolina": "USA", "nc": "USA",
    "the netherlands": "Netherlands",
    "korea": "South Korea",
    "czech republic": "Czechia",
    "zimbabwae": "Zimbabwe",          # misspelling in the source
}
# Not countries. Kept out of country-level figures rather than silently mapped to one.
NOT_A_COUNTRY = {"south america"}


def clean(v):
    v = (v or "").strip()
    return "" if v.lower() in NULL else v


def year(v):
    m = re.search(r"(19|20)\d{2}", v or "")
    return int(m.group(0)) if m else None


def country(v):
    """Normalised country from a location string, or "" if it is not a country.

    NCBI's convention is "Country: region, locality", so the colon split is authoritative when
    present; LLM and supplement values are free text like "Nanjing, China" where the country
    trails. The result is passed through COUNTRY_ALIAS so the three sources agree.
    """
    if not v:
        return ""
    raw = (v.split(":")[0] if ":" in v else v.split(",")[-1]).strip()
    if raw.lower() in NOT_A_COUNTRY:
        return ""
    return COUNTRY_ALIAS.get(raw.lower(), raw)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list-unresolved", action="store_true",
                    help="print the BioProjects with no location, with DOIs, for manual review")
    args = ap.parse_args()

    supp = {}
    if SUPP.exists():
        for r in csv.DictReader(open(SUPP), delimiter="\t"):
            # keep the richest row per BioSample: prefer one carrying a locality
            cur = supp.get(r["BioSample"])
            if not cur or (clean(r.get("location")) and not clean(cur.get("location"))):
                supp[r["BioSample"]] = r

    samples = {r["BioSample"]: r for r in csv.DictReader(open(SAMPLES), delimiter="\t")}
    analysed = {r["biosample"] for r in csv.DictReader(open(CLASSIFIED), delimiter="\t")
                if r["biosample"]}
    sub = [samples[b] for b in analysed if b in samples]

    rows = []
    for r in sub:
        sp = supp.get(r["BioSample"], {})
        sp_loc = clean(sp.get("location")) or clean(sp.get("country"))
        bs_geo = clean(r.get("geo_loc_name"))
        llm_geo = clean(r.get("llm_geographic_location"))
        loc = sp_loc or bs_geo or llm_geo
        src = ("supplement" if sp_loc else "biosample" if bs_geo else
               "llm" if llm_geo else "none")
        # normalise the supplement's country too, not just the parsed one: supplements use
        # their own spellings ("UK", "Czechia") and bypassing country() left them unaliased
        ctry = country(clean(sp.get("country"))) or country(loc)
        # Record WHERE the year came from, not just the year. `location_source` below
        # describes the LOCATION only, and the year is resolved independently, so a consumer
        # that reuses location_source for the year mislabels it. Measured 2026-10-05 while
        # importing this table into 01a_Literature: 52 of 1,843 rows disagreed, including 9
        # where the location is LLM-derived but the year is a sound BioSample value.
        sp_year = year(clean(sp.get("date")))
        bs_year = year(clean(r.get("collection_date")))
        cy = sp_year or bs_year
        ysrc = "supplement" if sp_year else ("biosample" if bs_year else "none")
        # year_upper_bound is a DEPOSIT-time proxy (submission or publication), set only when
        # no collection year exists. It is not an observation and must never be read as one.
        ub = year(r.get("submission_date")) or year(r.get("pub_date"))
        rows.append({
            "BioSample": r["BioSample"], "BioProject": r["BioProject"],
            "location": loc, "country": ctry, "location_source": src,
            "geo_agreement": r.get("geo_agreement", ""),
            "collection_year": cy or "", "year_source": ysrc,
            "year_upper_bound": "" if cy else (ub or ""),
            "year_upper_bound_kind": "" if cy else ("submission_or_pub_date" if ub else ""),
            "author_sample_type": clean(sp.get("sample_type")),
            "doi": r.get("doi", ""), "pmid": r.get("pmid", ""),
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    with open(OUT, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n = len(rows)
    src = collections.Counter(r["location_source"] for r in rows)
    print(f"{OUT}: {n} analysed BioSamples\n")
    print("LOCATION")
    for k in ("supplement", "biosample", "llm", "none"):
        print(f"  {k:10s} {src[k]:5d}  {src[k]/n*100:5.1f}%")
    resolved = n - src["none"]
    print(f"  resolved   {resolved:5d}  {resolved/n*100:5.1f}%")
    dis = sum(1 for r in rows if r["geo_agreement"] == "disagree")
    print(f"  (BioSample preferred over LLM on {dis} disagreeing rows)")
    at = collections.Counter(r["author_sample_type"] for r in rows if r["author_sample_type"])
    if at:
        print(f"  author-declared sample type, where a supplement gives one: {dict(at)}")
        lab = [r for r in rows if r["author_sample_type"] == "Lab"]
        print(f"    -> {len(lab)} samples the authors call LAB, currently in a field cohort")

    print(f"\n  countries ({len({r['country'] for r in rows if r['country']})} distinct):")
    for c, v in collections.Counter(r["country"] for r in rows if r["country"]).most_common(12):
        print(f"    {c:26s} {v:5d}  {v/n*100:5.1f}%")

    cy = [r["collection_year"] for r in rows if r["collection_year"]]
    print(f"\nCOLLECTION YEAR\n  known      {len(cy):5d}  {len(cy)/n*100:5.1f}%")
    ubonly = sum(1 for r in rows if not r["collection_year"] and r["year_upper_bound"])
    print(f"  bound only {ubonly:5d}  {ubonly/n*100:5.1f}%  (submission/publication year, an"
          f" upper bound on collection)")
    print(f"  neither    {n-len(cy)-ubonly:5d}")
    if cy:
        dec = collections.Counter((y // 10) * 10 for y in cy)
        print("  by decade: " + ", ".join(f"{d}s={v}" for d, v in sorted(dec.items())))

    both = [r for r in rows if r["country"] and r["collection_year"]]
    print(f"\n  with BOTH country and collection year: {len(both)} ({len(both)/n*100:.1f}%)"
          f"  <- the set any spatiotemporal claim rests on")

    if args.list_unresolved:
        un = [r for r in rows if r["location_source"] == "none"]
        byp = collections.defaultdict(list)
        for r in un:
            byp[r["BioProject"]].append(r)
        print(f"\n{'='*78}\nUNRESOLVED LOCATION: {len(un)} BioSamples in {len(byp)} BioProjects")
        print(f"{'='*78}")
        print(f"{'BioProject':14s} {'n':>4s}  {'DOI':42s} PMID")
        for p, rs in sorted(byp.items(), key=lambda kv: -len(kv[1])):
            d = next((x["doi"] for x in rs if x["doi"]), "") or "(no DOI)"
            pm = next((x["pmid"] for x in rs if x["pmid"]), "") or "-"
            print(f"{p:14s} {len(rs):4d}  {d:42s} {pm}")


if __name__ == "__main__":
    main()
