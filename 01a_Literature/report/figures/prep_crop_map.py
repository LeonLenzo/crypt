#!/usr/bin/env python3
"""prep_crop_map.py — locality x crop table for the field-cereal map.

Written 2026-10-08. The 02_literature map coloured country circles by collection period, a
single ORDERED value per country. Crop is categorical and several crops can share a locality,
so the unit here is one row per (locality, crop) and the figure plots points, not country
circles.

Crop comes from `gold/cohort_hosts.tsv`, which is written by `report/host_summary.py --per-run`
and is the only place the host is resolved correctly. Do NOT read `run_species.tsv` here: its
`scientific_name` is the PATHOGEN for every rust and FHB survey, which files infected wheat
leaves under Puccinia and reported 56 wheat runs instead of 1,012.

Coordinates come from `gold/cereal_localities.tsv`, written by
`02_literature/03_classify/geocode_localities.py --input gold/cohort_hosts.tsv`, which owns the
string-convention handling and the wrong-type / wrong-country rejections.

Its `precision` column is carried through and matters. The largest sites all name a research
facility no geocoder knows, so they are placed by dropping that name: "USA: Nebraska, Lincoln
(UNL Havelock Research Farm)" resolves as Lincoln, Nebraska. A point marked `region` is only
as good as a province centroid, so the figure reports how many samples sit at each precision
rather than implying every dot is a field.

MIN_RUNS drops crops too small to carry a co-infection rate. Leon's call on 2026-10-08: a
handful of samples has no power here, and plotting them implies a coverage the dataset does
not have. Barley is the one that hurts - 3 samples, both incidental pickups inside rust
surveys - and dropping it is the honest move precisely BECAUSE it is a real gap rather than a
thin measurement. What is dropped is printed and written to crop_totals.tsv with an
`included` flag, so the figure never silently loses a crop.

    python 01a_Literature/report/figures/prep_crop_map.py
"""

from __future__ import annotations

import collections, csv, re, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[2]), str(_HERE.parents[3])]

from _layout import GOLD

HOSTS = GOLD / "cohort_hosts.tsv"
LOCS  = GOLD / "cereal_localities.tsv"
OUT   = GOLD / "crop_map.tsv"
OUTC  = GOLD / "crop_totals.tsv"
OUTY  = GOLD / "crop_year.tsv"

MIN_RUNS = 100          # below this a crop cannot support a co-infection rate

# Order is plot order and legend order: by runs, except that Barley is held out of the
# grouped tail so its scarcity is visible.
CROP = [
    (r"^Triticum",                      "Wheat"),
    (r"^Zea mays$|^Zea mays subsp\.? mays|^Zea mays \(", "Maize"),
    (r"^Oryza sativa",                  "Rice"),
    (r"^Sorghum bicolor",               "Sorghum"),
    (r"^Hordeum",                       "Barley"),
    (r"^Secale cereale|^Triticale|^x?Triticosecale|^Avena|^Setaria|^Panicum|^Eleusine|^Pennisetum",
                                        "Other cereal"),
]
# The archive named these runs after the fungus and no host is recorded anywhere; keeping
# them visible stops the map from understating cereal coverage.
PATHOGEN = re.compile(r"Puccinia|Blumeria|Zymoseptoria|Fusarium|Pyrenophora|Rhizoctonia|^PDA$", re.I)
ORDER = ["Wheat", "Maize", "Rice", "Sorghum", "Barley", "Other cereal", "Cereal, host unresolved"]


def crop_of(host: str) -> str | None:
    h = (host or "").strip()
    for pat, label in CROP:
        if re.match(pat, h, re.I):
            return label
    if PATHOGEN.search(h):
        return "Cereal, host unresolved"
    return None                      # non-cereal: not plotted, but counted in the subtitle


# The same country arrives under several strings, because each came from a different
# submitter. Left unnormalised these draw two diamonds in the same place: "UK" and "United
# Kingdom" both appeared in the country-centroid set. "Zimbabwae" is a typo in the submitted
# BioSample metadata, mapped rather than corrected upstream because the archive string is
# what the curation records.
COUNTRY_CANON = {
    "uk": "United Kingdom", "u.k.": "United Kingdom", "great britain": "United Kingdom",
    "england": "United Kingdom", "scotland": "United Kingdom", "wales": "United Kingdom",
    "usa": "USA", "united states": "USA", "united states of america": "USA",
    "the netherlands": "Netherlands", "czech republic": "Czechia",
    "zimbabwae": "Zimbabwe", "turkiye": "Turkey",
}


def country_of(loc: str) -> str:
    """The country named by a location string, or "" if it names something finer.

    Only the bare-country case is wanted here. "USA: Nebraska, Lincoln" is a locality that
    failed to geocode for some other reason and must not be quietly demoted to a centroid;
    "United Kingdom" is a country and nothing more, and a centroid is the honest placement.
    """
    t = (loc or "").strip().rstrip(".")
    if not t or ":" in t or "," in t:
        return ""
    return COUNTRY_CANON.get(t.lower(), t)


def main() -> None:
    for p in (HOSTS, LOCS):
        if not p.exists():
            sys.exit(f"missing {p}; run host_summary.py --per-run and geocode_localities.py")

    runs = list(csv.DictReader(HOSTS.open(), delimiter="\t"))
    coords = {r["location"]: r for r in csv.DictReader(LOCS.open(), delimiter="\t")}

    cell = collections.defaultdict(lambda: dict(n=0, projects=set(), years=set()))
    ccell = collections.defaultdict(lambda: dict(n=0, projects=set(), years=set()))
    nocountry = 0
    tot = collections.Counter(); tot_loc = collections.defaultdict(set)
    tot_proj = collections.defaultdict(set)
    placed = unplaced = noloc = 0
    noloc_by_crop = collections.Counter()

    for r in runs:
        crop = crop_of(r["host"])
        if crop is None:
            continue
        tot[crop] += 1
        tot_proj[crop].add(r["BioProject"])
        loc = (r["location"] or "").strip()
        # A bare country is not a locality. Counting it as one inflated the localities panel
        # by roughly ten per crop and made "distinct localities" mean two different things
        # in the same figure.
        if loc and not country_of(loc):
            tot_loc[crop].add(loc)
        if not loc:
            noloc += 1; noloc_by_crop[crop] += 1; continue
        g = coords.get(loc)
        # A geocoded row whose precision is "country" was placed at a country centroid even
        # though the string named something finer. Route it to the centroid branch so it is
        # drawn as a diamond, not as a locality circle sitting where nobody sampled.
        if g is not None and g.get("precision") == "country":
            g = None
        if not g:
            # No geocoded point, but the string still names a country. Placing these at the
            # country centroid (leon, 2026-10-09) rather than dropping them: 311 samples,
            # almost all rust surveys where the submitter recorded the country an isolate
            # came from and no collection site. The centroid is a fiction about WHERE and
            # the figure must say so, which is why they carry precision 'country' and are
            # drawn with a different shape rather than silently joining the real points.
            unplaced += 1
            c = country_of(loc) or (loc.split(":", 1)[0].strip() if ":" in loc else "")
            c = COUNTRY_CANON.get(c.lower(), c) if c else ""
            if c:
                k = (c, crop)
                cc = ccell[k]
                cc["n"] += 1
                cc["projects"].add(r["BioProject"])
                if r["collection_date"][:4].isdigit():
                    cc["years"].add(r["collection_date"][:4])
                cc["country"] = c
            else:
                nocountry += 1
            continue
        placed += 1
        c = cell[(loc, crop)]
        c["n"] += 1
        c["projects"].add(r["BioProject"])
        if r["collection_date"][:4].isdigit():
            c["years"].add(r["collection_date"][:4])
        c["lat"], c["lon"], c["country"] = g["lat"], g["lon"], g["country"]
        c["precision"] = g.get("precision", "")

    keep = {c for c in ORDER if tot[c] >= MIN_RUNS}
    dropped = {c: tot[c] for c in ORDER if 0 < tot[c] < MIN_RUNS}

    cols = ["location", "country", "crop", "n", "lat", "lon", "precision", "projects",
            "year_min", "year_max"]
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for (loc, crop), c in sorted(cell.items(), key=lambda kv: -kv[1]["n"]):
            if crop not in keep:
                continue
            yrs = sorted(c["years"])
            w.writerow([loc, c["country"], crop, c["n"], c["lat"], c["lon"],
                        c.get("precision", ""), len(c["projects"]),
                        yrs[0] if yrs else "", yrs[-1] if yrs else ""])
        # Country-centroid rows carry no coordinates; the figure joins them from the map
        # data so there is one definition of a country's centre.
        for (ctry, crop), c in sorted(ccell.items(), key=lambda kv: -kv[1]["n"]):
            if crop not in keep:
                continue
            yrs = sorted(c["years"])
            w.writerow([ctry, ctry, crop, c["n"], "", "", "country", len(c["projects"]),
                        yrs[0] if yrs else "", yrs[-1] if yrs else ""])

    with OUTC.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["crop", "runs", "projects", "localities", "runs_country_only",
                    "runs_unplaced", "included"])
        country_n = collections.Counter()
        for (ctry, crop), c in ccell.items():
            country_n[crop] += c["n"]
        for crop in ORDER:
            if tot[crop]:
                w.writerow([crop, tot[crop], len(tot_proj[crop]), len(tot_loc[crop]),
                            country_n[crop], noloc_by_crop[crop],
                            "yes" if crop in keep else "no"])

    # Year x crop, as a COMPLETE grid. A stacked area chart draws a straight line across a
    # missing year, so a crop absent in 2018 would appear to taper through it rather than
    # stop. Zero-filling every year in the range makes the gaps honest. The range is
    # trimmed to where the data actually are: 6 samples predate 2010 and would otherwise
    # stretch the axis back to 1981 for a line of zeros.
    years = sorted({int(r["collection_date"][:4]) for r in runs
                    if crop_of(r["host"]) in keep and r["collection_date"][:4].isdigit()})
    lo = min(y for y in years if sum(1 for r in runs
                                     if r["collection_date"][:4] == str(y)) >= 10)
    span = [y for y in range(lo, max(years) + 1)]
    ycount = collections.Counter()
    undated = collections.Counter()
    for r in runs:
        c = crop_of(r["host"])
        if c not in keep:
            continue
        d = r["collection_date"][:4]
        if d.isdigit() and int(d) in span:
            ycount[(int(d), c)] += 1
        else:
            undated[c] += 1
    with OUTY.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["year", "crop", "n"])
        for y in span:
            for c in [x for x in ORDER if x in keep]:
                w.writerow([y, c, ycount[(y, c)]])
    print(f"wrote {OUTY}  ({span[0]}-{span[-1]}, "
          f"{sum(ycount.values()):,} dated, {sum(undated.values()):,} outside the range or undated)")

    print(f"wrote {OUT}  ({sum(1 for (l, c) in cell if c in keep)} locality x crop cells)")
    print(f"wrote {OUTC}")
    print(f"  cereal runs: {sum(tot.values()):,} in {len({r['BioProject'] for r in runs if crop_of(r['host'])})} projects")
    print(f"  plotted {placed:,} at a geocoded locality   "
          f"{sum(c['n'] for c in ccell.values()):,} at a country centroid   "
          f"no locality string {noloc:,}   neither {nocountry:,}")
    prec = collections.Counter()
    for c in cell.values():
        prec[c.get("precision", "?")] += c["n"]
    print("  placement precision: " + "  ".join(f"{k}:{v:,}" for k, v in prec.most_common()))
    for crop in ORDER:
        if tot[crop]:
            flag = "" if crop in keep else f"   DROPPED (< {MIN_RUNS} runs)"
            print(f"    {crop:<26}{tot[crop]:>6} runs  {len(tot_loc[crop]):>4} localities{flag}")
    if dropped:
        print(f"  dropped {sum(dropped.values())} runs across {len(dropped)} crops: "
              + ", ".join(f"{k} ({v})" for k, v in dropped.items()))


if __name__ == "__main__":
    main()
