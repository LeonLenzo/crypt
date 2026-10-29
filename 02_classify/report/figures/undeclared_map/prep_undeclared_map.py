#!/usr/bin/env python3
"""prep_undeclared_map.py -- one map, coloured by how many pathogens a library carries that
its own study did not report.

This is the companion to prep_pathogen_map.py and asks a different question. That script draws
one panel per species and answers "where does this species turn up". This one collapses the
species axis and answers "how many unreported pathogens does a sample from here carry".

## WHAT IS COLOURED, AND WHY IT IS NOT THE SPECIES COUNT

Colour is the MEDIAN NUMBER OF UNDECLARED PATHOGENS PER LIBRARY at that locality, not the
number of distinct undeclared species found there. The distinct count is written to the table
as `n_undeclared_species`, but it must not drive colour: a locality with forty libraries finds
more species than one with two for reasons of effort alone, and since effort also sets point
SIZE, the two channels would encode the same confound twice and the eye would read it as
biology. The per-library median has a denominator and is comparable between a site with two
libraries and a site with forty.

`n_undeclared` comes from detect.py via coinfection.tsv rather than being re-derived here, so
"undeclared" means exactly what it means everywhere else in the module: a species-rank plant
pathogen called in a library whose own BioProject did not declare it.

## ALL CROPS, NOT WHEAT ONLY

prep_pathogen_map.py is wheat-only because its six panels are the wheat foliar disease complex.
Restricting here would drop the maize sites, which as of 2026-10-23 are the heaviest co-infection
signal in the cohort, so the crop is carried as a column and drawn as a shape instead.

## THE DENOMINATOR IS STILL DRAWN

Localities where nothing undeclared was found are kept, with zero, and drawn as open points.
An undeclared-pathogen map that shows only the hits is a map of where people sampled.

Coordinate precision follows prep_pathogen_map.py: 'country' and 'region' precision name a
country and nothing finer, and are dropped rather than given a centroid that would invent a
field site at a place nobody sampled.

Writes undeclared_map.tsv and undeclared_map_regions.tsv.

    python 02_classify/report/figures/undeclared_map/prep_undeclared_map.py
"""

import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import COHORT_SAMPLES, GOLD  # noqa: E402
from _util import atomic_write  # noqa: E402

GOOD_PRECISION = {"locality", "locality-nofacility",
                  "corrected-high", "corrected-medium", "corrected-low"}
ALIASES = {"United Kingdom": "UK", "The Netherlands": "Netherlands", "Zimbabwae": "Zimbabwe"}

# Windows are drawn only where there is something to draw; the counts go in the caption so a
# reader can see an empty panel was a choice rather than an oversight.
# WHOLE CONTINENTS, not boxes cropped to the points. A reader needs to see that a continent
# was largely unsampled, and a box drawn around the points hides exactly that: cropping to the
# data makes sparse sampling look like thorough sampling of a small place. Leon's call,
# 2026-10-26. A point near a boundary can fall inside two windows and is drawn in both, which
# is correct for a map and costs nothing.
REGIONS = [
    ("Europe",            (-25.0, 45.0),  (34.0, 71.0)),
    ("North America",    (-168.0, -52.0), (7.0, 72.0)),
    ("South America",     (-82.0, -34.0), (-56.0, 13.0)),
    ("Africa",            (-19.0, 52.0),  (-36.0, 38.0)),
    ("Asia",               (40.0, 150.0), (5.0, 60.0)),
    ("Oceania",           (110.0, 180.0), (-48.0, -9.0)),
]

TSV = HERE / "undeclared_map.tsv"
REG_TSV = HERE / "undeclared_map_regions.tsv"


def main() -> None:
    coinf = {r["run"]: r for r in csv.DictReader(open(GOLD / "coinfection.tsv"), delimiter="\t")}
    prov = {r["Run"]: r for r in csv.DictReader(open(COHORT_SAMPLES), delimiter="\t")}

    runs, dropped = [], defaultdict(int)
    for run, c in coinf.items():
        p = prov.get(run, {})
        prec = (p.get("location_precision") or "").strip()
        lat, lon = (p.get("lat") or "").strip(), (p.get("lon") or "").strip()
        if not (lat and lon):
            dropped["no coordinates"] += 1
            continue
        if prec not in GOOD_PRECISION:
            dropped[f"precision '{prec or 'blank'}'"] += 1
            continue
        raw = (p.get("location") or "").split(":")[0].strip()
        und = [s for s in (c.get("undeclared") or "").split(";") if s]
        runs.append({
            "run": run, "proj": c["BioProject"], "crop": c["crop"],
            "loc": (p.get("location") or "").strip(),
            "country": ALIASES.get(raw, raw), "lat": float(lat), "lon": float(lon),
            "n_und": int(c["n_undeclared"] or 0), "und": und,
            "selection": c.get("sampling_selection", ""),
        })

    # Keyed on locality x CROP, not locality. Each crop is drawn as its own figure, so each
    # needs its own denominator: a site sampled for wheat and for maize is one grey point on
    # the wheat map and a separate one on the maize map, never a "mixed" point on both.
    sites = defaultdict(lambda: {"runs": [], "species": set(), "projs": set()})
    for r in runs:
        s = sites[(r["loc"], r["crop"])]
        s["runs"].append(r)
        s["species"].update(r["und"])
        s["projs"].add(r["proj"])
        s["lat"], s["lon"], s["country"] = r["lat"], r["lon"], r["country"]

    rows = []
    for (loc, crop), s in sorted(sites.items()):
        n = [r["n_und"] for r in s["runs"]]
        rows.append({
            "locality": loc, "country": s["country"],
            "lat": f"{s['lat']:.5f}", "lon": f"{s['lon']:.5f}",
            "crop": crop,
            "n_libraries": len(n), "n_projects": len(s["projs"]),
            "median_undeclared": f"{st.median(n):.1f}",
            "mean_undeclared": f"{sum(n) / len(n):.2f}",
            "max_undeclared": max(n),
            "n_undeclared_species": len(s["species"]),
            "libraries_with_any": sum(1 for x in n if x > 0),
            "top_species": ";".join(sorted(s["species"])[:6]),
        })
    cols = list(rows[0])
    with atomic_write(TSV) as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # One row per region x crop: the R script draws a continent panel only where that crop was
    # actually sampled, rather than printing empty continents for every crop.
    with atomic_write(REG_TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["region", "crop", "xmin", "xmax", "ymin", "ymax", "n_localities",
                    "n_libraries", "median_undeclared", "max_undeclared"])
        crops = sorted({r["crop"] for r in rows})
        for name, (x0, x1), (y0, y1) in REGIONS:
            for crop in crops:
                inside = [r for r in rows if r["crop"] == crop
                          and x0 <= float(r["lon"]) <= x1 and y0 <= float(r["lat"]) <= y1]
                if not inside:
                    continue
                libs = sum(r["n_libraries"] for r in inside)
                med = st.median([float(r["median_undeclared"]) for r in inside])
                w.writerow([name, crop, x0, x1, y0, y1, len(inside), libs,
                            f"{med:.1f}", max(r["max_undeclared"] for r in inside)])
        print("  region x crop coverage:")
        for crop in crops:
            sub = [r for r in rows if r["crop"] == crop]
            regs = [n for n, (x0, x1), (y0, y1) in REGIONS
                    if any(r["crop"] == crop and x0 <= float(r["lon"]) <= x1
                           and y0 <= float(r["lat"]) <= y1 for r in rows)]
            print(f"    {crop:<14}{len(sub):>4} locality-rows {sum(r['n_libraries'] for r in sub):>5} "
                  f"libraries   {', '.join(regs)}")

    tot = sum(r["n_libraries"] for r in rows)
    print(f"\n-> {TSV.name}: {len(rows)} localities, {tot} libraries with usable coordinates")
    print(f"   dropped: " + ", ".join(f"{v} {k}" for k, v in sorted(dropped.items(),
                                                                    key=lambda kv: -kv[1])))


main()
