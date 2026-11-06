#!/usr/bin/env python3
"""prep_cohort_map.py — one row per collection locality, for the cohort map.

    python release/01_cohort/scripts/prep_cohort_map.py

Reads  data/outputs/supp_data_s2_cohort.tsv    the run-level cohort, with coordinates
       data/curation/curated_attributes.tsv    how each coordinate was arrived at
Writes figures/cohort_map_points.tsv           one row per distinct point
       figures/cohort_map_summary.tsv          what is on the map and what is missing

BUILT FROM THE SUPPLEMENTARY TABLE, not from any working file, so the figure and the
supplement describe the same cohort. That rule is inherited from the crop map this
replaces and it is the only way to guarantee a reader can reconcile the two.

The crop panels this follows plotted four cereals and needed colour for crop and a bar
panel to stop the map reading as "mostly wheat". This cohort is wheat, so colour is free
and the comparison panels have nothing to compare. One map, one row per place, size by
how many runs came from it.

FIELD WHEAT ONLY. The map answers where the cohort was collected, so lab material and
off-host runs have no place on it even when they carry coordinates: an isolate deposit
records where its pustule came from, not where anything was sampled. 28 runs were being
drawn on that basis before this filter existed.

COUNTRY-ONLY RUNS ARE DRAWN, AS DIAMONDS. 321 field-wheat runs name a country and
nothing finer, and leaving them off a map that otherwise looks complete overstates its
coverage by a fifth. Every one was chased: three papers' supplements were read in full
and they record country for exactly these samples, so the resolution is the literature's,
not ours. They sit at a country centroid with `placement = Country`, the shape the crop
map used for the same situation.

PRECISION IS CARRIED AND IT MATTERS. A point is only as good as the string it was
geocoded from: `locality` is a named place, `region` is a province centroid, `country` is
a country centroid and is barely a point at all. The summary table reports the mix so the
figure can say what it is rather than implying every dot is a field.

WHAT IS STILL NOT ON THE MAP IS PART OF THE RESULT. A handful of runs name no place at
all, or name one no gazetteer and no country list resolves. The count goes in the
summary and belongs in the caption.
"""

import collections
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
COHORT = MODULE / "data" / "outputs" / "supp_data_s2_cohort.tsv"
CURATED = MODULE / "data" / "curation" / "curated_attributes.tsv"
FIGDIR = MODULE / "figures"
POINTS = FIGDIR / "cohort_map_points.tsv"
SUMMARY = FIGDIR / "cohort_map_summary.tsv"
CENTROIDS = FIGDIR / "country_centroids.tsv"

# Coordinates are rounded before grouping. Two runs from one field can differ in the
# sixth decimal, which is under a metre and would otherwise split one locality into two
# points drawn on top of each other. Three places is about 100 m.
PLACES = 3

# Centroid of each country's LARGEST polygon, not the midpoint of its bounding box. A
# bounding box is wrong for any country with distant territory: it puts France in the
# Atlantic and the USA in the Pacific. Taken from the same source the R map draws.
# Spellings differ between what papers write and what the world map calls a country.
# The map's name wins, since that is what the centroid is keyed on: it says "UK", not
# "United Kingdom". A name absent from the map is reported, never guessed.
COUNTRY_FIX = {"United Kingdom": "UK", "Czechia": "Czech Republic",
               "Türkiye": "Turkey", "Zimbabwae": "Zimbabwe",
               "Korea": "South Korea", "Russia": "Russia"}


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def field_wheat(run):
    """The material the map is about. Mirrors the network figure's test deliberately:
    the two panels have to be describing the same samples."""
    return (run.get("assay_in_scope", "1") == "1" and run["setting"] == "field"
            and run["host"].startswith("Triticum"))


def country_of(location):
    """The country a location string names, in the spelling the world map uses."""
    head = location.split(":", 1)[0].strip() if ":" in location else location.strip()
    if "," in head and ":" not in location:
        head = head.rsplit(",", 1)[-1].strip()
    return COUNTRY_FIX.get(head, head)


def main():
    cohort = [r for r in read(COHORT) if field_wheat(r)]
    precision = {r["run"]: r.get("location_precision", "") for r in read(CURATED)}
    centroid = {r["country"]: (float(r["lat"]), float(r["lon"]))
                for r in read(CENTROIDS)}

    placed = [r for r in cohort if r["lat"] and r["lon"]]
    rest = [r for r in cohort if not (r["lat"] and r["lon"])]

    rows = []
    points = collections.defaultdict(list)
    for run in placed:
        points[(round(float(run["lat"]), PLACES),
                round(float(run["lon"]), PLACES))].append(run)
    for (lat, lon), runs in sorted(points.items(), key=lambda kv: -len(kv[1])):
        mix = collections.Counter(precision.get(r["run"], "") or "uncurated"
                                  for r in runs)
        rows.append({"lat": f"{lat:.{PLACES}f}", "lon": f"{lon:.{PLACES}f}",
                     "runs": len(runs), "placement": "Locality",
                     "bioprojects": len({r["BioProject"] for r in runs}),
                     "precision": mix.most_common(1)[0][0],
                     "location": sorted({r["location"] for r in runs
                                         if r["location"]})[0]
                     if any(r["location"] for r in runs) else ""})

    # One diamond per country, carrying every run from it that has no finer place.
    by_country = collections.defaultdict(list)
    homeless = []
    for run in rest:
        name = country_of(run["location"]) if run["location"] else ""
        (by_country[name] if name in centroid else homeless).append(run)
    for name, runs in sorted(by_country.items(), key=lambda kv: -len(kv[1])):
        lat, lon = centroid[name]
        rows.append({"lat": f"{lat:.{PLACES}f}", "lon": f"{lon:.{PLACES}f}",
                     "runs": len(runs), "placement": "Country",
                     "bioprojects": len({r["BioProject"] for r in runs}),
                     "precision": "country", "location": name})

    FIGDIR.mkdir(exist_ok=True)
    with open(POINTS, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    mix = collections.Counter(precision.get(r["run"], "") or "uncurated"
                              for r in placed)
    summary = [
        {"measure": "field-wheat runs in cohort", "value": len(cohort)},
        {"measure": "drawn at a locality", "value": len(placed)},
        {"measure": "drawn at a country centroid",
         "value": sum(len(v) for v in by_country.values())},
        {"measure": "not drawn: no usable place", "value": len(homeless)},
        {"measure": "locality points", "value": sum(1 for r in rows
                                                    if r["placement"] == "Locality")},
        {"measure": "country points", "value": len(by_country)},
    ] + [{"measure": f"locality runs at precision: {k}", "value": v}
         for k, v in mix.most_common()]
    with open(SUMMARY, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["measure", "value"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(summary)

    print(f"-> {POINTS.name}  {len(rows)} points")
    for row in summary:
        print(f"   {row['measure']:<42}{row['value']:>6}")
    if homeless:
        print("   unplaceable location strings: "
              + str(dict(collections.Counter(r["location"] or "(blank)"
                                             for r in homeless))))


if __name__ == "__main__":
    main()
