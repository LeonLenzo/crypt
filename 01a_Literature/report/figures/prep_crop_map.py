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


def main() -> None:
    for p in (HOSTS, LOCS):
        if not p.exists():
            sys.exit(f"missing {p}; run host_summary.py --per-run and geocode_localities.py")

    runs = list(csv.DictReader(HOSTS.open(), delimiter="\t"))
    coords = {r["location"]: r for r in csv.DictReader(LOCS.open(), delimiter="\t")}

    cell = collections.defaultdict(lambda: dict(n=0, projects=set(), years=set()))
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
        if loc:
            tot_loc[crop].add(loc)
        if not loc:
            noloc += 1; noloc_by_crop[crop] += 1; continue
        g = coords.get(loc)
        if not g:
            unplaced += 1; noloc_by_crop[crop] += 0; continue
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

    with OUTC.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["crop", "runs", "projects", "localities", "runs_unplaced", "included"])
        for crop in ORDER:
            if tot[crop]:
                w.writerow([crop, tot[crop], len(tot_proj[crop]), len(tot_loc[crop]),
                            noloc_by_crop[crop], "yes" if crop in keep else "no"])

    print(f"wrote {OUT}  ({sum(1 for (l, c) in cell if c in keep)} locality x crop cells)")
    print(f"wrote {OUTC}")
    print(f"  cereal runs: {sum(tot.values()):,} in {len({r['BioProject'] for r in runs if crop_of(r['host'])})} projects")
    print(f"  plotted {placed:,}   no locality string {noloc:,}   locality not geocoded {unplaced:,}")
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
