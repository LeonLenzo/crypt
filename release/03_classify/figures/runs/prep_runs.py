#!/usr/bin/env python3
"""prep_runs.py — the decision for a handful of runs, one point per candidate taxon.

    python release/03_classify/figures/runs/prep_runs.py          # on Setonix

Reads  ../../data/raw/candidates_db2.tsv        every scored row; SCRATCH ONLY, ~200 MB
       ../../data/inputs/panel_runs.tsv         which runs the panels draw
       ../../data/inputs/thresholds.tsv         the operating point, for the drawn floors
       ../../../02_reference/data/outputs/supp_data_s3_reference.tsv   the labels
Writes runs_points.tsv   one row per candidate in those runs, with its colour and verdict
       runs_floors.tsv   the floors and the null curve, so the R script invents neither

MUST RUN ON SETONIX. The candidate table is one row per run and taxon and never enters the
repo; what comes back is the few thousand rows these panels draw.

## What the panel is for

samples.png is 1,347 columns one pixel wide and a reader has no way to know what a column
IS. These are four of those columns opened up, on the axes the call is actually made on:
redundancy available against coverage observed, both normalised by the taxon's own db2
minimizer total so the null is ONE curve for every taxon rather than one per taxon.

The ladder is the control set, so the panels are a test and not an illustration: a spore
library should show one organism and no host, an uninoculated plant should show no fungus,
and the field sample should show what the chapter is about. Written beside the expectation
in panel_runs.tsv, which was fixed before the results were read.

## Colour is the release palette, not a figure-local one

Points take the reference matrix slot, the same hue as 02_reference's matrix and tree and
03_classify's samples panel, so a category is one colour everywhere in the chapter. The
research-tree version of this figure used its own four-way FUNGuild scheme (animal pathogen
/ endophyte / plant pathogen / unlabelled) and that scheme is retired: FUNGuild genus priors
are what labelled chytrids and industrial yeasts as wheat pathogens.
"""

import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parents[1]
RELEASE = MODULE.parent
sys.path.insert(0, str(RELEASE / "02_reference/scripts"))
sys.path.insert(0, str(MODULE / "scripts"))
from _palette import MATRIX, darken                      # noqa: E402
from detect import load_thresholds, taxon_class, declaration, COHORT   # noqa: E402

# samples.png's two non-matrix colours, so the composite reads as one figure
DECLARED = "#7a7a7a"          # mirrors prep_samples.DECLARED_COLOUR
HOST_COLOUR = "#9ecae1"       # _palette.HAS_GENOME, the scope funnel's light blue

CAND = MODULE / "data/raw/candidates_db2.tsv"
PANEL = MODULE / "data/inputs/panel_runs.tsv"
THRESH = MODULE / "data/inputs/thresholds.tsv"
REF = RELEASE / "02_reference/data/outputs/supp_data_s3_reference.tsv"
PHI = {"none": "none", "plant_pathogen": "vagrant", "wheat_pathogen": "wheat_path"}
csv.field_size_limit(10 ** 7)


def read(p, comment="#"):
    with open(p) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith(comment)),
                                   delimiter="\t"))


# ── the five categories the panels distinguish ───────────────────────────────
# Leon's scheme, 2026-11-25. The earlier version keyed on the reference MATRIX slot alone,
# which is the right language for 02_reference's own figures but the wrong one here: on a
# detection-space panel the question is what each point IS, and a declared pathogen and an
# undeclared one are the same slot and the whole point of the chapter. So the pathogen slot
# SPLITS on the declaration and the three resident slots COLLAPSE to two.
#
# The colours are samples.png's, so the composite reads as one figure: declared grey is
# DECLARED_COLOUR, undeclared is the wheat-pathogen orange, the residents are the matrix
# greens. The host is the one category samples.png never draws, and it takes the scope
# funnel's light blue -- a hue from the release palette that no other category here uses.
CATS = [
    ("Wheat Resident",      MATRIX[("wheat_resident", "none")]),
    ("Plant Resident",      MATRIX[("plant_resident", "none")]),
    ("Declared Pathogen",   DECLARED),
    ("Undeclared Pathogen", MATRIX[("wheat_resident", "wheat_path")]),
    ("Host",                HOST_COLOUR),
]
CAT_FILL = dict(CATS)


def category(klass, rr, crow, taxon, genus):
    """-> (category, fill). The point's identity, which is slot AND declaration."""
    if klass == "host":
        return "Host", CAT_FILL["Host"]
    if klass in ("wheat_pathogen", "plant_pathogen"):
        state = declaration(crow, taxon, genus, False)[0] if crow else "unknown"
        name = "Declared Pathogen" if state == "declared" else "Undeclared Pathogen"
        return name, CAT_FILL[name]
    # residents: the USDA row is the distinction that survives, wheat against other plants
    name = ("Wheat Resident" if rr.get("usda_row") == "wheat_resident"
            else "Plant Resident")
    return name, CAT_FILL[name]


def main():
    panel = read(PANEL)
    want = {r["run"]: r for r in panel}
    th = load_thresholds(THRESH, "db2")
    # the cohort row is declaration()'s first argument, exactly as in prep_samples
    crow = {r["run"]: r for r in read(COHORT)}
    ref = {r["taxon"]: r for r in read(REF)}

    pts = []
    with open(CAND) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["run"] not in want:
                continue
            try:
                red, cov = float(r["redundancy"]), float(r["coverage"])
            except (ValueError, TypeError):
                continue
            if red <= 0 or cov <= 0:
                continue
            rr = ref.get(r["taxon"], {})
            klass = taxon_class(rr, False)
            cat, fill = category(klass, rr, crow.get(r["run"], {}), r["taxon"],
                                 r.get("genus", ""))
            pts.append({"run": r["run"], "rung": want[r["run"]]["rung"],
                        "taxon": r["taxon"], "taxon_class": klass,
                        "redundancy": f"{red:.6g}", "coverage": f"{cov:.6g}",
                        "distinct": r["distinct"], "call": r["call"],
                        "reason": r["reason"], "fill": fill, "edge": darken(fill),
                        "category": cat,
                        # label only what the reader needs to check the verdict: the calls,
                        # plus anything that cleared the absolute floor and still failed
                        "label": r["taxon"] if (r["call"] != "rejected"
                                                or r["reason"] != "distinct") else ""})

    pts.sort(key=lambda p: (int(p["rung"]), -float(p["coverage"])))
    with open(HERE / "runs_points.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(pts[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(pts)

    with open(HERE / "runs_floors.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["item", "value"])
        for k in ("distinct_floor", "coverage_floor", "saturation_floor", "min_redundancy"):
            w.writerow([k, th[k]])
        # the null, y = 1-exp(-x), sampled so R draws the same curve the floors are read on
        for i in range(241):
            x = 10 ** (-5 + i * 7 / 240)
            w.writerow([f"null:{x:.6g}", f"{1 - math.exp(-x):.6g}"])

    print(f"-> {HERE}/runs_points.tsv  {len(pts)} points over {len({p['run'] for p in pts})} runs")
    for r in panel:
        sub = [p for p in pts if p["run"] == r["run"]]
        called = [p for p in sub if p["call"] == "detected"]
        print(f"   rung {r['rung']} {r['run']:12s} {len(sub):5d} candidates, "
              f"{len(called)} called: {', '.join(p['taxon'] for p in called) or 'NONE'}")


if __name__ == "__main__":
    main()
