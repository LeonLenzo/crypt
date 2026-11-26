#!/usr/bin/env python3
"""Prepare the purpose test: what the 135 predicted references did to the rank profile.

    python release/03_classify/figures/resolution/prep_resolution.py

THIS SCRIPT DECIDES; plot_resolution.R DRAWS. The axis transform is part of the decision
and is applied here, so R receives positions and tick labels and invents neither.

Reads  ../../data/outputs/arm_resolution.tsv   one row per run, both arms, by rank bucket
Writes resolution_curves.tsv   one point per run per bucket: position, cumulative share
       resolution_ticks.tsv    where the axis labels go, in transformed space
       resolution_summary.tsv  per bucket: medians, how many runs moved, total reads

## Two panels, and they are different quantities on purpose

`resolution_reads.png` is the headline in three numbers: the net READS that came to rest at
each rank in db2 against db1p, across the whole cohort. `resolution.png` is the per-run
evidence that the direction is not a handful of deep libraries: the distribution of each
run's change in SHARE. The summary table carries both, and they do not agree to the run,
because a run can gain reads at a rank while that rank's share of its fungal reads falls.
The caption says so rather than leaving a reader to find it.

## What the figure is for

`db2` is `db1p` plus 135 gene-predicted genomes and nothing else, so any difference
between the arms is attributable to those references. The module README states the result
in prose: the predictions take reads OFF species and put them ON the genus node, which is
tiers A and C2 doing what they were built to do and not tier B failing. This draws it.

## An ECDF, because the change is small in most runs and large in a few

The first attempt was a paired scatter, each run a point at (db1p share, db2 share) with
the 1:1 line, and it cannot work on this data: the median run moves 0.19 of a percentage
point while the tail reaches 29, so the axis range is set by a handful of runs and the
systematic shift is about one percent of the panel width. An ECDF spends its axis on
QUANTILES instead of on magnitudes, so the median shift and the tail are both legible, and
it still says one thing per run rather than pooling them.

## The x axis is symmetric-log, and the transform is here rather than in the scale

Changes span 0.001 to 29 percentage points in both directions and include exact zero, so
neither a linear nor a log axis can hold them. The transform is

    t(v) = sign(v) * log10(1 + |v| / LIN)

which is linear within about LIN of zero and logarithmic outside it, with t(0) = 0 exactly.
LIN is quoted in the caption, because a reader has to know where the axis stops being
linear to read the middle of it.

## Fungal reads only, and every classified run

The host is most of a plant RNA-seq library and resolves at species, so it sets the
all-reads denominator and buries a fungal shift below the second decimal place; the
fungal restriction is compare_arms.py's lineage test on the Fungi node. All 1,560 runs
classified in both arms are drawn, not the analysed set the other panels narrow to: the
question here is about the DATABASE, and every run is evidence about it.

## Three buckets, because they sum to one

Reads that leave species have to arrive somewhere. Drawing species and genus alone would
leave a reader to assume the remainder went nowhere.
"""

import csv
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
FIG = Path(__file__).resolve().parent
RES = HERE / "data/outputs/arm_resolution.tsv"

ARMS = ("db1p", "db2")
# Ordered as the rank ladder, specific to broad, which is also how the legend reads.
BUCKETS = (("f_species", "species"), ("f_genus", "genus"),
           ("f_above_genus", "above genus"))
# The scope funnel's blues, broad -> specific, borrowed from _palette for the same reason
# the other figures borrow their hues: one chapter, one set of colours. Here the ladder is
# taxonomic rank rather than database scope, and the ordering is the same kind of thing.
COLOUR = {"f_species": "#6baed6", "f_genus": "#2171b5", "f_above_genus": "#08306b"}
LIN = 0.05              # percentage points; inside this the axis is linear
TICKS = (0.1, 1, 10)    # decades labelled either side of zero, plus zero itself


def t(v):
    """The symmetric-log position of a change of v percentage points."""
    return math.copysign(math.log10(1 + abs(v) / LIN), v)


def main():
    with open(RES) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    curves, summary = [], []
    for key, label in BUCKETS:
        pts = []
        for r in rows:
            d1, d2 = int(r["db1p_f_reads"]), int(r["db2_f_reads"])
            if not d1 or not d2:
                continue            # no fungal read came to rest anywhere in this run
            s1, s2 = int(r[f"db1p_{key}"]) / d1, int(r[f"db2_{key}"]) / d2
            pts.append((r["run"], s1, s2, (s2 - s1) * 100))
        pts.sort(key=lambda p: p[3])
        n = len(pts)
        for i, (run, s1, s2, d) in enumerate(pts):
            curves.append({"bucket": label, "run": run,
                           # the step goes at the TOP of each point, so the curve reads as
                           # "this share of runs changed by at most this much"
                           "y": f"{(i + 1) / n:.6f}", "delta_pp": f"{d:.6g}",
                           # no colour column: the bucket's hex is written ONCE, in
                           # resolution_summary.tsv, and 4,680 copies of it would be a
                           # quarter of this file
                           "x": f"{t(d):.6f}"})
        deltas = [p[3] for p in pts]
        summary.append({
            "bucket": label, "key": key, "colour": COLOUR[key], "n": n,
            "median_db1p": f"{statistics.median(p[1] for p in pts):.6f}",
            "median_db2": f"{statistics.median(p[2] for p in pts):.6f}",
            "median_delta_pp": f"{statistics.median(deltas):.4f}",
            "n_up": sum(1 for d in deltas if d > 0),
            "n_down": sum(1 for d in deltas if d < 0),
            "n_same": sum(1 for d in deltas if d == 0),
            # the same question asked of READS rather than of the share, which is what the
            # bar panel draws: a run can gain reads at a rank and still lose share there
            "n_up_reads": sum(1 for r in rows
                              if int(r[f"db2_{key}"]) > int(r[f"db1p_{key}"])),
            "n_down_reads": sum(1 for r in rows
                                if int(r[f"db2_{key}"]) < int(r[f"db1p_{key}"])),
            "min_pp": f"{min(deltas):.4f}", "max_pp": f"{max(deltas):.4f}",
            "total_delta_reads": sum(int(r[f"db2_{key}"]) - int(r[f"db1p_{key}"])
                                     for r in rows)})

    ticks = [{"value": 0, "x": "0.000000", "label": "0"}]
    for v in TICKS:
        for s in (-1, 1):
            ticks.append({"value": s * v, "x": f"{t(s * v):.6f}",
                          "label": ("−" if s < 0 else "+") + f"{v:g}"})
    ticks.sort(key=lambda r: float(r["x"]))

    def dump(name, rows_):
        with open(FIG / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows_[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows_)

    dump("resolution_curves.tsv", curves)
    dump("resolution_ticks.tsv", ticks)
    dump("resolution_summary.tsv", summary)

    print(f"wrote {FIG}/resolution_{{curves,ticks,summary}}.tsv")
    print(f"  {summary[0]['n']:,} runs classified in both arms, "
          f"{len(rows) - summary[0]['n']} with no fungal read excluded")
    for s in summary:
        print(f"  {s['bucket']:12s} median {float(s['median_db1p']):7.2%} -> "
              f"{float(s['median_db2']):7.2%}  delta {float(s['median_delta_pp']):+6.3f} pp "
              f"(up {s['n_up']:,}, down {s['n_down']:,}, unchanged {s['n_same']:,}), "
              f"total {s['total_delta_reads']:+,d} reads")


if __name__ == "__main__":
    main()
