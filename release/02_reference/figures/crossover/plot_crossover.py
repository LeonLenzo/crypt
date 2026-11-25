#!/usr/bin/env python3
"""plot_crossover.py - the minimizer-crossover supplementary figure (and its legend).

Rows/cols ordered and bracketed by the NCBI taxon tree (phylum -> ... -> species) as a
cladogram top and left; a category strip per axis using the matrix-figure palette; cells =
crossover (shared / smaller taxon's minimizers). Each of the 25 wheat pathogens is labelled on
the right, and the ones sharing a large fraction of their minimizers with another taxon (hard
to resolve) are highlighted. Heatmap and legend are written as SEPARATE images: a colorbar
attached to the heatmap axis shrinks it out of register with the tree and strips.

    python release/02_reference/figures/crossover/plot_crossover.py

Reads  crossover_matrix.tsv   (beside this script) the 254x254 shared-minimizer matrix,
       computed by compute_matrix.py on Setonix from the db1p CDS (provenance in that file)
       ../../data/outputs/supp_data_s3_reference.tsv   for the per-taxon matrix category
Writes crossover.png, crossover_legend.png  (beside this script)
"""
import argparse, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
S3 = Path(__file__).resolve().parents[2] / "data/outputs/supp_data_s3_reference.tsv"
INK, MUTED = "#1a1a1a", "#707070"
CAT_COL = {"none": "#1aa64b", "plant_pathogen": "#f5a800", "wheat_pathogen": "#e05a26"}
CAT_LAB = {"none": "recorded on wheat (resident)",
           "plant_pathogen": "plant pathogen, not wheat",
           "wheat_pathogen": "wheat pathogen"}
HILITE = "#c0202a"                 # problematic wheat pathogens
PROBLEM_THRESH = 0.20              # max crossover to another taxon to be "problematic"
CMAP = LinearSegmentedColormap.from_list(
    "crossblue", ["#ffffff", "#d0e1f2", "#6baed6", "#2171b5", "#08306b"])
CMAP.set_bad("#e9e9e9")
NORM = PowerNorm(0.35, vmin=0, vmax=1)
CBLABEL = "crossover\n(shared / smaller set)"


def load_categories():
    cat = {}
    with open(S3) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cat[r["ncbi_taxid"]] = r["phibase_col"]
    return cat


def tree_nodes(rank_rows):
    n = len(rank_rows)
    nodes = []
    for level in range(6):
        i = 0
        while i < n:
            key = tuple(rank_rows[i][:level + 1])
            j = i
            while j < n and tuple(rank_rows[j][:level + 1]) == key:
                j += 1
            nodes.append((level, i, j))
            i = j
    return nodes


def draw_tree(ax, nodes, n, orient):
    # leaf i centred at i+0.5 to match the heatmap/strip cell centres
    center = {nd: (nd[1] + nd[2]) / 2.0 for nd in nodes}
    by_level = {}
    for nd in nodes:
        by_level.setdefault(nd[0], []).append(nd)
    segs = []
    for level in range(6):
        for (lv, i0, i1) in by_level[level]:
            kids = ([(l2, a, b) for (l2, a, b) in by_level[level + 1] if a >= i0 and b <= i1]
                    if level < 5 else [])
            cpar = center[(lv, i0, i1)]
            kc = [center[k] for k in kids] if kids else [cpar]
            segs.append(((level, min(kc)), (level, max(kc))))
            for c in kc:
                segs.append(((level, c), (level + 1, c)))
            if not kids:
                segs.append(((level, cpar), (6, cpar)))
    pc = [center[p] for p in by_level[0]]
    segs.append(((0, min(pc)), (0, max(pc))))
    for (d0, p0), (d1, p1) in segs:
        if orient == "left":
            ax.plot([d0, d1], [p0, p1], color=MUTED, lw=0.5)
        else:
            ax.plot([p0, p1], [6 - d0, 6 - d1], color=MUTED, lw=0.5)


def deoverlap(ys, gap, lo, hi):
    """Push a sorted list of desired y positions apart by >=gap, keep within [lo,hi]."""
    ys = sorted(ys)
    out = list(ys)
    for i in range(1, len(out)):
        if out[i] < out[i - 1] + gap:
            out[i] = out[i - 1] + gap
    # if overflow past hi, shift the whole stack up
    if out and out[-1] > hi:
        shift = out[-1] - hi
        out = [y - shift for y in out]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix", nargs="?", default=str(HERE / "crossover_matrix.tsv"))
    ap.add_argument("-o", "--out", default=str(HERE / "crossover.png"))
    ap.add_argument("--legend-out", default=None)
    args = ap.parse_args()
    out = Path(args.out)
    legout = Path(args.legend_out) if args.legend_out else out.with_name(out.stem + "_legend.png")

    rows = list(csv.reader(open(args.matrix), delimiter="\t"))
    data = rows[1:]
    n = len(data)
    taxid = [r[0] for r in data]
    rankrows = [r[1:7] for r in data]
    name = [r[6] for r in data]
    dist = np.array([int(r[7]) for r in data], float)
    M = np.array([[int(x) for x in r[8:8 + n]] for r in data], float)
    cats = [load_categories().get(t, "none") for t in taxid]

    mn = np.minimum.outer(dist, dist)
    with np.errstate(invalid="ignore", divide="ignore"):
        A = np.where(mn > 0, M / mn, 0.0)
    Aoff = A.copy(); np.fill_diagonal(Aoff, 0.0)
    max_cross = Aoff.max(axis=1)                 # per taxon, max crossover to any other
    np.fill_diagonal(A, np.nan)
    nodes = tree_nodes(rankrows)

    wp = [i for i in range(n) if cats[i] == "wheat_pathogen"]
    problem = {i for i in wp if max_cross[i] >= PROBLEM_THRESH}

    fig = plt.figure(figsize=(12, 9.4))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.0, 0.09, 9], height_ratios=[1.0, 0.09, 9],
                          wspace=0.008, hspace=0.008, left=0.03, right=0.72, top=0.99, bottom=0.03)
    ax_top = fig.add_subplot(gs[0, 2]); ax_left = fig.add_subplot(gs[2, 0])
    ax_hcs = fig.add_subplot(gs[1, 2]); ax_vcs = fig.add_subplot(gs[2, 1])
    ax = fig.add_subplot(gs[2, 2])

    ax.imshow(A, cmap=CMAP, norm=NORM, interpolation="nearest", extent=[0, n, n, 0], aspect="auto")
    ax.set_xlim(0, n); ax.set_ylim(n, 0)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(MUTED)

    for i, c in enumerate(cats):
        ax_hcs.add_patch(plt.Rectangle((i, 0), 1, 1, color=CAT_COL[c], lw=0))
        ax_vcs.add_patch(plt.Rectangle((0, i), 1, 1, color=CAT_COL[c], lw=0))
    ax_hcs.set_xlim(0, n); ax_hcs.set_ylim(0, 1); ax_hcs.axis("off")
    ax_vcs.set_xlim(0, 1); ax_vcs.set_ylim(n, 0); ax_vcs.axis("off")

    draw_tree(ax_left, nodes, n, "left")
    ax_left.set_xlim(0, 6); ax_left.set_ylim(n, 0); ax_left.axis("off")
    draw_tree(ax_top, nodes, n, "top")
    ax_top.set_xlim(0, n); ax_top.set_ylim(0, 6); ax_top.axis("off")

    # highlight problematic WP: ring their diagonal cell
    for i in problem:
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fill=False, edgecolor=HILITE, lw=1.4))

    # right-side labels: every wheat pathogen, de-overlapped, with leader lines
    row_y = [i + 0.5 for i in wp]
    lab_y = deoverlap(row_y, gap=n / 42.0, lo=0, hi=n)
    xlab = n + n * 0.015
    for i, yr, yl in zip(wp, row_y, lab_y):
        hot = i in problem
        ax.plot([n, xlab], [yr, yl], color=HILITE if hot else MUTED, lw=0.5, clip_on=False)
        ax.text(xlab + n * 0.004, yl, name[i], va="center", ha="left", fontsize=6.6,
                style="italic", color=HILITE if hot else INK,
                fontweight="bold" if hot else "normal", clip_on=False)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"wrote {out}  ({len(wp)} wheat pathogens, {len(problem)} highlighted)")
    for i in sorted(problem, key=lambda i: -max_cross[i]):
        j = int(Aoff[i].argmax())
        print(f"  problematic: {name[i]:<34} max crossover {max_cross[i]:.2f} with {name[j]}")

    # ---- separate legend ----
    lf = plt.figure(figsize=(3.9, 4.6))
    cax = lf.add_axes([0.12, 0.46, 0.12, 0.48])
    cb = lf.colorbar(ScalarMappable(norm=NORM, cmap=CMAP), cax=cax)
    cb.set_label(CBLABEL, fontsize=9, color=INK); cb.ax.tick_params(labelsize=8, colors=MUTED)
    handles = [Patch(facecolor=CAT_COL[k], label=CAT_LAB[k]) for k in
               ("none", "plant_pathogen", "wheat_pathogen")]
    handles.append(Patch(facecolor="none", edgecolor=HILITE, lw=1.4,
                         label=f"hard to resolve (crossover ≥ {PROBLEM_THRESH:.0%})"))
    lf.legend(handles=handles, loc="lower left", fontsize=9, frameon=False,
              bbox_to_anchor=(0.05, 0.04), title="matrix type (USDA x PHI-base)",
              title_fontsize=9)
    lf.savefig(legout, dpi=170, bbox_inches="tight")
    print(f"wrote {legout}")


if __name__ == "__main__":
    main()
