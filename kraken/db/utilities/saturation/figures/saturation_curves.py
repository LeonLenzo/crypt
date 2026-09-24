#!/usr/bin/env python3
"""
saturation_curves.py — figure for kraken_db_saturation.py.

Three panels answering "how many assemblies per species does the Kraken2 build need":
  A  k-mer accumulation, as fold increase over a single assembly
  B  marginal gain, the % of new k-mers contributed by the nth assembly
  C  whether extra assemblies cost species-level resolution, measured as the
     fraction of a species' k-mers shared with its congeners in the same build

Run from crypt/: python kraken/db/figures/saturation_curves.py
Output: kraken/db/utilities/saturation/figures/saturation_curves.pdf + .png
"""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = Path("kraken/db/utilities/saturation/data")
OUT_DIR  = Path("kraken/db/utilities/saturation/figures")

BLUE, ORANGE = "#2a78d6", "#eb6834"
GREY, FAINT, INK, MUTED = "#8a8985", "#d8d7d2", "#0b0b0b", "#52514e"

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def load():
    curve = list(csv.DictReader(open(DATA_DIR / "saturation_curve.tsv"), delimiter="\t"))
    summ  = list(csv.DictReader(open(DATA_DIR / "saturation_summary.tsv"), delimiter="\t"))
    by = defaultdict(list)
    for r in curve:
        by[r["organism_name"]].append(r)
    for k in by:
        by[k].sort(key=lambda r: int(r["n_assemblies"]))
    return by, summ


def panel_a(ax, by, summ):
    ratio = {r["organism_name"]: float(r["pangenome_ratio"]) for r in summ}
    named = sorted(ratio, key=lambda k: -ratio[k])[:3]

    for name, rows in by.items():
        if name in named:
            continue
        y = np.array([float(r["mean_kmers"]) for r in rows])
        y = y / y[0]
        ax.plot(np.arange(1, len(y) + 1), y, color=GREY, lw=0.7, alpha=0.45,
                solid_capstyle="round")

    # Median over a FIXED cohort. Taking it over "every taxon with at least n
    # assemblies" makes the set shrink with n, which puts steps in the line that
    # are set membership changing, not k-mer content accumulating.
    DEPTH = 8
    cohort = [v for v in by.values() if len(v) >= DEPTH]
    med = [(n, np.median([float(v[n - 1]["mean_kmers"]) / float(v[0]["mean_kmers"])
                          for v in cohort])) for n in range(1, DEPTH + 1)]
    ax.plot([m[0] for m in med], [m[1] for m in med], color=BLUE, lw=2.0,
            solid_capstyle="round", zorder=5,
            label="median, {} species with $\\geq${} assemblies".format(len(cohort), DEPTH))

    for name in named:
        y = np.array([float(r["mean_kmers"]) for r in by[name]])
        y = y / y[0]
        x = np.arange(1, len(y) + 1)
        ax.plot(x, y, color=ORANGE, lw=1.4, solid_capstyle="round", zorder=4)
        sp = name.split()
        ax.annotate("{}. {}".format(sp[0][0], " ".join(sp[1:])),
                    xy=(x[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=7, style="italic", color=ORANGE)

    # Heaps' law: N(n) = k.n^gamma is a straight line on log-log, so the guide line
    # makes the median exponent readable as a slope rather than leaving the panel
    # looking like the curves fail to bend.
    gam = []
    for rows in by.values():
        n = np.arange(1, len(rows) + 1)
        N = np.array([float(r["mean_kmers"]) for r in rows])
        gam.append(np.polyfit(np.log(n), np.log(N), 1)[0])
    gmed = float(np.median(gam))
    gx = np.array([1, 300.0])
    ax.plot(gx, gx ** gmed, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=2)
    ax.annotate(r"median $N \propto n^{%.2f}$, so nothing plateaus" % gmed,
                xy=(300, 300 ** gmed), xytext=(-2, -12), textcoords="offset points",
                ha="right", va="top", fontsize=7, color=MUTED)

    ax.set_yscale("log", base=2)
    ax.set_yticks([1, 2, 4, 8, 16])
    ax.set_yticklabels(["1x", "2x", "4x", "8x", "16x"])
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 3, 5, 10, 20, 50, 150])
    ax.set_xticklabels(["1", "2", "3", "5", "10", "20", "50", "150"])
    ax.set_xlim(0.93, 400)
    ax.set_xlabel("assemblies included")
    ax.set_ylabel("distinct 31-mers, relative to one assembly")
    ax.set_title("A   Pan-genome accumulation, log-log", loc="left", color=INK, pad=8)
    ax.grid(axis="y", color=FAINT, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right", handlelength=1.4,
              borderaxespad=0.2, bbox_to_anchor=(1.0, -0.02))


def panel_b(ax, by):
    ns, meds, q1, q3 = [], [], [], []
    for n in range(2, 16):
        vals = [float(v[n - 1]["mean_frac_new"]) * 100 for v in by.values() if len(v) >= n]
        if len(vals) < 5:
            break
        ns.append(n); meds.append(np.median(vals))
        q1.append(np.percentile(vals, 25)); q3.append(np.percentile(vals, 75))

    ax.axhline(5, color=MUTED, lw=0.7, ls=(0, (3, 3)), zorder=1)
    ax.annotate("5% of current content", xy=(ns[-1], 5), xytext=(-2, 4),
                textcoords="offset points", ha="right", fontsize=7, color=MUTED)
    ax.fill_between(ns, q1, q3, color=BLUE, alpha=0.16, lw=0, zorder=2)
    ax.plot(ns, meds, color=BLUE, lw=2.0, solid_capstyle="round", zorder=3)
    ax.plot(ns, meds, "o", color=BLUE, ms=4.5, mec="white", mew=1.2, zorder=4)
    for n, m in zip(ns, meds):
        if n <= 5:
            ax.annotate("{:.0f}%".format(m), xy=(n, m), xytext=(7, 6),
                        textcoords="offset points", ha="left", fontsize=7.5, color=INK)

    ax.set_xticks(ns)
    ax.set_ylim(0, max(q3) * 1.28)
    ax.set_xlabel("nth assembly added")
    ax.set_ylabel("new 31-mers contributed (%)")
    ax.set_title("B   Marginal gain per assembly", loc="left", color=INK, pad=8)
    ax.annotate("falls as $\\gamma/n$, never to zero", xy=(0.97, 0.87),
                xycoords="axes fraction", ha="right", va="top", fontsize=7, color=MUTED)
    ax.grid(axis="y", color=FAINT, lw=0.6)
    ax.set_axisbelow(True)
    ax.annotate("median, shaded IQR", xy=(0.97, 0.95), xycoords="axes fraction",
                ha="right", va="top", fontsize=7, color=MUTED)


def panel_c(ax, summ):
    pts = [(float(r["frac_shared_1_assembly"]) * 100, float(r["frac_shared_all"]) * 100,
            r["organism_name"]) for r in summ if int(r["n_congener_taxa"]) > 0]
    ax.plot([0, 100], [0, 100], color=MUTED, lw=0.7, ls=(0, (3, 3)), zorder=1)
    worse  = [p for p in pts if p[1] > p[0] + 0.5]
    better = [p for p in pts if p[1] <= p[0] + 0.5]
    ax.scatter([p[0] for p in better], [p[1] for p in better], s=26, color=BLUE,
               edgecolor="white", linewidth=0.8, zorder=3)
    ax.scatter([p[0] for p in worse], [p[1] for p in worse], s=30, color=ORANGE,
               edgecolor="white", linewidth=0.8, zorder=4)
    for x, y, n in worse:
        sp = n.split()
        ax.annotate("{}. {}".format(sp[0][0], " ".join(sp[1:])), xy=(x, y),
                    xytext=(7, -1), textcoords="offset points", fontsize=7,
                    style="italic", color=ORANGE)

    ax.annotate("shared fraction fell\nin {} of {} species".format(len(better), len(pts)),
                xy=(0.96, 0.08), xycoords="axes fraction", ha="right", fontsize=7.5,
                color=MUTED, linespacing=1.5)
    ax.set_xlim(-4, 104); ax.set_ylim(-4, 104)
    ax.set_xticks([0, 25, 50, 75, 100]); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("shared with congeners, 1 assembly (%)")
    ax.set_ylabel("shared with congeners, all (%)")
    ax.set_title("C   Cost to species-level resolution", loc="left", color=INK, pad=8)
    ax.grid(color=FAINT, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_aspect("equal", adjustable="box")


def main():
    by, summ = load()
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    panel_a(axes[0], by, summ)
    panel_b(axes[1], by)
    panel_c(axes[2], summ)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.86, bottom=0.16, wspace=0.36)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        p = OUT_DIR / "saturation_curves.{}".format(ext)
        fig.savefig(p, dpi=300)
        print("wrote", p)


if __name__ == "__main__":
    main()
