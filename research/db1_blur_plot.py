#!/usr/bin/env python3
"""Scatter: per-genus blur against the number of reference species in the genus.

Shows that blur does not track species count (if it did, points would climb left to right).
Point area is the genus's total minimizers, so a big faint point is a lot of sequence.
A few genera that carry the argument are labelled; the rest stay unlabelled on purpose.

    python3 db1_blur_plot.py <db_inspect.tsv> [-o out.png]
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt       # noqa: E402

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("db1_blur", HERE / "db1_blur.py")
blur_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(blur_mod)

INK, MUTED = "#1a1a1a", "#707070"
PUCC = "#2171b5"     # clean: blue
BLUR = "#e05a26"     # collapsed: orange
LABEL = {"Puccinia", "Zymoseptoria", "Parastagonospora", "Alternaria",
         "Fusarium", "Aspergillus", "Rhynchosporium", "Colletotrichum", "Ceratocystis"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inspect")
    ap.add_argument("-o", "--out", default=str(HERE / "figures/db1_blur.png"))
    ap.add_argument("--min-mins", type=float, default=5e6)
    args = ap.parse_args()

    genera = blur_mod.blur_by_genus(blur_mod.parse(args.inspect))
    rows = [g for g in genera.values()
            if g["total_mins"] >= args.min_mins and g["n_species"] >= 1]

    x = [g["n_species"] for g in rows]
    y = [g["blur"] * 100 for g in rows]
    size = [max(8, g["total_mins"] / 3e5) for g in rows]
    colour = [BLUR if g["blur"] > 0.20 else PUCC for g in rows]

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.scatter(x, y, s=size, c=colour, alpha=0.45, edgecolors="none", zorder=2)

    # label the argument-carrying genera, nudged off their point
    for g in rows:
        if g["name"] in LABEL:
            ax.annotate(g["name"], (g["n_species"], g["blur"] * 100),
                        xytext=(6, 4), textcoords="offset points",
                        fontsize=8.5, color=INK, style="italic")

    ax.axhline(20, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(ax.get_xlim()[1], 20.8, "20% blur", ha="right", va="bottom",
            fontsize=7.5, color=MUTED)

    ax.set_xscale("log")
    ax.set_xlabel("reference species in genus")
    ax.set_ylabel("blur  (% of distinctive sequence shared below genus)")
    ax.set_ylim(-2, max(y) * 1.08)
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    ax.set_title(f"{len(rows)} genera  (>= {args.min_mins/1e6:.0f}M minimizers)",
                 fontsize=9, color=MUTED, loc="left")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}  ({len(rows)} genera)")


if __name__ == "__main__":
    import matplotlib.ticker       # noqa: E402,F401
    main()
