#!/usr/bin/env python3
"""Blur that matters: per-genus blur against how often the genus is actually CALLED.

A genus can blur badly and still cost nothing if it never appears in the results. This joins
the blur index to the detection table and keeps only genera that are both called and
multi-species, so the top-right corner is the real worry-list: species calls that are common
AND unreliable. A high-blur genus that is never called sits at the far left and can be ignored.

    python3 db1_blur_impact.py <db_inspect.tsv> <detections.tsv> [-o out.png]
"""

import argparse
import csv
import importlib.util
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt       # noqa: E402

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("db1_blur", HERE / "db1_blur.py")
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)

INK, MUTED = "#1a1a1a", "#707070"
CLEAN, WORRY = "#2171b5", "#e05a26"
BLUR_CUT = 0.20                       # above this, species calls lean on genus-shared sequence


def wheat_genus_calls(detections):
    """Genus -> number of wheat, species-rank, positive calls. Fungi only."""
    n = Counter()
    with open(detections) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("call") != "yes":
                continue
            if r.get("kingdom") == "Viridiplantae":
                continue
            host = r.get("host", "").lower()
            if "triticum" not in host and "wheat" not in host:
                continue
            sp = r.get("species", "").split()
            if len(sp) < 2:              # need a binomial to have a genus-vs-species question
                continue
            n[sp[0]] += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inspect")
    ap.add_argument("detections")
    ap.add_argument("-o", "--out", default=str(HERE / "figures/db1_blur_impact.png"))
    args = ap.parse_args()

    genera = bm.blur_by_genus(bm.parse(args.inspect))
    by_name = {g["name"]: g for g in genera.values()}
    calls = wheat_genus_calls(args.detections)

    rows = []
    for name, c in calls.items():
        g = by_name.get(name)
        if not g or g["n_species"] < 2:        # single-species genera have no blur question
            continue
        rows.append((name, c, g["blur"], g["n_species"]))

    if not rows:
        raise SystemExit("no multi-species genera with wheat calls found")

    xs = [c for _, c, _, _ in rows]
    ys = [b * 100 for _, _, b, _ in rows]
    size = [max(14, n * 7) for _, _, _, n in rows]
    colour = [WORRY if b > BLUR_CUT else CLEAN for _, _, b, _ in rows]

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.scatter(xs, ys, s=size, c=colour, alpha=0.5, edgecolors="none", zorder=2)

    # label anything that is either frequently called or badly blurred
    xmax = max(xs)
    for name, c, b, n in rows:
        if c >= 0.04 * xmax or b > 0.30:
            ax.annotate(name, (c, b * 100), xytext=(6, 3), textcoords="offset points",
                        fontsize=8.5, color=INK, style="italic")

    ax.axhline(BLUR_CUT * 100, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(xmax, BLUR_CUT * 100 + 0.8, f"{BLUR_CUT:.0%} blur",
            ha="right", va="bottom", fontsize=7.5, color=MUTED)

    ax.set_xscale("log")
    ax.set_xlabel("wheat species-level calls for the genus  (log)")
    ax.set_ylabel("blur  (% of distinctive sequence shared below genus)")
    ax.set_ylim(-2, max(ys) * 1.1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_title("point size = reference species in genus;  "
                 "top-right = common AND unreliable",
                 fontsize=8.5, color=MUTED, loc="left")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")

    worry = sorted((r for r in rows if r[2] > BLUR_CUT), key=lambda r: -r[1])
    print(f"wrote {args.out}  ({len(rows)} multi-species genera called on wheat)\n")
    print(f"the worry-list (blur > {BLUR_CUT:.0%}), by wheat calls:\n")
    print(f"  {'genus':<20}{'calls':>7}{'blur':>8}{'spp':>5}")
    for name, c, b, n in worry:
        print(f"  {name:<20}{c:>7}{b:>8.1%}{n:>5}")


if __name__ == "__main__":
    import matplotlib.ticker       # noqa: F401
    main()
