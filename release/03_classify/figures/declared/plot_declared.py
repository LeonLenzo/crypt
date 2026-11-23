#!/usr/bin/env python3
"""Panel A: what each study declared against what its reads contain.

    python release/03_classify/figures/declared/plot_declared.py

Reads  ../../data/outputs/supp_data_s5_detections.tsv    the calls
       ../../../01_cohort/data/outputs/supp_data_s2_cohort.tsv   runs and declarations
Writes declared.png, declared_legend.png, declared_cells.tsv  (beside this script)

One row per BioProject, one column per wheat pathogen, cell shaded by the share of that
project's runs carrying the call. **The outlined cell is the pathogen the study declared.**

The panel carries two claims at once and that is why it is the lead. The outlined cells are
the positive control: the screen recovers the organism the study went looking for, in nearly
every run of nearly every project, with no tuning per project. Everything else in the row is
the result, and it is read as a SPREAD across projects rather than a pooled rate
(`crypt-prevalence-demoted`): the project is the unit of analysis.

ONLY THE 12 PROJECTS WITH A CURATED DECLARATION ARE DRAWN. For the other 18 (344 runs)
nobody has established what the study reported, so their detections are unknown rather than
undeclared and an outline would be a claim the curation does not support. The caption says
so; the supplementary table carries them.

Zero cells are not drawn. An outlined cell with no fill would be a declared pathogen the
screen failed to find, and the absence of any such cell is itself the control passing.

Every number is read at draw time; none is transcribed.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _palette import WHEAT_PATH, INK, MUTED             # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from _cohort import is_analysed, describe               # noqa: E402

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"
FIG = Path(__file__).resolve().parent
csv.field_size_limit(10 ** 7)

# One hue, white to the reference figure's wheat-pathogen orange, so panel and reference read
# as one scheme. Sequential data gets one hue (leon-figure-style).
# A drawn cell means the call was made, so the palest call must be visibly a call: the ramp
# starts at a warm tint rather than at white, or a 3-in-550 cell renders white and reads as an
# empty box. The mapping is otherwise LINEAR. A square root was tried and rejected: it lifts
# the 1-3% tail to the weight of a 10% cell, and the panel then reads as signal everywhere.
RAMP = LinearSegmentedColormap.from_list(
    "wp", ["#fdeee1", "#fbd5b5", WHEAT_PATH, "#8c3310"])
FLOOR = 0.06
tint = lambda s: RAMP(FLOOR + (1 - FLOOR) * s)
TICKS = (0.0, 0.25, 0.50, 0.75, 1.00)

# A row is drawn if the taxon is called in at least this many of the curated runs. Below it a
# row is one or two cells in 1,208 runs and costs a tenth of the panel's height to say so; the
# caption names what fell out and `declared_cells.tsv` keeps every number. A DECLARED taxon is
# never cut, whatever its count -- the control must be complete.
MIN_RUNS = 10


def abbrev(name):
    """'Puccinia striiformis' -> 'P. striiformis'. Genus-only names pass through."""
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


def load():
    """-> (projects, taxa, share, declared, runs)

    `taxa` is ordered by the number of runs the taxon is called in. **`projects` is ordered by
    its declared target's position in that order**, size breaking the tie, so the outlined
    cells descend as a staircase instead of scattering: the eye reads the control as one shape
    and everything off it as the result. `declared` holds (project, taxon) pairs, matched on the
    SPECIES BINOMIAL because a declaration is curated at forma specialis ('Puccinia
    striiformis f. sp. tritici') and neither database holds a sub-species node, so a call can
    only ever be the species. Matching on the full string would mark every target undeclared.
    """
    runs, declared_raw, status = defaultdict(set), defaultdict(set), {}
    with open(S2) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not is_analysed(r):          # lab libraries and the host contrast; see _cohort
                continue
            runs[r["BioProject"]].add(r["run"])
            status[r["BioProject"]] = r["declared_status"]
            for t in filter(None, (x.strip() for x in r["declared_pathogens"].split(";"))):
                declared_raw[r["BioProject"]].add(" ".join(t.split()[:2]))

    keep = {p for p, s in status.items() if s == "curated"}
    hit = defaultdict(set)
    with open(S5) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r["is_wheat_pathogen"] == "yes" and r["BioProject"] in keep
                    and r["run"] in runs[r["BioProject"]]):
                hit[(r["BioProject"], r["taxon"])].add(r["run"])

    n_called = lambda t: len({r for (p, tt), v in hit.items() if tt == t for r in v})
    declared_any = {t for p in keep for t in declared_raw[p]}
    taxa = {t for _, t in hit} | declared_any
    dropped = {t: n_called(t) for t in taxa
               if t not in declared_any and n_called(t) < MIN_RUNS}
    # NAME BREAKS EVERY TIE. Without it the order falls out of set iteration, and two taxa
    # or two projects on the same count reordered between runs: declared_cells.tsv stopped
    # being reproducible the moment the curation added projects that tied.
    taxa = sorted(taxa - set(dropped), key=lambda t: (-n_called(t), t))
    rank = {t: i for i, t in enumerate(taxa)}
    projects = sorted(keep, key=lambda p: (min((rank[t] for t in declared_raw[p]), default=99),
                                           -len(runs[p]), p))
    share = {k: len(v) / len(runs[k[0]]) for k, v in hit.items()}
    declared = {(p, t) for p in keep for t in declared_raw[p]}
    return (projects, taxa, share, declared, {p: len(runs[p]) for p in projects}, hit,
            # name breaks the tie, or two taxa on one run reorder between invocations and
            # the TSV's row order stops being reproducible
            dict(sorted(dropped.items(), key=lambda kv: (-kv[1], kv[0]))))


def main():
    projects, taxa, share, declared, n_runs, hit, dropped = load()
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.facecolor": "white",
                         "figure.facecolor": "white"})

    # TAXA ARE ROWS. Species names are long and set italic, so as rows they read horizontally
    # at full size; as columns they have to be rotated and they then set the panel's height.
    # The project accessions are short and uniform, so rotating those costs almost nothing.
    CW, CH, GAP = 0.38, 0.29, 0.035         # cell width, height, space between, in inches
    X0, Y0 = 1.30, 0.40                     # grid origin measured from the figure corner
    W = X0 + len(projects) * (CW + GAP) + 0.25
    H = Y0 + len(taxa) * (CH + GAP) + 1.05
    fig, ax = plt.subplots(figsize=(W, H))

    for i, t in enumerate(taxa):
        y = Y0 + (len(taxa) - 1 - i) * (CH + GAP)
        ax.text(X0 - 0.11, y + CH / 2, abbrev(t), ha="right", va="center", fontsize=7.4,
                style="italic", color=INK)
        for j, p in enumerate(projects):
            x = X0 + j * (CW + GAP)
            sh = share.get((p, t), 0.0)
            is_dec = (p, t) in declared
            if not sh and not is_dec:
                continue                    # the gap is the statement; see the caption
            ax.add_patch(Rectangle(
                (x, y), CW, CH, facecolor=tint(sh) if sh else "#ffffff",
                edgecolor=INK if is_dec else "#dddddd",
                linewidth=1.7 if is_dec else 0.6, zorder=3 if is_dec else 2))
            # Only the cells a reader would otherwise have to estimate carry a number.
            # Below 10% the ramp is enough and the digits are chartjunk (leon-figure-style).
            if sh >= 0.10:
                ax.text(x + CW / 2, y + CH / 2, f"{round(sh * 100):d}", ha="center",
                        va="center", fontsize=6.6, zorder=4,
                        color="#ffffff" if sh > 0.60 else INK)

    gh = len(taxa) * (CH + GAP) - GAP
    for j, p in enumerate(projects):
        x = X0 + j * (CW + GAP) + CW / 2
        ax.text(x, Y0 + gh + 0.09, p, rotation=55, rotation_mode="anchor", ha="left",
                va="baseline", fontsize=7.0, color=INK)
        ax.text(x, Y0 - 0.12, f"{n_runs[p]:,}", ha="center", va="top", fontsize=6.6,
                color=MUTED)
    ax.text(X0 - 0.11, Y0 - 0.12, "runs", ha="right", va="top", fontsize=6.6, color=MUTED)

    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "declared.png", dpi=300)

    legend()

    with open(FIG / "declared_cells.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["BioProject", "n_runs", "taxon", "n_runs_with_call", "share", "declared"])
        for p in projects:
            for t in list(taxa) + list(dropped):      # dropped rows stay in the table
                k = (p, t)
                if k in hit or k in declared:
                    w.writerow([p, n_runs[p], t, len(hit.get(k, ())),
                                f"{share.get(k, 0.0):.4f}", int(k in declared)])

    d_sh = [share.get(k, 0.0) for k in declared]
    print(f"wrote {FIG / 'declared.png'}")
    print(f"  {len(projects)} projects, {sum(n_runs.values()):,} {describe()} runs, "
          f"{len(taxa)} taxa drawn")
    print(f"  below MIN_RUNS={MIN_RUNS}, not drawn: " +
          ", ".join(f"{t} ({n})" for t, n in dropped.items()))
    print(f"  {len(declared)} declared cells, recovery min {min(d_sh):.2f} "
          f"median {sorted(d_sh)[len(d_sh) // 2]:.2f}")


def legend():
    """The key is its own image so the panel can be placed without reserving space for it.

    The bar is drawn through `tint`, the same function the cells use, and the ticks are placed
    LINEARLY because the mapping is linear. They were briefly left at `t ** 0.5` from a ramp
    that had been reverted, which put the 25% tick at the bar's midpoint: a legend that
    disagrees with its panel is worse than no legend.
    """
    W, H = 3.00, 0.78
    BX, BY, BW, BH = 0.10, 0.40, 1.95, 0.17      # colour bar, in inches
    fig, ax = plt.subplots(figsize=(W, H))
    for i in range(160):
        ax.add_patch(Rectangle((BX + i * BW / 160, BY), BW / 160 + 0.002, BH,
                               facecolor=tint(i / 159), edgecolor="none"))
    ax.add_patch(Rectangle((BX, BY), BW, BH, facecolor="none", edgecolor="#d8d8d8",
                           linewidth=0.6))
    for t in TICKS:
        ax.text(BX + t * BW, BY - 0.045, f"{round(t * 100)}", ha="center", va="top",
                fontsize=7.0, color=INK)
    ax.text(BX, BY + BH + 0.05, "% of the project's runs with the call", ha="left",
            va="bottom", fontsize=7.4, color=INK)
    ax.add_patch(Rectangle((BX, 0.06), 0.19, 0.14, facecolor="#ffffff", edgecolor=INK,
                           linewidth=1.7))
    ax.text(BX + 0.27, 0.13, "declared target of the study", ha="left", va="center",
            fontsize=7.4, color=INK)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(FIG / "declared_legend.png", dpi=300, facecolor="white")


if __name__ == "__main__":
    main()
