#!/usr/bin/env python3
"""Panel A: the reference as a crossed-axis matrix of USDA record x PHI-base score.

    python release/02_reference/figures/matrix/plot_matrix.py

Reads  ../../data/outputs/supp_data_s3_reference.tsv   the published reference table
Writes matrix.png, matrix_numbers.tsv  (beside this script; the figure owns its directory)

A taxon carries two INDEPENDENT facts: where USDA recorded it, and what PHI-base scores it
as. They do not nest, so drawn as a Venn the awkward combination (a plant pathogen USDA
records on wheat, which PHI-base never scores on wheat) needs an overlap lens with no name,
and every box ends up wearing an inclusive set name while holding an exclusive count.
Crossed axes draw as a MATRIX, where that combination is simply a cell.

    rows      USDA Plant-Fungal Database: Wheat Resident / Plant Resident
    columns   PHI-base: Non Pathogen / Plant Pathogen / Wheat Pathogen

THE ROW AXIS IS WHERE THE TAXON IS NAMED, NOT HOW IT ENTERED. USDA either records this
species on wheat or it does not. Whether it reached the reference because USDA named the
species or because USDA named its genus is a fact about the database build, it is in the
`route` column, and drawing it here made the figure assert something about the ORGANISM that
the route cannot support.

THE EMPTY CELLS ARE NOT DRAWN, and their absence is the result. The grid is a staircase
because three of the nine combinations do not occur, for two different reasons that the
caption states it: no wheat pathogen in scope is off-wheat in USDA, two independent sources
agreeing, which is the design test passing.

Colour carries both axes so the figure needs no legend: hue is the column, saturation is the
row. See `_palette.MATRIX`. The unit is a TAXON, not a described species -- ~40% of the
genus-only row is NCBI `sp.` nodes, unidentified isolates that are real Kraken2 competitors
but not species in any taxonomic sense.

THE PLANT ROW HOLDS THE GENUS-ROUTE TAXA and does not separate them, deliberately. USDA
holding no record under a species name is a fact about USDA's coverage of a NAME and never
evidence about the organism, and the 3,465 taxa in that position carry 15.0% of all
field-wheat detections. See README, "The genus tier is not a weak tier".

PRED. IS A GENE MODEL WE MADE, not one NCBI carries, and it is a separate count for that
reason. Every predicted taxon had no annotated assembly at all, so the Pred. number is the
taxa each cell gained that the database could not have been built from otherwise. 101 of
the 135 are in the wheat rows, which is where the worklist aimed them.

Every number is read from the curation tables at draw time; none is transcribed.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parents[2]          # release/02_reference
sys.path.insert(0, str(HERE / "scripts"))           # the shared palette lives with the pipeline
from _palette import MATRIX, darken, INK            # noqa: E402

S3 = HERE / "data/outputs/supp_data_s3_reference.tsv"
FIG = Path(__file__).resolve().parent
csv.field_size_limit(10 ** 7)

ROWS = ["wheat_resident", "plant_resident"]
COLS = ["none", "vagrant", "wheat_path"]
# the column keys in S3 are the full names; the short ones above are this script's own
S3_COL = {"none": "none", "vagrant": "plant_pathogen", "wheat_path": "wheat_pathogen"}
ROW_LABEL = {"wheat_resident": "Wheat Resident", "plant_resident": "Plant Resident"}
COL_LABEL = {"none": "Non\nPathogen", "vagrant": "Plant\nPathogen",
             "wheat_path": "Wheat\nPathogen"}
STATS = [("Taxa", 0), ("Asm.", 1), ("CDS", 2), ("Pred.", 3)]


def grid():
    """-> (cells, totals): cells[(row, col)] = (taxa, assemblies, cds), straight from S3.

    One pass over the published table and nothing else. The panel used to re-derive its
    cells from three curation tables plus the taxdump, which meant the figure and the
    deliverable could disagree; now the deliverable IS the figure's input.
    """
    cells = {(r, c): [0, 0, 0, 0] for r in ROWS for c in COLS}
    inv = {v: k for k, v in S3_COL.items()}
    with open(S3) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            # S3 now publishes the HOST SINK as well as the fungal scope (route=host,
        # kingdom=Viridiplantae, empty fungal columns), because a reference table that
        # omitted the wheat in the database would misdescribe what was built. This panel
        # is the fungal scope, so it filters. Without the filter a host row keys on an
        # empty phibase_col and raises.
            if r["kingdom"] != "Fungi":
                continue
            key = (r["usda_row"], inv[r["phibase_col"]])
            cell = cells[key]
            cell[0] += 1
            cell[1] += int(r["n_assemblies"])
            cell[2] += int(r["n_annotated"])
            cell[3] += int(r["n_predicted"])
    cells = {k: tuple(v) for k, v in cells.items()}
    add = lambda ks: tuple(sum(cells[k][i] for k in ks) for i in range(4))
    totals = {("row", r): add([(r, c) for c in COLS]) for r in ROWS}
    totals.update({("col", c): add([(r, c) for r in ROWS]) for c in COLS})
    totals[("all", "cohort")] = add(list(cells))
    return cells, totals


def main():
    C, T = grid()
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.facecolor": "white",
                         "figure.facecolor": "white"})

    # Sized to be PLACED, not cropped later. Data units are inches, so every constant below
    # is a real measurement on the page.
    #
    # The page is A4 PORTRAIT (8.27 x 11.69 in) with the tree taking the full width on top
    # and this panel beneath it, beside the legend.
    #
    # THE AXIS TITLES ARE GONE, Leon's call: with the row and column labels naming the
    # source in words ("Wheat species", "Wheat Pathogen"), "USDA Plant-Fungal Database" down
    # the side and "PHI-Base" across the top were restating what the caption says anyway,
    # and the rotated title was the only reason the left gutter had to be wide.
    #
    # The panel is 3.85 x 2.17: two rows of four-line pills under the column labels. It is
    # sized to be PLACED, so it states its own size rather than assuming an older layout.
    CW, CH, GAP = 0.87, 0.88, 0.07          # pill width, height, and the space between
    X0, YTOP = 1.10, 1.83                   # grid origin: left edge, top edge
    PAD, ROUND = 0.085, 0.075               # inside a pill, and its corner radius
    LEFT = 0.06                             # margin left of the longest row label
    fig, ax = plt.subplots(figsize=(X0 + 3 * CW + 2 * GAP + 0.06 - LEFT, YTOP + 0.34))

    for i, rn in enumerate(ROWS):
        y = YTOP - (i + 1) * CH - i * GAP
        ax.text(X0 - 0.13, y + CH / 2, ROW_LABEL[rn], ha="right", va="center",
                fontsize=8.5, color=INK, linespacing=1.3)
        for j, cn in enumerate(COLS):
            s_, a_, c_, p_ = C[(rn, cn)]
            if not s_:
                continue                    # the gap IS the statement; see the caption
            x = X0 + j * (CW + GAP)
            fc = MATRIX[(rn, cn)]
            ax.add_patch(FancyBboxPatch(
                (x, y), CW, CH, boxstyle=f"round,pad=0,rounding_size={ROUND}",
                facecolor=fc, edgecolor=darken(fc), linewidth=1.0, zorder=2))
            for k, (lab, idx) in enumerate(STATS):
                ty = y + CH - CH * (k + 0.5) / len(STATS)
                ax.text(x + PAD, ty, lab, ha="left", va="center", fontsize=6.8,
                        style="italic", color=INK, zorder=3)
                ax.text(x + CW - PAD, ty, f"{(s_, a_, c_, p_)[idx]:,}", ha="right",
                        va="center", fontsize=8.0, color=INK, zorder=3)

    for j, cn in enumerate(COLS):
        ax.text(X0 + j * (CW + GAP) + CW / 2, YTOP + 0.07, COL_LABEL[cn], ha="center",
                va="bottom", fontsize=7.8, color=INK, linespacing=1.3)
    gw = 3 * CW + 2 * GAP
    nr = len(ROWS)

    ax.set_xlim(LEFT, X0 + gw + 0.06)
    ax.set_ylim(YTOP - nr * CH - (nr - 1) * GAP - 0.06, YTOP + 0.34)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "matrix.png", dpi=300)

    with open(FIG / "matrix_numbers.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["row", "col", "taxa", "assemblies", "cds", "predicted"])
        for rn in ROWS:
            for cn in COLS:
                w.writerow([rn, cn, *C[(rn, cn)]])
        for (kind, k), v in T.items():
            w.writerow([kind, k, *v])
    drawn = sum(1 for rn in ROWS for cn in COLS if C[(rn, cn)][0])
    tot = sum(C[(rn, cn)][0] for rn in ROWS for cn in COLS)
    print(f"wrote {FIG / 'matrix.png'}")
    print(f"  {drawn} cells drawn, {len(ROWS) * len(COLS) - drawn} empty; {tot:,} taxa, "
          f"cohort is {T[('all', 'cohort')][0]:,}")


if __name__ == "__main__":
    main()
