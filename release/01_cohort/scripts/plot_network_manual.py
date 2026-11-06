#!/usr/bin/env python3
"""Draw the cohort network from hand-placed node positions. No solver.

    python release/01_cohort/scripts/plot_network_manual.py

Reads  figures/cohort_network_nodes.tsv   one row per node: angle_deg and radius
       data/outputs/supp_data_s1_studies.tsv, supp_data_s2_cohort.tsv,
       data/curation/screen_fulltext.tsv, relations.tsv   via plot_cohort_network
Writes figures/cohort_network_manual.png
       figures/cohort_network_manual_legend.png

The companion to plot_cohort_network.py, which solves the layout: that one spends
minutes deciding where 31 BioProjects go so their edges stay radial and clear each
other's samples. This one does none of that. Every paper and BioProject sits exactly
where `cohort_network_nodes.tsv` says, and the only thing computed is the cloud of runs
around each BioProject, because a sample has no position worth arguing about.

Get a starting table from `plot_cohort_network.py --map`, which writes the solver's
answer in this format, then move whatever you want. Nothing here checks your placement:
clouds may overlap and edges may run through them. That is the trade. A solver that
cannot be overruled is worse than one you can ignore.

WHY THE FIGURE IS SIZED FROM ITS CONTENTS. The table is in points, because a point is
what marker sizes and clearances are in, and the only way a radius of 236 points means
236 points on the page is for one data unit to be one point. That fixes the canvas at
the extent of what is drawn rather than at a chosen width, which is why there is no
FIGSIZE here.
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_cohort_network as net

NODES = net.FIGDIR / "cohort_network_nodes.tsv"

# ----------------------------------------------------------------- PAINT ----
# This figure's appearance, in full and on purpose. It is NOT imported from
# plot_cohort_network: the two drawings have different jobs, the solved one for checking
# the layout and this one for the page, and tying them together means every change to
# either has to suit both. Copy a value across when you want them to match.

MARGIN = 24.0       # points of air around the drawn extent
DPI = 300
BACKGROUND = "white"

# Marker area in points squared. Shared with the packer below, which is why it is here
# and not merely a colour: the radii in cohort_network_nodes.tsv were measured against
# these, so enlarging a marker moves nothing but does crowd what is already placed.
SIZE = {"search": 700, "paper": 110, "no_deposit": 80, "project": 90,
        "field_wheat": 6, "other_run": 6}
SQUARE = {"search"}          # cleared on the diagonal rather than the inscribed radius
RUN_PAD = 0.35               # between two runs inside a cloud, in points
PAD = 2.0                    # between one cloud and the next, for the overlap report

FILL = {"search": "white", "paper": "#EF946C", "no_deposit": "white",
        "project": "#067BC2", "project_unused": "white",
        "field_wheat": "#41823f", "other_run": "white"}
STROKE = {"search": ("black", 2.0), "paper": ("black", 1.2),
          "no_deposit": ("black", 0.8), "project": ("black", 1.2),
          "project_unused": ("black", 0.8), "field_wheat": ("none", 0.80),
          "other_run": ("black", 0.8)}
MARKER = {"search": "D", "paper": "o", "no_deposit": "o", "project": "o",
          "project_unused": "o", "field_wheat": "o", "other_run": "o"}

EDGE = {
    "search_paper":           dict(colour="black", width=1.2, style="solid"),
    "search_paper_unused":    dict(colour="black",   width=0.8, style="solid"),
    "paper_project":          dict(colour="black", width=1.2, style="solid"),
    "paper_project_unused":   dict(colour="black",  width=0.8, style="solid"),
    "project_run":            dict(colour="0.7",  width=0.1, style="solid"),
}

LABEL = {
    "search": "search",
    "paper": "paper, reads in cohort",
    "no_deposit": "paper, field survey, reads not deposited",
    "project": "BioProject, field wheat inside",
    "project_unused": "BioProject, none",
    "field_wheat": "run, field wheat",
    "other_run": "run, lab, off-host or not RNA-seq",
}
LEGEND = ["search", "paper", "no_deposit", "project", "project_unused",
          "field_wheat", "other_run"]
LEGEND_SIZE = (3.4, 2.4)   # inches
LEGEND_FONT = 12
LEGEND_MS = {"search": 13, "paper": 11, "no_deposit": 11, "project": 10,
             "project_unused": 10, "field_wheat": 6, "other_run": 6}

# The run packer lives in plot_cohort_network and reads its marker sizes from the module
# it is defined in, so point it at this file's sizes. Done here, once and visibly,
# rather than by copying twenty lines of sunflower geometry to change two numbers.
net.SIZE, net.SQUARE = SIZE, SQUARE
# ------------------------------------------------------------- END PAINT ----


def write_legend(path):
    fig, ax = plt.subplots(figsize=LEGEND_SIZE)
    fig.patch.set_facecolor(BACKGROUND)
    ax.set_axis_off()
    ax.legend(handles=[Line2D([], [], ls="", marker=MARKER[k], mfc=FILL[k],
                              mec=STROKE[k][0], mew=STROKE[k][1], ms=LEGEND_MS[k],
                              label=LABEL[k]) for k in LEGEND],
              frameon=False, fontsize=LEGEND_FONT, loc="center", labelspacing=1.0,
              handletextpad=1.0, borderpad=0)
    fig.savefig(path, dpi=DPI, facecolor=BACKGROUND, bbox_inches="tight")
    plt.close(fig)


def read_positions():
    """node -> (angle in radians, radius in points), from the hand-editable table."""
    if not NODES.is_file():
        sys.exit(f"no {NODES.name}. Run plot_cohort_network.py --map to write one.")
    placed = {}
    with open(NODES, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            node = tuple(row["node"].split("|", 1))
            placed[node] = (math.radians(float(row["angle_deg"])),
                            float(row["radius"]))
    return placed


def main():
    data = net.Cohort()
    placed = read_positions()

    missing = [n for n in net.placeable(data) if n not in placed]
    if missing:
        sys.exit(f"{len(missing)} nodes are in the figure but not in {NODES.name}, so "
                 f"they have nowhere to go. Re-run plot_cohort_network.py --map after "
                 f"the cohort changes. First few: "
                 f"{', '.join('|'.join(n) for n in missing[:5])}")
    extra = [n for n in placed if n not in set(net.placeable(data))]
    if extra:
        print(f"note: {len(extra)} rows name nodes the cohort no longer has; ignored: "
              f"{', '.join('|'.join(n) for n in extra[:5])}")

    # One data unit is one point. See the docstring: the table is in points, so this is
    # what makes the numbers in it mean what they say.
    pos = {n: np.array([r * math.cos(a), r * math.sin(a)])
           for n, (a, r) in placed.items() if n in set(net.placeable(data))}
    by_hand = len(pos)

    # Runs are the one thing still computed, because nobody wants to place 1,435 dots.
    for b in data.projects:
        members = data.members[b]
        if members:
            pos.update(net.pack(members, pos[("B", b)], 1.0, "field_wheat",
                                RUN_PAD, net.radius("project"))[0])

    reach = max(np.linalg.norm(p) for p in pos.values()) + MARGIN
    side = 2 * reach / 72                      # inches, so that a point is a point
    fig, ax = plt.subplots(figsize=(side, side))
    fig.patch.set_facecolor(BACKGROUND)
    ax.set_facecolor(BACKGROUND)
    ax.set_axis_off()
    ax.set_xlim(-reach, reach)
    ax.set_ylim(-reach, reach)
    ax.set_aspect("equal")

    g = data.graph()
    for group, drawn in (
            ("project_run", data.run_edges),
            ("paper_project",
             [e for e in data.project_edges if e[1][1] in data.used]),
            ("paper_project_unused",
             [e for e in data.project_edges if e[1][1] not in data.used]),
            ("search_paper",
             [e for e in data.paper_edges if e[1][1] in data.depositing]),
            ("search_paper_unused",
             [e for e in data.paper_edges if e[1][1] not in data.depositing])):
        if drawn:
            style = EDGE[group]
            nx.draw_networkx_edges(g, pos, ax=ax, edgelist=drawn,
                                   edge_color=style["colour"], width=style["width"],
                                   style=style["style"])

    def scatter(nodes, kind, zorder):
        if nodes:
            ax.scatter(*np.array([pos[n] for n in nodes]).T,
                       s=SIZE.get(kind, SIZE["project"]), c=FILL[kind],
                       marker=MARKER[kind], ec=STROKE[kind][0],
                       lw=STROKE[kind][1], zorder=zorder)

    shown = [c for c in data.cohort if c["BioProject"] in data.used]
    scatter([("R", c["run"]) for c in shown if not net.field_wheat(c)], "other_run", 2)
    scatter([("R", c["run"]) for c in shown if net.field_wheat(c)], "field_wheat", 2)
    scatter([("B", b) for b in data.projects if b in data.used], "project", 3)
    scatter([("B", b) for b in data.projects if b not in data.used],
            "project_unused", 3)
    scatter([("P", p) for p in data.papers if p not in data.depositing],
            "no_deposit", 4)
    scatter([("P", p) for p in data.papers if p in data.depositing], "paper", 4)
    scatter([data.search], "search", 5)

    fig.tight_layout()
    out = net.FIGDIR / "cohort_network_manual.png"
    fig.savefig(out, dpi=DPI, facecolor=BACKGROUND)
    plt.close(fig)
    legend_out = net.FIGDIR / "cohort_network_manual_legend.png"
    write_legend(legend_out)

    # No placement checks are run, but the reader should still be told what they got.
    overlaps = 0
    clouds = {("B", b): (net.pack(data.members[b], np.zeros(2), 1.0, "field_wheat",
                                 RUN_PAD, net.radius("project"))[1] + PAD)
              for b in data.projects}
    nodes = list(clouds)
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            # A tenth of a point of tolerance: the table rounds angle and radius to
            # one decimal, so two clouds the solver left exactly touching come back
            # fractionally inside each other and would read as a fault.
            if np.linalg.norm(pos[a] - pos[b]) < clouds[a] + clouds[b] - 0.1:
                overlaps += 1
    print(f"wrote {out}\nwrote {legend_out}")
    print(f"{side * 72:.0f} pt square, {by_hand} nodes placed by hand, "
          f"{len(data.cohort):,} runs auto-placed around them")
    print(f"{overlaps} pair(s) of BioProject clouds overlap at these positions"
          + ("" if overlaps else "; nothing was moved to achieve that"))


if __name__ == "__main__":
    main()
