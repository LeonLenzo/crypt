#!/usr/bin/env python3
"""Draw the cohort network from hand-placed node positions. No solver.

    python release/01_cohort/figures/network/plot_network.py
    python release/01_cohort/figures/network/plot_network.py --declared

Reads  network_nodes.tsv   one row per node: angle_deg and radius
       data/outputs/supp_data_s1_studies.tsv, supp_data_s2_cohort.tsv,
       data/curation/screen_fulltext.tsv, relations.tsv   via solve_network
       data/curation/declared_pathogens.tsv   --declared: what each deposit named
Writes network.png, network_legend.png
       network_declared.png, network_declared_legend.png   --declared

TWO DRAWINGS, ONE SET OF POSITIONS. Nothing below moves a node: --declared repaints the
outer ring and leaves out the papers that deposited nothing, so the two figures can be
read against each other. What it changes and why is in the DECLARED PAINT block.

The companion to solve_network.py, which solves the layout: that one spends
minutes deciding where 31 BioProjects go so their edges stay radial and clear each
other's samples. This one does none of that. Every paper and BioProject sits exactly
where `network_nodes.tsv` says, and the only thing computed is the cloud of runs
around each BioProject, because a sample has no position worth arguing about.

Get a starting table from `solve_network.py --map`, which writes the solver's
answer in this format, then move whatever you want. Nothing here checks your placement:
clouds may overlap and edges may run through them. That is the trade. A solver that
cannot be overruled is worse than one you can ignore.

WHY THE FIGURE IS SIZED FROM ITS CONTENTS. The table is in points, because a point is
what marker sizes and clearances are in, and the only way a radius of 236 points means
236 points on the page is for one data unit to be one point. That fixes the canvas at
the extent of what is drawn rather than at a chosen width, which is why there is no
FIGSIZE here.
"""

import argparse
import collections
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
import solve_network as net

NODES = net.FIG / "network_nodes.tsv"

# ----------------------------------------------------------------- PAINT ----
# This figure's appearance, in full and on purpose. It is NOT imported from
# solve_network: the two drawings have different jobs, the solved one for checking
# the layout and this one for the page, and tying them together means every change to
# either has to suit both. Copy a value across when you want them to match.

MARGIN = 24.0       # points of air around the drawn extent
DPI = 300
BACKGROUND = "white"

# Marker area in points squared. Shared with the packer below, which is why it is here
# and not merely a colour: the radii in network_nodes.tsv were measured against
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

# Three different nodes are a white circle with a black ring, and on the figure they
# are told apart by which ring they sit on rather than by the mark. The key cannot
# reproduce that, so each label names its own level: paper, BioProject, run.
LABEL = {
    "search": "search",
    "paper": "paper, reads in the cohort",
    "no_deposit": "paper, field survey, no usable deposit",
    "project": "BioProject, field wheat inside",
    "project_unused": "BioProject, no field wheat inside",
    "field_wheat": "run, field wheat",
    "other_run": "run, lab, off-host or not RNA-seq",
}
LEGEND = ["search", "paper", "no_deposit", "project", "project_unused",
          "field_wheat", "other_run"]
LEGEND_SIZE = (3.4, 2.4)   # inches
LEGEND_FONT = 12
LEGEND_MS = {"search": 13, "paper": 11, "no_deposit": 11, "project": 10,
             "project_unused": 10, "field_wheat": 6, "other_run": 6}

# The run packer lives in solve_network and reads its marker sizes from the module
# it is defined in, so point it at this file's sizes. Done here, once and visibly,
# rather than by copying twenty lines of sunflower geometry to change two numbers.
net.SIZE, net.SQUARE = SIZE, SQUARE
# ------------------------------------------------------------- END PAINT ----


# -------------------------------------------------------- DECLARED PAINT ----
# The variant drawn by --declared. Same positions, same runs, different question of the
# outer ring: not "is there field wheat in here", which the presence of a sample cloud
# already answers, but "what did the study say it was looking at". That frees the fill
# the plain figure spends on a redundancy.
#
# Three channels on the BioProject, each carrying a different fact:
#   fill    the declared pathogen, from data/curation/declared_pathogens.tsv
#   shape   circle if the deposit holds any field run, square if it is lab only
#   stroke  heavy if the deposit holds field wheat, light if not (the same thing the
#           absence of a cloud says; kept as reinforcement, not as the only cue)
#
# FIELD AND LAB ARE MARKED AT BOTH LEVELS, because neither level can say it alone and
# they cover different runs. Four deposits hold both kinds: a square there would be
# false and a circle alone says nothing about the lab runs inside, so the RUN is marked
# too. But 69 of the 86 lab runs are in the four LAB-ONLY deposits, whose clouds are
# not drawn at all, so a run mark alone would show 17 of 86 and silently drop the rest.
# A circle therefore means "this deposit holds field material", never "all field", and
# the caption says so.
#
# The run mark is a FILL and not a shape. A run marker is 6 points squared, about a
# millimetre across in print, and shape does not survive that: the map learned the same
# thing when diamonds for country centroids had to be swapped for circles.

DECLARATIONS = net.MODULE / "data" / "curation" / "declared_pathogens.tsv"

# Chosen around what is already spoken for: the paper salmon and the run green. Pst
# keeps the blue the plain figure gave every BioProject, so the fifteen commonest nodes
# do not change colour between the two drawings.
PATHOGEN_FILL = {
    "Puccinia striiformis f. sp. tritici": "#067BC2",
    "Puccinia triticina": "#B2182B",
    "Puccinia graminis f. sp. tritici": "#7B3294",
    "Fusarium graminearum": "#01665E",
    "Blumeria graminis": "#DFC27D",
}
# Not grey. These nine deposits declared something, it just was not a fungus, and they
# are the only part of the cohort that did not condition on fungal disease. Grey would
# read as a missing value and lose exactly the contrast worth drawing.
UNDECLARED_FILL = "#CC79A7"
UNDECLARED_LABEL = "virus survey, no fungus declared"

# A square of the same area reads heavier than a circle, so it is trimmed.
PROJECT_MARKER = {"field": "o", "lab": "s"}
PROJECT_SIZE = {"field": SIZE["project"], "lab": SIZE["project"] * 0.87}

# The three run states, which the plain figure pools into two. Off-host keeps the plain
# figure's hollow mark; lab is filled dark, because at this size a filled mark against
# a hollow one is the only difference that reads.
RUN_FILL = {"field_wheat": FILL["field_wheat"], "lab_run": "#3d3d3d",
            "offhost_run": "white"}
RUN_STROKE = {"field_wheat": ("none", 0.0), "lab_run": ("none", 0.0),
              "offhost_run": ("black", 0.8)}
RUN_LABEL = {"field_wheat": "run, field wheat",
             "lab_run": "run, lab material",
             "offhost_run": "run, field but off-host or not RNA-seq"}

DECLARED_LEGEND_SIZE = (4.0, 4.4)   # inches; two blocks, so taller than the plain key
# ---------------------------------------------------------- END DECLARED ----


def declarations():
    """BioProject -> the curated declaration, as a list of pathogen names.

    Read from the curation table and not from the cohort's own `declared_pathogens`
    column, which is that table inherited onto every run. The declaration is a property
    of the study, adjudicated once per deposit against a quoted sentence, and reading it
    per run would invite drawing it per run.
    """
    out = {}
    with open(DECLARATIONS, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            named = [p.strip() for p in row["declared_pathogens"].split(";") if p.strip()]
            out[row["BioProject"]] = named
    return out


def project_settings(data):
    """BioProject -> 'field' if any run in it came out of a field, else 'lab'.

    Deliberately not a majority: four deposits hold both, all of them field surveys with
    a handful of lab runs alongside, and their lab runs are already drawn inside the
    cloud. Shape here is about what kind of study made the deposit.
    """
    kinds = collections.defaultdict(set)
    for c in data.cohort:
        kinds[c["BioProject"]].add(c["setting"])
    return {b: ("field" if "field" in s else "lab") for b, s in kinds.items()}


def italic(name):
    """A species name as mathtext, with the rank abbreviations left roman."""
    parts, out = name.split(" f. sp. "), []
    for part in parts:
        out.append(r"$\it{" + part.replace(" ", r"\ ") + "}$")
    return " f. sp. ".join(out)


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


def write_declared_legend(path, shown):
    """Two blocks, because the key now has two independent dimensions.

    `shown` is the declarations actually drawn, in the order they are to be listed, so
    the key cannot name a pathogen the figure does not contain.
    """
    structure = [Line2D([], [], ls="", marker=MARKER[k], mfc=FILL[k], mec=STROKE[k][0],
                        mew=STROKE[k][1], ms=LEGEND_MS[k], label=LABEL[k])
                 for k in ("search", "paper")]
    structure += [
        Line2D([], [], ls="", marker="o", mfc=RUN_FILL[k], mec=RUN_STROKE[k][0],
               mew=max(RUN_STROKE[k][1], 0.8), ms=LEGEND_MS["field_wheat"],
               label=RUN_LABEL[k])
        for k in ("field_wheat", "lab_run", "offhost_run")]
    structure += [
        Line2D([], [], ls="", marker=PROJECT_MARKER["field"], mfc="white", mec="black",
               mew=1.2, ms=LEGEND_MS["project"],
               label="BioProject, holds field material"),
        Line2D([], [], ls="", marker=PROJECT_MARKER["lab"], mfc="white", mec="black",
               mew=1.2, ms=LEGEND_MS["project"], label="BioProject, lab only"),
        Line2D([], [], ls="", marker=PROJECT_MARKER["field"], mfc="white", mec="black",
               mew=0.6, ms=LEGEND_MS["project"],
               label="light outline: no field wheat inside"),
    ]
    swatches = []
    for named in shown:
        fills = [PATHOGEN_FILL[p] for p in named] or [UNDECLARED_FILL]
        label = ", ".join(italic(p) for p in named) or UNDECLARED_LABEL
        swatches.append(Line2D([], [], ls="", marker="o", fillstyle="full"
                               if len(fills) == 1 else "left",
                               mfc=fills[0], mfcalt=fills[-1], mec="black", mew=1.2,
                               ms=LEGEND_MS["project"], label=label))

    # ONE legend, not two. Two figure-level legends have to be anchored by hand, and
    # the lower block climbed into the upper one the moment the pathogen list changed
    # length. A blank entry separates the groups and a handle-less entry heads the
    # second, so matplotlib lays out a single column and nothing can collide.
    blank = Line2D([], [], ls="", marker="", label=" ")
    header = Line2D([], [], ls="", marker="", label="fill: declared pathogen")

    fig, ax = plt.subplots(figsize=DECLARED_LEGEND_SIZE)
    fig.patch.set_facecolor(BACKGROUND)
    ax.set_axis_off()
    ax.legend(handles=structure + [blank, header] + swatches, frameon=False,
              fontsize=LEGEND_FONT, loc="center", labelspacing=1.0,
              handletextpad=1.0, borderpad=0, alignment="left")
    fig.savefig(path, dpi=DPI, facecolor=BACKGROUND, bbox_inches="tight")
    plt.close(fig)


def read_positions():
    """node -> (angle in radians, radius in points), from the hand-editable table."""
    if not NODES.is_file():
        sys.exit(f"no {NODES.name}. Run solve_network.py --map to write one.")
    placed = {}
    with open(NODES, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            node = tuple(row["node"].split("|", 1))
            placed[node] = (math.radians(float(row["angle_deg"])),
                            float(row["radius"]))
    return placed


def main(declared=False):
    data = net.Cohort()
    placed = read_positions()
    # --declared drops the papers that deposited nothing. They are a real attrition and
    # the plain figure is where they are drawn; here they would be 26 unexplained white
    # nodes in the middle of a figure whose subject is what the deposits declared.
    drawn_papers = ([p for p in data.papers if p in data.depositing] if declared
                    else list(data.papers))

    missing = [n for n in net.placeable(data) if n not in placed]
    if missing:
        sys.exit(f"{len(missing)} nodes are in the figure but not in {NODES.name}, so "
                 f"they have nowhere to go. Re-run solve_network.py --map after "
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
             [] if declared else
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
    if declared:
        def run_state(c):
            return ("field_wheat" if net.field_wheat(c)
                    else "lab_run" if c["setting"] == "lab" else "offhost_run")
        for state in ("offhost_run", "lab_run", "field_wheat"):
            nodes = [("R", c["run"]) for c in shown if run_state(c) == state]
            if nodes:
                ax.scatter(*np.array([pos[n] for n in nodes]).T,
                           s=SIZE["field_wheat"], c=RUN_FILL[state], marker="o",
                           ec=RUN_STROKE[state][0], lw=RUN_STROKE[state][1], zorder=2)
    else:
        scatter([("R", c["run"]) for c in shown if not net.field_wheat(c)],
                "other_run", 2)
        scatter([("R", c["run"]) for c in shown if net.field_wheat(c)],
                "field_wheat", 2)

    listed = []
    if declared:
        named_by = declarations()
        kind_by = project_settings(data)
        # Drawn one at a time, because fill, shape and outline weight are three
        # per-deposit facts and scatter takes one marker for the whole call. 31 nodes.
        order = sorted(set(map(tuple, named_by.values())),
                       key=lambda t: (-sum(1 for v in named_by.values()
                                           if tuple(v) == t), t))
        listed = [list(t) for t in order]
        for b in data.projects:
            named = named_by.get(b, [])
            fills = [PATHOGEN_FILL[p] for p in named] or [UNDECLARED_FILL]
            kind = kind_by[b]
            ax.plot(*pos[("B", b)], ls="", marker=PROJECT_MARKER[kind],
                    fillstyle="full" if len(fills) == 1 else "left",
                    mfc=fills[0], mfcalt=fills[-1], mec="black",
                    mew=1.2 if b in data.used else 0.6,
                    ms=math.sqrt(PROJECT_SIZE[kind]), zorder=3)
    else:
        scatter([("B", b) for b in data.projects if b in data.used], "project", 3)
        scatter([("B", b) for b in data.projects if b not in data.used],
                "project_unused", 3)

    scatter([("P", p) for p in drawn_papers if p not in data.depositing],
            "no_deposit", 4)
    scatter([("P", p) for p in drawn_papers if p in data.depositing], "paper", 4)
    scatter([data.search], "search", 5)

    fig.tight_layout()
    stem = "network_declared" if declared else "network"
    out = net.FIG / f"{stem}.png"
    fig.savefig(out, dpi=DPI, facecolor=BACKGROUND)
    plt.close(fig)
    legend_out = net.FIG / f"{stem}_legend.png"
    if declared:
        write_declared_legend(legend_out, listed)
    else:
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
    print(f"{side * 72:.0f} pt square, {by_hand} nodes placed, {len(drawn_papers)} "
          f"papers drawn, {sum(len(v) for v in data.members.values()):,} runs "
          f"auto-placed around them")
    if declared:
        kinds = collections.Counter(project_settings(data)[b] for b in data.projects)
        print(f"{kinds['field']} field-survey deposits, {kinds['lab']} lab-only; "
              f"{len(listed)} declarations in the key")
    print(f"{overlaps} pair(s) of BioProject clouds overlap at these positions"
          + ("" if overlaps else "; nothing was moved to achieve that"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--declared", action="store_true",
                    help="the variant: BioProjects coloured by declared pathogen and "
                         "shaped by field or lab, papers with no deposit left out")
    main(**vars(ap.parse_args()))
