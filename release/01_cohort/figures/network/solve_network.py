#!/usr/bin/env python3
"""Draw the cohort as a network of papers, BioProjects and runs.

    python release/01_cohort/figures/network/solve_network.py          # layout if stale, draw
    python release/01_cohort/figures/network/solve_network.py --render # draw from the cache
    python release/01_cohort/figures/network/solve_network.py --layout # recompute and cache
    python release/01_cohort/figures/network/solve_network.py --report # measure, draw nothing

Reads  data/outputs/supp_data_s1_studies.tsv   studies, their verdict and BioProjects
       data/outputs/supp_data_s2_cohort.tsv     the run-level cohort
       data/curation/screen_fulltext.tsv        what the full text said each study did
       data/curation/relations.tsv                authored: who generated or reanalysed
       data/curation/pinned_nodes.tsv             authored: nodes placed by hand
Writes network_solved.png        the figure, no legend
       network_solved_legend.png the legend on its own, to place in assembly
       network_layout.json       the computed positions, so drawing is cheap
       network_nodes.tsv         --map: every placeable node and where it sits
       network_solved_map.png    --map: the same figure with those nodes named

This is the SOLVER, and its output is a check on the layout rather than the page figure.
The page figure is plot_network.py, which draws from network_nodes.tsv and no solver.

WHY THIS IS SPLIT. Solving the layout takes minutes, because every candidate position for
a BioProject is tested against the other clouds AND against every edge already drawn.
Drawing takes under a second. Changing a colour used to pay the layout cost anyway, so
the two are separated and the positions are cached.

The cache key covers everything geometric: the GEOMETRY block below, and the size and
mtime of the three input tables. Change a colour and `--render` reuses the cache; change
a marker size, a clearance or the data and the cache is stale and the layout is redone.
PAINT never enters the key, GEOMETRY always does, so the split has to be respected when
adding a constant: a size belongs in GEOMETRY even though it is also an appearance.

WHAT THE LAYOUT IS DOING, since none of it is obvious:

Each run hangs off the BioProject that holds it, each BioProject off the papers that cite
it, and every included paper off a single node standing for the search.

A paper is drawn on one of two grounds. Orange: its reads are in the cohort. Grey: its
full text reports sequencing of field material, so it is a survey this cohort would have
used, and it deposited nothing usable. A paper that deposited nothing and ran no field
sequencing is a different kind of study and is left out.

Only a BioProject holding field wheat shows its samples. The rest are drawn as a hollow
node: their runs are still in the cohort and still counted, but spelling out a cloud of
lab or off-host samples spends a quarter of the canvas on material the analysis does not
use.

Four things the plain spring layout gets wrong, and what is done about each:

The search node holds the picture together. Without it the graph is a dozen unconnected
components, which a spring layout drives to the rim.

The layout is stopped at 60 iterations. Run to convergence it collapses each cloud onto
its BioProject.

A spring layout places dimensionless points and knows nothing about marker size, so the
markers drawn on them overlap. Runs are packed into a cloud around their BioProject and
papers into rings around the search node, then those clouds are relaxed apart as bodies.

Edges have to stay in open space. Testing only that a cloud's centre clears its
neighbours leaves a deposit sitting legitimately beyond a bigger one with its edge drawn
straight through that neighbour's samples, which is the thing a reader follows.
"""

import argparse
import collections
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

MODULE = Path(__file__).resolve().parents[2]    # release/01_cohort
FIG = Path(__file__).resolve().parent           # this figure's own directory
OUTPUTS = MODULE / "data" / "outputs"
CURATION = MODULE / "data" / "curation"
STUDIES = OUTPUTS / "supp_data_s1_studies.tsv"
COHORT = OUTPUTS / "supp_data_s2_cohort.tsv"
SCREEN = CURATION / "screen_fulltext.tsv"
RELATIONS = CURATION / "relations.tsv"
PINNED = CURATION / "pinned_nodes.tsv"
NODE_TABLE = FIG / "network_nodes.tsv"
CACHE = FIG / "network_layout.json"

# -------------------------------------------------------------- GEOMETRY ----
# Everything here changes where things sit, so everything here is in the cache key.

FIGSIZE = 8.5          # inches, square

# Marker area in points squared. Size is geometric: the clearances are computed from it.
SIZE = {"search": 700, "paper": 110, "no_deposit": 80, "project": 90,
        "field_wheat": 4, "other_run": 4}
# Markers that need clearing on their diagonal rather than their inscribed radius.
SQUARE = {"search"}

PAD = 2.0          # between one cloud and the next, in points
RUN_PAD = 0.35     # between two runs inside a cloud
PAPER_PAD = 0.8    # between two papers inside the search cloud

# Area per paper, as a multiple of the tightest packing. At 1.0 each tier is squeezed
# into a band one marker deep and the spiral degenerates into a ring, which reads as a
# drawn wheel rather than a cloud.
PAPER_SPREAD = 3.0

# BioProjects sharing a citing paper all want one bearing; fan them around it or they
# queue along it and the edge to the far ones runs through the near.
PROJECT_FAN = 0.30      # angular step between them, radians
PROJECT_FAN_MAX = 0.95  # half-width the fan may not exceed

# How far a BioProject may swing off its paper's bearing to avoid being pushed out past a
# bigger neighbour. Taken only when it brings the cloud closer in by more than its own
# radius, so the bearing is given up only for a real gain.
PROJECT_SWING = 0.15

# Extra radius for named BioProjects, as a multiple of the nearest they could sit. The
# 550-run deposit is wide enough that anything sharing its sector has to be drawn across
# it; sending it further out leaves that sector open. One long edge is a better trade
# than an edge through a neighbour's samples, and naming it here makes it a stated choice
# about one deposit rather than something the placement inferred.
PROJECT_PUSH = {}

LAYOUT_K = 0.12          # spring layout repulsion
LAYOUT_ITERATIONS = 60   # see the docstring; converging ruins it
CONCORD_ROUNDS = 8       # rounds of settling papers and BioProjects onto each other
UNCROSS_ROUNDS = 6
GOLDEN = math.pi * (3 - math.sqrt(5))
# ---------------------------------------------------------- END GEOMETRY ----

# ----------------------------------------------------------------- PAINT ----
# Appearance only. Nothing here moves anything, so none of it invalidates the cache.

DPI = 300
BACKGROUND = "white"
LEGEND_ON_FIGURE = False

FILL = {"search": "white", "paper": "#eb6834", "no_deposit": "0.9",
        "project": "#2a78d6", "project_unused": "white",
        "field_wheat": "#4a9b34", "other_run": "white"}
STROKE = {"search": ("black", 2.0), "paper": ("black", 1.2),
          "no_deposit": ("black", 1.2), "project": ("black", 1.0),
          "project_unused": ("black", 1.0), "field_wheat": ("none", 0.0),
          "other_run": ("black", 0.2)}
MARKER = {"search": "D", "paper": "o", "no_deposit": "o", "project": "o",
          "project_unused": "o", "field_wheat": "o", "other_run": "o"}

# `style` is a matplotlib linestyle. Reanalysis edges are dashed because they are a
# different claim, not a weaker one: the paper used data someone else deposited, so the
# edge has no reason to be short and will cross things. That is the honest picture of
# reuse and the reason the default figure leaves them out.
EDGE = {
    "search_paper":           dict(colour="black", width=0.6, style="solid"),
    "search_paper_unused":    dict(colour="0.5",   width=0.3, style="solid"),
    "paper_project":          dict(colour="black", width=0.6, style="solid"),
    "paper_project_unused":   dict(colour="0.88",  width=0.3, style="solid"),
    "paper_project_reanalysed": dict(colour="0.45", width=0.5, style=(0, (4, 2))),
    "project_run":            dict(colour="0.80",  width=0.1, style="solid"),
}

LABEL = {
    "search": "search",
    "paper": "paper, reads in cohort",
    "no_deposit": "paper, field survey, reads not deposited",
    "project": "BioProject, field wheat inside",
    "project_unused": "BioProject, none",
    "field_wheat": "run, field wheat",
    "other_run": "run, lab, off-host or not RNA-seq",
    "reanalysed": "edge: paper reanalysed this deposit (curated)",
}
LEGEND = ["search", "paper", "no_deposit", "project", "project_unused",
          "field_wheat", "other_run"]
LEGEND_SIZE = (3.4, 2.4)   # inches
LEGEND_FONT = 12
# Marker diameters in the legend, in points. Set here rather than taken from SIZE: the
# search node and a run differ 175-fold in area, which is right on the figure and
# illegible in a key.
LEGEND_MS = {"search": 13, "paper": 11, "no_deposit": 11, "project": 10,
             "project_unused": 10, "field_wheat": 6, "other_run": 6}

# Curvature of the reanalysis edges. Drawn straight they are chords across a layout that
# deliberately ignores them, so they cut over everything and read as noise. Bowing them
# all the same way makes them a layer the eye can follow or discount as a group, and
# keeps them off the straight primary edges they would otherwise run along.
REANALYSIS_ARC = 0.22

# The --map overlay: node names drawn on the figure, to read against the node table.
MAP_FONT = 4.5
MAP_COLOUR = "0.25"
# ------------------------------------------------------------- END PAINT ----


# ------------------------------------------------------------------ data ----

def read_relations():
    """(doi, BioProject) -> relation, from the authored overlay. Comments start with #.

    Reanalysis has to be read from here and cannot be inferred. The tempting inference,
    "named in the text but not in the deposit section", is what the figure used first,
    and it labelled 17 edges reanalysis of which only 6 had ever been adjudicated: the
    other 11 are papers that cite an accession in passing. A mention is the weakest
    relation there is, and drawing it as reanalysis repeats in the figure the mistake
    that bioprojects_own made in the table.
    """
    if not RELATIONS.is_file():
        return {}
    lines = [l for l in RELATIONS.read_text(encoding="utf-8").splitlines()
             if l and not l.startswith("#")]
    return {(r["doi"].strip().lower(), r["BioProject"].strip()): r["relation"].strip()
            for r in csv.DictReader(lines, delimiter="\t")}


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def accessions(field):
    return [a.strip() for a in (field or "").replace(";", ",").split(",") if a.strip()]


def field_wheat(run):
    """A field-collected Triticum sample the classifier can actually use.

    Triticale is a wheat-rye hybrid and is counted off-host, with the barley and the rye.
    The assay test matters as much as the host: a deposit's WGS reference genomes and
    amplicon diagnostics are no more usable here than an off-host run.
    """
    return (run.get("assay_in_scope", "1") == "1" and run["setting"] == "field"
            and run["host"].startswith("Triticum"))


def surveyed_in_field(screen):
    """The full text reports sequencing, and the material came out of a field.

    This is the whole basis on which a paper with no deposit is drawn. Without it the
    grey set collects anything that passed screening and happened not to deposit,
    including first-report notes that sequenced nothing at all.
    """
    return screen.get("seq_in_body") == "1" and screen.get("field") == "1"


class Cohort:
    """The tables, and the sets the figure is built from."""

    def __init__(self):
        studies = read(STUDIES)
        self.cohort = read(COHORT)
        screened = {s["doi"]: s for s in read(SCREEN)}

        self.runs = collections.Counter(c["BioProject"] for c in self.cohort)
        edges = []
        for s in studies:
            if s["verdict"] != "keep":
                continue
            own = set(accessions(s["bioprojects_own"]))
            # Union, not just the scanned list: a relations.tsv row asserts a link the
            # text scan could not see, which is the only evidence available for a paper
            # whose body could not be retrieved.
            for bp in sorted(set(accessions(s["bioprojects"])) | own):
                if bp in self.runs:
                    edges.append((s["doi"], bp,
                                  "generated" if bp in own else "reanalysed"))
        # The layout, and every metric reported against it, uses the generated edges
        # only. Reanalysis is carried separately and drawn over the top on request, so
        # asking for it does not move anything and does not need a second solve.
        self.edges = [e for e in edges if e[2] == "generated"]
        relation = read_relations()
        self.reanalysed = [e for e in edges if e[2] == "reanalysed"
                           and relation.get((e[0].lower(), e[1])) == "reanalysed"]
        self.mentioned = [e for e in edges if e[2] == "reanalysed"
                          and relation.get((e[0].lower(), e[1])) != "reanalysed"]

        self.projects = sorted(self.runs, key=lambda b: -self.runs[b])
        self.depositing = {e[0] for e in self.edges}
        degree = collections.Counter(e[0] for e in self.edges)
        self.papers = sorted(
            (s["doi"] for s in studies if s["verdict"] == "keep"
             and (s["doi"] in self.depositing
                  or surveyed_in_field(screened.get(s["doi"], {})))),
            key=lambda d: (d not in self.depositing, -degree[d]))

        self.search = ("S", "search")
        self.sourced = [("P", p) for p in self.papers if p in self.depositing]
        self.dead_end = [("P", p) for p in self.papers if p not in self.depositing]
        self.cited = {p: [b for q, b, _ in edges if q == p] for p in self.depositing}

        # A BioProject counts as used when it holds at least one field-wheat run.
        self.used = {c["BioProject"] for c in self.cohort if field_wheat(c)}
        self.members = collections.defaultdict(list)
        for c in self.cohort:
            if c["BioProject"] in self.used:
                self.members[c["BioProject"]].append(("R", c["run"]))

        # A BioProject with no generating paper has no edge to aim at. It still has to go
        # somewhere, so it keeps the bearing of whoever cites it, and shows no edge.
        self.fallback = collections.defaultdict(list)
        for s in studies:
            if s["verdict"] != "keep" or s["doi"] not in self.depositing:
                continue
            for b in accessions(s["bioprojects"]):
                if b in self.runs:
                    self.fallback[b].append(s["doi"])

        self.paper_edges = [(self.search, ("P", p)) for p in self.papers]
        self.project_edges = [(("P", p), ("B", b)) for p, b, _ in self.edges]
        # Only between nodes the figure already draws; a reanalysis by a paper that is
        # not in the figure has no node to start from.
        drawn = set(self.papers)
        self.reanalysis_edges = [(("P", p), ("B", b)) for p, b, _ in self.reanalysed
                                 if p in drawn]
        self.run_edges = [(("B", c["BioProject"]), ("R", c["run"]))
                          for c in self.cohort if c["BioProject"] in self.used]

    def graph(self):
        g = nx.Graph()
        g.add_nodes_from(("B", b) for b in self.projects)
        # Reanalysis edges are deliberately absent. They are drawn over the result,
        # but adding them here feeds them to spring_layout, which is steering by the
        # back door: it moved every node and took crossings from 0 to 4. Nothing in the
        # drawing needs them in the graph, since they join nodes that are already in it
        # and are passed to the renderer as an explicit edgelist.
        g.add_edges_from(self.paper_edges + self.project_edges + self.run_edges)
        return g


# -------------------------------------------------------------- geometry ----

def radius(kind):
    """Radius of a marker in points, circumscribed so a square clears on its diagonal.

    Circles take the inscribed radius, which is the real one; charging a run the diagonal
    would inflate every cloud by 40%.
    """
    return math.sqrt(SIZE[kind]) / 2 * (math.sqrt(2) if kind in SQUARE else 1)


def bearings(points, origin):
    delta = np.asarray(points) - origin
    return np.arctan2(delta[:, 1], delta[:, 0])


def circular_mean(angles):
    return math.atan2(np.mean(np.sin(angles)), np.mean(np.cos(angles)))


def clears(a, b, centre, keep_out):
    """True if the segment a->b stays outside a disc of `keep_out` about `centre`."""
    ab = b - a
    span = float(np.dot(ab, ab))
    t = 0.0 if span == 0 else float(np.clip(np.dot(centre - a, ab) / span, 0, 1))
    return float(np.linalg.norm(centre - (a + t * ab))) >= keep_out


def pack(members, centre, scale, kind, pad, hole, bearing=None, spread=1.0, phase=0.0):
    """Lay members out on a sunflower spiral around the node that owns them.

    Successive slots sit a golden angle apart at a radius that holds the area per slot
    constant, so neighbours are evenly spaced at every radius. The spiral starts clear of
    `hole`, a radius in points: the owning marker for the first ring, the previous ring
    for any after it. Given `bearing`, members are assigned to slots by direction rather
    than in order, so each sits on the side its own edges leave towards. `spread`
    multiplies the area each member gets, widening the band. `phase` rotates the spiral,
    which is what keeps one tier of papers off the radial edges of the tier outside it:
    two tiers of the same size otherwise get the same golden-angle sequence and every
    inner member lands exactly under an outer one.

    Returns the positions and the radius of the cloud, in data units.
    """
    spacing = 2 * (radius(kind) + pad) * scale
    own = hole * scale
    hole = (hole + radius(kind) + pad) * scale
    if not members:
        # Nothing to ring, so claim only the owning marker. Returning the hole instead
        # reserved a run's width and RUN_PAD around a node with no runs, which on a
        # hollow BioProject is 8.1pt of clearance around a 4.7pt marker.
        return {}, own
    cloud = math.sqrt(hole ** 2 + spacing ** 2 * len(members) * spread)
    k = np.arange(len(members)) + 0.5
    r = np.sqrt(hole ** 2 + (cloud ** 2 - hole ** 2) * k / len(members))
    theta = k * GOLDEN + phase
    slots = np.c_[r * np.cos(theta), r * np.sin(theta)]
    if bearing is not None:
        # Both sides must be wrapped the same way before they are ranked against each
        # other. Sorting raw bearings, which run -pi to pi, against slot angles taken
        # modulo 2*pi offsets every member by half a turn from where it was asked for.
        slots = slots[np.argsort(np.mod(theta, 2 * np.pi))]
        members = [m for _, m in sorted(zip(np.mod(bearing, 2 * np.pi), members),
                                        key=lambda mb: mb[0])]
    return dict(zip(members, centre + slots)), cloud


def place_outward(origin, angle, cloud, floor, step, alternates=None, anchors=None,
                  obstacles=None, swing=PROJECT_SWING):
    """Put each cloud near its own bearing, as close in as it will go without touching.

    Crowding is mostly resolved by moving a cloud further out rather than sideways: the
    bearing is the point, since a BioProject off its paper's bearing makes an edge that
    runs across the figure instead of onward from the paper.

    Pure outward movement alone stacks a small cloud beyond a big one whenever the two
    share a bearing, and the edge to it then runs the whole length of the big one. So a
    cloud may also swing up to `swing` radians, taking it only where that brings it
    closer in by more than its own radius.

    Largest clouds are placed first, so a big one takes the near slot and the small ones
    fit around it rather than the reverse.
    """
    alternates = alternates or [()] * len(angle)
    anchors = anchors if anchors is not None else [origin] * len(angle)
    # Papers are fixed obstacles. A deposit used to be floored outside the whole paper
    # band, which left a deposit whose own paper sits deep inside the band stranded
    # about 29pt out with nothing in between. Flooring it against its own paper instead
    # needs the other papers to become things it can collide with.
    obstacles = obstacles if obstacles is not None else []
    where = np.zeros((len(angle), 2))
    done, drawn = [], []

    def reach(i, bearing):
        """Smallest radius on this bearing that clears everything already placed.

        Both the cloud and the edge running out to it have to clear, in both directions:
        this cloud must not sit on an edge already drawn, and its own edge must not cross
        a cloud already placed. Testing only the centre, or only one direction, leaves a
        deposit with its edge through a neighbour's samples.
        """
        direction = np.array([math.cos(bearing), math.sin(bearing)])
        # Its own paper is not an obstacle to it: the edge between them is the one line
        # that is supposed to touch both.
        mine = [(q, c) for q, c in obstacles
                if np.linalg.norm(q - anchors[i]) > 1e-9]
        r = floor[i]
        while True:
            centre = origin + r * direction
            if (all(np.linalg.norm(centre - q) >= cloud[i] + c
                    and clears(anchors[i], centre, q, c) for q, c in done)
                    and all(np.linalg.norm(centre - q) >= cloud[i] + c
                            and clears(anchors[i], centre, q, c) for q, c in mine)
                    and all(clears(a, b, centre, cloud[i]) for a, b in drawn)):
                return r, direction
            r += step
            if r > 20 * floor[i]:        # unplaceable; take the bearing and move on
                return floor[i], direction

    for i in np.argsort(-cloud):
        r, direction = reach(i, angle[i])
        for off in [k * swing / 8 * sign for k in range(1, 9) for sign in (1, -1)]:
            swung, swung_dir = reach(i, angle[i] + off)
            if r - swung > cloud[i]:      # worth giving up the bearing for
                r, direction = swung, swung_dir
        # A deposit cited by more than one paper is aimed at the midpoint between them,
        # which is often the one blocked direction. Each citing paper's own bearing is an
        # equally true place for it, so try those too. Sitting beside one of its papers
        # beats sitting deep between both.
        for alt in alternates[i]:
            cand, cand_dir = reach(i, alt)
            if cand < r - cloud[i]:
                r, direction = cand, cand_dir
        # Letting a boxed-in cloud off its bearing entirely was tried and reverted:
        # freeing one deposit moved it early in the placement order, which changed what
        # was left for everything after it and reordered the papers through the settle
        # loop. Two others ended up deeper than the one fixed, crossings went 2 to 4 and
        # the median turn 5 to 10 degrees. Placement is sequential and fed back, so no
        # fix here stays local to one cloud. PROJECT_PUSH is the stated-exception route.
        where[i] = origin + r * direction
        done.append((where[i], cloud[i]))
        drawn.append((anchors[i], where[i]))
    return where


def settle(data, pos, scale, clouds, rounds=CONCORD_ROUNDS):
    """Place the papers and swing each BioProject onto the bearing of its own papers.

    A paper-to-BioProject edge should run outward, roughly along the search-to-paper
    edge, not across the figure. That needs the BioProject at the paper's bearing. An
    earlier version only permuted the projects among the positions the spring layout had
    already chosen, so a project could get no nearer than the closest position that
    happened to exist.

    Here the angle comes from the papers and the radius is the shortest that clears the
    paper ring and the clouds already placed, so spacing varies with how much each
    BioProject holds. Angle and papers depend on each other, so they alternate: place the
    papers from the current bearings, read their angles back, move each project onto the
    mean bearing of the papers citing it, and let any that collide slide outward.

    Returns the paper positions, the outer paper radius, and the project positions.
    """
    origin = np.array(pos[data.search])
    order = list(data.projects)
    angle = bearings([pos[("B", b)] for b in order], origin)
    cloud = np.array([clouds[b] for b in order])
    citing = {b: [p for _, p in data.sourced if b in data.cited[p]]
              or data.fallback.get(b, []) for b in order}
    outer_hole = radius("search")

    def place_papers(angle):
        """Lay the two paper tiers out for the BioProject bearings given."""
        at = dict(zip(order, angle))
        bearing = [circular_mean([at[b] for b in data.cited[p]])
                   for _, p in data.sourced]
        # The grey papers have only the one edge, so their bearing is free. Putting them
        # between the orange bearings keeps a grey paper off the radial edge out to an
        # orange one, which otherwise reads as the two being connected to each other.
        ranked = np.sort(bearing)
        midpoint = np.angle(np.exp(1j * ranked) + np.exp(1j * np.roll(ranked, -1)))
        inner = [midpoint[k % len(midpoint)] for k in range(len(data.dead_end))]
        inside, hole = pack(data.dead_end, origin, scale, "no_deposit", PAPER_PAD,
                            outer_hole, bearing=inner, spread=PAPER_SPREAD,
                            phase=math.pi / max(len(data.dead_end), 1))
        outside, edge = pack(data.sourced, origin, scale, "paper", PAPER_PAD,
                             hole / scale, bearing=bearing, spread=PAPER_SPREAD)
        return {**inside, **outside}, edge

    def floors(outer, placed=None):
        """Nearest each deposit could sit: clear of the whole paper band.

        Flooring against each deposit's OWN paper instead was tried, to close the ~29pt
        of empty space a paper sitting deep inside the band leaves in front of it. It
        was reverted: the deposits come in, but the ring they sat on is what kept the
        figure readable outward. Crossings went 0 to 4, edges through a neighbour's
        samples 1 to 4, and the median turn 5 to 21 degrees. The gap is the price of
        the radial reading, not waste.
        """
        return np.array([(outer + cloud[i] + PAD * scale) * PROJECT_PUSH.get(b, 1.0)
                         for i, b in enumerate(order)])

    placed, outer, where = {}, 0.0, None
    for _ in range(rounds):
        placed, outer = place_papers(angle)
        paper_angle = {p: math.atan2(*reversed(placed[("P", p)] - origin))
                       for _, p in data.sourced}
        # A BioProject nothing placed points at keeps the bearing it already had, rather
        # than an average over an empty set.
        want = np.array([circular_mean([paper_angle[p] for p in citing[b]])
                         if citing[b] else angle[i] for i, b in enumerate(order)])
        together = collections.defaultdict(list)
        for i, b in enumerate(order):
            together[tuple(sorted(citing[b]))].append(i)
        for group in together.values():
            if len(group) > 1:
                half = min(PROJECT_FAN * (len(group) - 1) / 2, PROJECT_FAN_MAX)
                for i, d in zip(sorted(group, key=lambda i: cloud[i]),
                                np.linspace(-half, half, len(group))):
                    want[i] += d

        alternates = [[paper_angle[p] for p in citing[b]] if len(citing[b]) > 1 else []
                      for b in order]
        anchors = [np.array(placed.get(("P", citing[b][0]), origin)) if citing[b]
                   else origin for b in order]
        where = place_outward(origin, want, cloud + PAD * scale, floors(outer),
                              step=0.02 * outer, alternates=alternates, anchors=anchors)
        # Feed the realised bearings back. A cloud that swung to get past a bigger
        # neighbour sits off the bearing its paper was placed on, and the edge between
        # them goes diagonal unless the paper is allowed to follow it.
        angle = bearings(where, origin)

    # Answer the last swing, then place the clouds once more against those final paper
    # positions: clearance was computed against where the papers used to be, and moving
    # them puts the edges back through the samples the placement just avoided.
    placed, outer = place_papers(angle)
    anchors = [np.array(placed.get(("P", citing[b][0]), origin)) if citing[b] else origin
               for b in order]
    where = place_outward(origin, bearings(where, origin), cloud + PAD * scale,
                          floors(outer), step=0.02 * outer, anchors=anchors)
    return placed, outer, {("B", b): xy for b, xy in zip(order, where)}


def separate(pos, radii, carried, scale, iterations=3000):
    """Push the given clouds apart until none overlap, holding each near where it was.

    The spring back to the layout position is annealed to zero. Held constant it reaches
    equilibrium against the push and leaves a few points of overlap unresolved.

    Tangency is not overlap, so the test carries a hundredth of a point of tolerance; a
    float epsilon leaves clouds that have settled against each other reported as
    failures. A cloud's clearance radius bounds it as a circle, so two clouds can
    intersect while no marker in either touches a marker in the other; any residual is
    reported in points to make that judgeable rather than alarming.
    """
    nodes = list(radii)
    p = np.array([pos[n] for n in nodes], dtype=float)
    anchored = p.copy()
    r = np.array([radii[n] for n in nodes])

    for step in range(iterations):
        pairs = cKDTree(p).query_pairs(2 * r.max(), output_type="ndarray")
        i, j = pairs[:, 0], pairs[:, 1]
        delta = p[i] - p[j]
        dist = np.hypot(delta[:, 0], delta[:, 1])
        overlap = r[i] + r[j] - dist
        hit = overlap > 0.01 * scale
        if not hit.any():
            break
        i, j, delta = i[hit], j[hit], delta[hit]
        dist = np.where(dist[hit] == 0, 1e-9, dist[hit])
        push = delta / dist[:, None] * (r[i] + r[j] - dist)[:, None] * 0.5
        move = np.zeros_like(p)
        np.add.at(move, i, push)
        np.add.at(move, j, -push)
        p += 0.5 * move + 0.02 * (1 - step / iterations) ** 2 * (anchored - p)
    else:
        print(f"note: {hit.sum()} cloud bounds still intersect after {iterations} "
              f"iterations, worst by {overlap.max() / scale:.1f} pt")

    for node, old, new in zip(nodes, anchored, p):
        pos[node] = new
        for child in carried.get(node, ()):
            pos[child] = pos[child] + (new - old)


def crossings(pos, pairs):
    """How many of these edges cross each other, counted as segment intersections."""
    def side(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    segments = [(np.asarray(pos[u]), np.asarray(pos[v])) for u, v in pairs]
    total = 0
    for i, (a, b) in enumerate(segments):
        for c, d in segments[i + 1:]:
            if side(a, b, c) * side(a, b, d) < 0 and side(c, d, a) * side(c, d, b) < 0:
                total += 1
    return total


def uncross(pos, papers, watched, rounds=UNCROSS_ROUNDS):
    """Swap papers between their own positions while that removes edge crossings.

    `watched` must list EVERY edge that can cross, the spokes out to the papers as well
    as the edges on to the BioProjects. Counting only the latter made this actively
    harmful: it removed 25 paper-to-BioProject crossings by creating 132 against the
    spokes, tripling the total, because placing a paper on its BioProject's bearing had
    already made spoke crossings almost impossible and the swaps gave that up.

    What it buys over the bearing order alone is the exchanges a barycentre ordering
    cannot reach, where two papers improve together and neither improves alone. Only
    papers move, and only onto each other's positions.
    """
    nodes = [("P", p) for p in papers]
    best = crossings(pos, watched)
    for _ in range(rounds):
        improved = False
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                pos[nodes[i]], pos[nodes[j]] = pos[nodes[j]], pos[nodes[i]]
                now = crossings(pos, watched)
                if now < best:
                    best, improved = now, True
                else:
                    pos[nodes[i]], pos[nodes[j]] = pos[nodes[j]], pos[nodes[i]]
        if not improved:
            break
    return best


def data_per_point(fig, ax):
    """Size of one typographic point in data units, with the axes already scaled."""
    fig.canvas.draw()
    inverse = ax.transData.inverted()
    origin = inverse.transform((0, 0))
    step = inverse.transform((fig.dpi / 72, 0))
    return abs(step[0] - origin[0])


# ----------------------------------------------------------------- cache ----

def placeable(data):
    """The nodes a person can place: the search node, the papers, the BioProjects.

    Runs are deliberately absent. A run has no position of its own worth arguing about;
    it sits in the cloud its BioProject carries, and giving 1,435 of them rows would
    bury the handful of decisions that matter in a file nobody could read.
    """
    return ([data.search] + [("P", p) for p in data.papers]
            + [("B", b) for b in data.projects])


def apply_pins(pos, carried, scale, search):
    """Move hand-placed nodes where they were asked for, carrying their samples.

    Applied after everything else and checked against nothing, so a pin always wins.
    That is the point: the solver optimises what it can measure, and a pin is for the
    reasons it cannot. Anything a pin breaks is reported rather than corrected, because
    silently undoing an instruction is worse than an overlap.
    """
    if not PINNED.is_file():
        return []
    lines = [l for l in PINNED.read_text(encoding="utf-8").splitlines()
             if l and not l.startswith("#")]
    origin = np.array(pos[search])
    moved = []
    for row in csv.DictReader(lines, delimiter="\t"):
        node = tuple(row["node"].strip().split("|", 1))
        if node not in pos:
            print(f"warning: pinned_nodes.tsv names {row['node']}, which is not a node "
                  f"in this figure; ignored")
            continue
        here = np.array(pos[node]) - origin
        angle = (math.radians(float(row["angle_deg"]))
                 if row.get("angle_deg", "").strip()
                 else math.atan2(here[1], here[0]))
        length = (float(row["radius"]) * scale if row.get("radius", "").strip()
                  else float(np.linalg.norm(here)))
        target = origin + length * np.array([math.cos(angle), math.sin(angle)])
        shift = target - np.array(pos[node])
        pos[node] = target
        for child in carried.get(node, ()):
            pos[child] = pos[child] + shift
        moved.append((row["node"], float(np.linalg.norm(shift)) / scale))
    return moved


def cache_key():
    """Everything that changes where things sit: the GEOMETRY block, and the inputs."""
    stamp = {}
    for path in (STUDIES, COHORT, SCREEN, RELATIONS, PINNED):
        if not path.is_file():
            continue
        st = path.stat()
        stamp[path.name] = [st.st_size, int(st.st_mtime)]
    return {
        "inputs": stamp,
        "geometry": {
            "FIGSIZE": FIGSIZE, "SIZE": SIZE,
            "SQUARE": sorted(SQUARE), "PAD": PAD, "RUN_PAD": RUN_PAD,
            "PAPER_PAD": PAPER_PAD, "PAPER_SPREAD": PAPER_SPREAD,
            "PROJECT_FAN": PROJECT_FAN, "PROJECT_FAN_MAX": PROJECT_FAN_MAX,
            "PROJECT_SWING": PROJECT_SWING, "PROJECT_PUSH": PROJECT_PUSH,
            "LAYOUT_K": LAYOUT_K, "LAYOUT_ITERATIONS": LAYOUT_ITERATIONS,
            "CONCORD_ROUNDS": CONCORD_ROUNDS, "UNCROSS_ROUNDS": UNCROSS_ROUNDS,
        },
    }


def save_layout(pos, radii, scale, outer_ring, window):
    CACHE.parent.mkdir(exist_ok=True)
    CACHE.write_text(json.dumps({
        "key": cache_key(),
        "scale": scale,
        "outer_ring": outer_ring,
        "window": list(window),
        "positions": {f"{k[0]}|{k[1]}": [float(v[0]), float(v[1])]
                      for k, v in pos.items()},
        "radii": {f"{k[0]}|{k[1]}": float(v) for k, v in radii.items()},
    }), encoding="utf-8")


def load_layout():
    """Cached positions, or None if absent or computed under different geometry."""
    if not CACHE.is_file():
        return None
    try:
        blob = json.loads(CACHE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    if blob.get("key") != cache_key():
        return None
    unkey = lambda s: tuple(s.split("|", 1))
    blob["positions"] = {unkey(k): np.array(v) for k, v in blob["positions"].items()}
    blob["radii"] = {unkey(k): v for k, v in blob["radii"].items()}
    return blob


def compute_layout(data):
    """Solve the layout. This is the slow half: minutes, not seconds."""
    pos = nx.spring_layout(data.graph(), k=LAYOUT_K, iterations=LAYOUT_ITERATIONS,
                           seed=1)

    # A throwaway figure, only to learn how big a point is in data units. Marker sizes
    # are in points and every clearance is derived from them, so the layout cannot be
    # solved without knowing the mapping.
    fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE))
    ax.set_axis_off()
    corner = np.abs(np.array(list(pos.values()))).max() * 1.06
    ax.set_xlim(-corner, corner)
    ax.set_ylim(-corner, corner)
    ax.set_aspect("equal")
    scale = data_per_point(fig, ax)
    plt.close(fig)

    clouds = {b: pack(data.members[b], np.zeros(2), scale, "field_wheat", RUN_PAD,
                      radius("project"))[1] for b in data.projects}
    placed, outer, projected = settle(data, pos, scale, clouds)
    pos.update(placed)
    pos.update(projected)

    carried, radii = collections.defaultdict(list), {}
    for b in data.projects:
        node = ("B", b)
        carried[node] = data.members[b]
        pos.update(pack(data.members[b], pos[node], scale, "field_wheat", RUN_PAD,
                        radius("project"))[0])
        radii[node] = clouds[b] + PAD * scale
    carried[data.search] = data.dead_end + data.sourced
    radii[data.search] = outer + PAD * scale
    outer_ring = radii[data.search]

    separate(pos, radii, carried, scale)
    watched = data.project_edges + data.paper_edges
    before = crossings(pos, watched)
    after = uncross(pos, [p for p in data.papers if p in data.depositing], watched)
    print(f"edge crossings (BioProject edges and spokes): {before} -> {after}")

    for node, distance in apply_pins(pos, carried, scale, data.search):
        print(f"   pinned {node}, moved {distance:.0f} pt by hand")

    # Placing the BioProjects on their papers' bearings leaves the figure off centre as
    # well as a different size, so fit the window to what will actually be drawn.
    drawn = np.array(list(pos.values()))
    low, high = drawn.min(axis=0), drawn.max(axis=0)
    margin = max(radii.values()) * 0.10
    middle, half = (low + high) / 2, (high - low).max() / 2 + margin
    window = (middle[0] - half, middle[0] + half, middle[1] - half, middle[1] + half)
    return pos, radii, scale, outer_ring, window


# --------------------------------------------------------------- measure ----

def report(data, pos, radii, outer_ring):
    """Everything the layout is trying to achieve, stated as a number."""
    deep = sorted(((np.linalg.norm(pos[("B", b)] - pos[data.search])
                    / (radii[("B", b)] + outer_ring), b) for b in data.projects),
                  reverse=True)
    for ratio, b in deep[:5]:
        print(f"   {b:<14} sits {ratio:.2f}x its nearest possible radius")

    through = []
    for paper, bp, _ in data.edges:
        a, b = pos[("P", paper)], pos[("B", bp)]
        for other in data.projects:
            if other != bp and data.members[other] and not clears(
                    a, b, pos[("B", other)], radii[("B", other)]):
                through.append((paper, bp, other))
    print(f"   edges passing through another BioProject's samples: {len(through)}")
    for paper, bp, other in through[:4]:
        print(f"      {bp} <- {paper} crosses {other}")

    # The figure is meant to be read outward, so report how far from that it ended up:
    # the turn between the search-to-paper edge and each paper-to-BioProject edge.
    turn = []
    for paper, bp, _ in data.edges:
        out = pos[("P", paper)] - pos[data.search]
        on = pos[("B", bp)] - pos[("P", paper)]
        turn.append(abs(math.degrees(math.atan2(out[0] * on[1] - out[1] * on[0],
                                                float(np.dot(out, on))))))
    turn = np.array(turn)
    print(f"paper to BioProject edges turn a median {np.median(turn):.0f} deg off "
          f"radial, worst {turn.max():.0f} deg, {(turn > 45).sum()} of {len(turn)} "
          f"over 45 deg")

    hollow = [b for b in data.projects if b not in data.used]
    unknown = [b for b in hollow if not any(c["setting"] for c in data.cohort
                                            if c["BioProject"] == b)]
    print(f"{len(hollow)} BioProjects drawn hollow, holding no field-wheat run; "
          f"{len(unknown)} of those only because setting is blank: "
          f"{', '.join(unknown) or 'none'}")

    green = sum(1 for c in data.cohort if field_wheat(c))
    shared = sum(1 for b in data.projects
                 if len({e[0] for e in data.edges if e[1] == b}) > 1)
    relation = collections.Counter(e[2] for e in data.edges)
    print(f"{len(data.papers)} papers ({len(data.depositing)} depositing, "
          f"{len(data.papers) - len(data.depositing)} field surveys that deposited "
          f"nothing), {len(data.projects)} BioProjects, {len(data.cohort):,} runs "
          f"({green:,} field wheat, {len(data.cohort) - green} lab or off-host), "
          f"{len(data.edges)} edges ({relation['generated']} generated, "
          f"{relation['reanalysed']} reanalysed), {shared} cited by >1 paper")


# ---------------------------------------------------------------- render ----

def handles(reanalysis=False):
    keys = [Line2D([], [], ls="", marker=MARKER[k], mfc=FILL[k], mec=STROKE[k][0],
                   mew=STROKE[k][1], ms=LEGEND_MS[k], label=LABEL[k]) for k in LEGEND]
    if reanalysis:
        style = EDGE["paper_project_reanalysed"]
        keys.append(Line2D([], [], color=style["colour"], lw=style["width"],
                           ls=style["style"], label=LABEL["reanalysed"]))
    return keys


def write_legend(path, reanalysis=False):
    fig, ax = plt.subplots(figsize=LEGEND_SIZE)
    fig.patch.set_facecolor(BACKGROUND)
    ax.set_axis_off()
    ax.legend(handles=handles(reanalysis), frameon=False, fontsize=LEGEND_FONT,
              loc="center",
              labelspacing=1.0, handletextpad=1.0, borderpad=0)
    fig.savefig(path, dpi=DPI, facecolor=BACKGROUND, bbox_inches="tight")
    plt.close(fig)


def seed_pins(rows):
    """Fill pinned_nodes.tsv with a row per node, if it has none yet.

    An overlay of comments and a header is not a table anyone can use: you cannot see
    what is pinnable without running something. Seeded, every node has a row with its
    solved angle and radius left blank, so pinning is editing a cell rather than
    composing a line. Blank means keep what the solver chose, so a fully seeded file
    changes nothing until a number is typed into it.

    Only ever written when the file holds no rows. Once there are decisions in it, it
    is a curation file like any other and this does not touch it.
    """
    text = PINNED.read_text(encoding="utf-8") if PINNED.is_file() else ""
    body = [l for l in text.splitlines() if l and not l.startswith("#")]
    if len(body) > 1:                      # header plus at least one decision
        return 0
    header = [l for l in text.splitlines() if l.startswith("#")]
    out = header + ["\t".join(["id", "node", "angle_deg", "radius", "reason"])]
    out += ["\t".join([r["id"], r["node"], "", "", ""]) for r in rows]
    PINNED.write_text("\n".join(out) + "\n", encoding="utf-8")
    return len(rows)


def write_node_table(data, pos, scale):
    """Every placeable node, named, with where it currently sits.

    This is the companion to pinned_nodes.tsv: the `node` column is exactly the key that
    file expects, and `angle_deg` and `radius` are in the units it takes, so moving
    something is copy the row, change a number. Samples are not listed, because a run
    has no position worth arguing about; it rides in its BioProject's cloud.
    """
    pinned = set()
    if PINNED.is_file():
        lines = [l for l in PINNED.read_text(encoding="utf-8").splitlines()
                 if l and not l.startswith("#")]
        pinned = {r["node"].strip() for r in csv.DictReader(lines, delimiter="\t")}

    origin = np.array(pos[data.search])
    rows = []
    for n, node in enumerate(placeable(data)):
        here = np.array(pos[node]) - origin
        key = f"{node[0]}|{node[1]}"
        rows.append({
            "id": LABEL_FOR(data, node, n),
            "node": key,
            "kind": {"S": "search", "P": "paper", "B": "bioproject"}[node[0]],
            "angle_deg": f"{math.degrees(math.atan2(here[1], here[0])) % 360:.1f}",
            "radius": f"{np.linalg.norm(here) / scale:.1f}",
            "runs": len(data.members.get(node[1], [])) if node[0] == "B" else "",
            "pinned": "yes" if key in pinned else "",
        })
    with open(NODE_TABLE, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


def LABEL_FOR(data, node, n):
    """Short name for a node. BioProjects use their accession, which is already short
    and meaningful; papers get a number, because a DOI on the figure is unreadable."""
    if node[0] == "S":
        return "search"
    if node[0] == "B":
        return node[1]
    return f"P{data.papers.index(node[1]) + 1:02d}"


def render(data, pos, window, reanalysis=False, labels=False):
    fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE))
    fig.patch.set_facecolor(BACKGROUND)
    ax.set_facecolor(BACKGROUND)
    ax.set_axis_off()
    ax.set_xlim(window[0], window[1])
    ax.set_ylim(window[2], window[3])
    ax.set_aspect("equal")

    g = data.graph()
    if reanalysis and data.reanalysis_edges:
        style = EDGE["paper_project_reanalysed"]
        nx.draw_networkx_edges(
            g, pos, ax=ax, edgelist=data.reanalysis_edges,
            edge_color=style["colour"], width=style["width"], style=style["style"],
            arrows=True, arrowstyle="-",
            connectionstyle=f"arc3,rad={REANALYSIS_ARC}")

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
            nx.draw_networkx_edges(g, pos, ax=ax, edgelist=drawn,
                                   edge_color=EDGE[group]["colour"],
                                   width=EDGE[group]["width"],
                                   style=EDGE[group]["style"])

    def scatter(nodes, kind, zorder):
        if not nodes:
            return
        xy = np.array([pos[n] for n in nodes]).T
        ax.scatter(*xy, s=SIZE.get(kind, SIZE["project"]), c=FILL[kind],
                   marker=MARKER[kind], ec=STROKE[kind][0], lw=STROKE[kind][1],
                   zorder=zorder)

    shown = [c for c in data.cohort if c["BioProject"] in data.used]
    scatter([("R", c["run"]) for c in shown if not field_wheat(c)], "other_run", 2)
    scatter([("R", c["run"]) for c in shown if field_wheat(c)], "field_wheat", 2)
    scatter([("B", b) for b in data.projects if b in data.used], "project", 3)
    scatter([("B", b) for b in data.projects if b not in data.used],
            "project_unused", 3)
    scatter([("P", p) for p in data.papers if p not in data.depositing],
            "no_deposit", 4)
    scatter([("P", p) for p in data.papers if p in data.depositing], "paper", 4)
    scatter([data.search], "search", 5)

    if labels:
        for node in placeable(data):
            ax.annotate(LABEL_FOR(data, node, 0), pos[node],
                        textcoords="offset points", xytext=(0, 7), ha="center",
                        fontsize=MAP_FONT, color=MAP_COLOUR, zorder=6,
                        clip_on=True)

    if LEGEND_ON_FIGURE:
        ax.legend(handles=handles(reanalysis), frameon=False, fontsize=LEGEND_FONT,
                  loc="lower left")
    fig.tight_layout()
    # The two variants are separate files. The reanalysis one is an alternative read of
    # the same layout, not a replacement, so drawing it must not cost the default.
    stem = ("network_solved_map" if labels
            else "network_solved_reanalysis" if reanalysis else "network_solved")
    out = FIG / f"{stem}.png"
    fig.savefig(out, dpi=DPI, facecolor=BACKGROUND)
    plt.close(fig)
    legend_out = FIG / f"{stem}_legend.png"
    write_legend(legend_out, reanalysis)
    print(f"wrote {out}\nwrote {legend_out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--layout", action="store_true",
                    help="recompute the layout and cache it, then draw")
    ap.add_argument("--render", action="store_true",
                    help="draw from the cache; fail rather than recompute")
    ap.add_argument("--report", action="store_true",
                    help="measure the cached layout, draw nothing")
    ap.add_argument("--map", action="store_true", dest="node_map",
                    help="write network_nodes.tsv and a copy of the "
                         "figure with every placeable node named, to work out what to "
                         "put in curation/pinned_nodes.tsv")
    ap.add_argument("--reanalysis", action="store_true",
                    help="also draw, dashed, the edge to data a paper only reanalysed. "
                         "Writes network_solved_reanalysis.png and reuses the same "
                         "layout, since these edges do not steer it")
    args = ap.parse_args()

    FIG.mkdir(exist_ok=True)
    data = Cohort()
    cached = None if args.layout else load_layout()

    if cached is None:
        if args.render or args.report:
            sys.exit(f"no usable layout cached at {CACHE.name}: it is missing, or the "
                     f"geometry or the input tables changed since it was written. "
                     f"Run with --layout.")
        pos, radii, scale, outer_ring, window = compute_layout(data)
        save_layout(pos, radii, scale, outer_ring, window)
        print(f"wrote {CACHE}")
    else:
        pos = cached["positions"]
        radii = cached["radii"]
        outer_ring = cached["outer_ring"]
        window = cached["window"]

    report(data, pos, radii, outer_ring)
    if args.reanalysis:
        print(f"   plus {len(data.reanalysis_edges)} curated reanalysis edges drawn "
              f"dashed. {len(data.mentioned)} further paper-to-deposit links are "
              f"mentions in the text that no one has adjudicated, and are NOT drawn")
    if not args.report:
        render(data, pos, window, reanalysis=args.reanalysis)
    if args.node_map:
        rows = write_node_table(data, pos, cached["scale"] if cached else scale)
        render(data, pos, window, labels=True)
        print(f"wrote {NODE_TABLE}: {len(rows)} placeable nodes "
              f"({sum(1 for r in rows if r['pinned'])} pinned by hand)")
        seeded = seed_pins(rows)
        if seeded:
            print(f"seeded {PINNED} with {seeded} empty rows, one per node")


if __name__ == "__main__":
    main()
