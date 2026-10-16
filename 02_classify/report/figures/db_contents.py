#!/usr/bin/env python3
"""db_contents.py — what db_v4 can see, as a circular taxonomic cladogram.

    python 02_classify/report/figures/db_contents.py

One tip per species the survey found with at least one assembly (3,774). Tips are coloured by
the ROUTE that put the species in scope, and each carries a two-layer radial bar: the pale bar
is every assembly that exists, the solid bar is the subset carrying a gene model. Solid over
pale follows `01_literature`'s `e_design` panel, where solid is the real quantity and pale is
the remainder.

## It is a TAXONOMY, not a phylogeny, and the figure must say so

The topology comes from holobase's NCBI taxonomy: parent/child links only, no branch lengths,
and polytomies wherever NCBI did not resolve a rank. That is the right idiom for a
database-contents figure and the wrong thing to call phylogenetic. A real phylogeny over 3,774
fungal species needs marker-gene alignment and tree inference, which is a separate project and
would not change what this figure is for. Radius encodes taxonomic DEPTH, nothing else.

## The bar is log, and it is not an additive stack

Assemblies per species run from 1 to 1,247 with a median of 2, so a linear bar would be a few
visible species and 3,700 invisible ones.

    pale bar  = log10(total assemblies + 1)
    solid bar = log10(annotated assemblies + 1)

Both are drawn from the same origin, solid over pale, so the visible pale remainder is
log10(total/annotated): the log-fold gap between what exists and what we can use. It is NOT
`log(a) + log(b)`, because that would be meaningless, and the caption says so. A species with
no gene model at all has an entirely pale bar, which is the point: the 2,172 all-pale tips ARE
the miniprot worklist, and on a ring they read as the gaps in the database.

## Colour precedence, because a species is usually in by several routes

Nine route combinations collapse to five colours by taking the strongest evidence:

    host > target > host_associate > plant_associate > congener

`congener` is last and deliberately GREY. Those species are in the database only because they
share a genus with a known pathogen, they exist to give Kraken2's LCA somewhere to park a read,
and `05_detect` must never report a detection on a congener-only taxon the way it reports a
target. Grey says de-emphasised, which is exactly their status.

Written in matplotlib rather than ggtree because the two-layer log bar is not an additive stack
and `geom_fruit` would fight that, and because the rest of the module is Python.
"""

import argparse
import csv
import math
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _layout import SILVER, FIGURES, HOLOBASE, HOSTS   # noqa: E402
from _seeds import read_seed_tsv      # noqa: E402
import _taxonomy                      # noqa: E402

SURVEY = SILVER / "assembly_survey.tsv"
OUT = FIGURES / "db_contents"

# Okabe-Ito, colourblind-safe. congener is grey ON PURPOSE: see the docstring.
# `target` is SPLIT by whether the pathogen is recorded on one of the four cohort cereals.
# Two weights of the same hue rather than two hues, because they are one category at two
# levels: a known pathogen detected on a crop it is already recorded from is an EXPECTED
# finding, while a known pathogen detected on a crop it has never been recorded from is a
# potential new host association, which is what the chapter is looking for.
ROUTE_COLOUR = {
    "host":            "#21531f",
    "target_cereal":   "#8c0404",
    "target_other":    "#be5656",
    "host_associate":  "#22127e",
    "plant_associate": "#545CC5",
    "congener":        "#c0bbbb",
}
PRECEDENCE = ("host", "target_cereal", "target_other", "host_associate",
              "plant_associate", "congener")
ROUTE_LABEL = {
    "host":            "host (plant)",
    "target_cereal":   "pathogen recorded on wheat, rice, maize or sorghum",
    "target_other":    "plant pathogen, not recorded on those crops",
    "host_associate":  "host associate: seen in a phytobiome study",
    "plant_associate": "plant associate: recorded on a plant",
    "congener":        "congener: shares a genus with a pathogen",
}

# GEOMETRY, and the one relation that matters.
#
# The angular pitch per tip is fixed by the tip count: 360 - GAP degrees over however many
# species are drawn, and nothing here can change it. What CAN change is how much of the page
# that pitch is worth, because the figure is scaled to fit its canvas. A tip's bar is
#
#     width_on_page  =  pitch_radians * (canvas_inches / 2) * (R_BAR0 / lim)
#
# so the bars widen only as R_BAR0 approaches `lim`, the outer extent of the drawing. Two
# levers push the same way: a LARGER tree radius moves the bar origin out, and a SHORTER bar
# brings the outer extent in. Both were turned up on 2026-10-16 because the bars still crowded
# at R_TREE 1.0 / BAR_SCALE 0.30, where the ratio was 0.56.
#
# The ratio cannot exceed 1 (that would be bars of zero length), so this is a tradeoff with a
# ceiling rather than a free parameter. draw() prints the width it achieved, so the next
# adjustment is made against a number instead of an impression.
R_TREE = 1.80      # tips of the cladogram sit here
R_INNER = 1.00     # ...and the ROOT sits here, not at the centre. See below.
R_BAR0 = R_TREE + 0.045    # bars start just outside
BAR_SCALE = 0.16   # radial length of one log10 decade

# BRANCH LENGTH is (R_TREE - R_INNER) / 7, because the seven ranks are spread evenly across
# that band. RAISING R_INNER SHORTENS EVERY BRANCH and opens the hole in the middle; dropping
# it to 0 puts the root at the centre and gives the longest branches. The structure worth
# seeing is all in the outer ranks -- family, genus, species -- while the inner ones are a
# handful of long spokes carrying no information, so compressing the band costs nothing and
# buys a tidier tree. At R_TREE 1.80 / R_INNER 1.00 a branch is 0.114 units, against 0.257
# when the tree ran from the centre.

# BRANCH APPEARANCE. Terminal and internal edges are drawn separately and deliberately differ:
# there are as many terminal edges as there are species, 0.22 degrees apart, so at the weight
# of the internal edges they merge into a solid ring and bury the structure. Edit these here
# rather than hunting through draw().
BRANCH_COLOUR = "#8a8a8a"      # internal edges: rank-to-rank, and the arcs joining siblings
BRANCH_WIDTH = 0.30
TIP_BRANCH_COLOUR = "#dcdcdc"  # terminal edges: the genus ring out to the species ring
TIP_BRANCH_WIDTH = 0.12
GAP_DEG = 14.0     # wedge left open at the top: the ring needs a start, and the scale sits in it
LABEL_MIN = 40     # a class needs this many species to earn a label

# Radius is TAXONOMIC RANK, not raw depth in the tree. Raw depth was the first attempt and it
# produced a solid black disc: NCBI nests unevenly, up to 20+ levels, so every terminal branch
# piled up against the rim and 3,744 radial segments 0.1 degrees apart merged into fill. Ranks
# give seven clean concentric rings, and the radius then means something a reader can state.
RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
RANK_LEVEL = {r: i + 1 for i, r in enumerate(RANKS)}


def cereal_recorded(tax) -> set:
    """Internal taxids with ANY recorded host association to one of the four cohort cereals.

    This is a LABEL, not a filter, and the distinction is the whole reason it is allowed here.
    `seeds/README.md` forbids scoping the database by the cohort's own hosts, because the study
    exists to find UNREPORTED host associations and a filter on known ones makes them
    undetectable by construction. Colouring after selection removes nothing: every species is
    still in the index and still classifiable. Same discipline as pathogen status in
    `05_detect`, which is a column rather than a gate.

    **Absence of a record is not evidence that a fungus cannot infect a cereal.** USDA records
    what has been observed and written down, so the unmarked 86% is precisely the population in
    which an unreported association could be found. The caption says so.

    A record against any descendant of a cereal genus counts (cultivars, subspecies), and so
    does a record held against a sub-taxon of the fungus.
    """
    crops = {}
    for row in read_seed_tsv(HOSTS, require=("ncbi_taxid", "species", "role")):
        if (row.get("role") or "").strip() != "cohort":
            continue          # the two outgroups are not cohort crops
        tid = tax.resolve(row["ncbi_taxid"])
        genus = tax.ancestor_at(tid, "genus") if tid else 0
        for d in tax.descendants(genus or tid):
            crops[d] = row["species"]
    con = sqlite3.connect(HOLOBASE)
    direct = {f for f, h in con.execute(
        "SELECT fungus_id, host_id FROM host_association WHERE host_id IS NOT NULL")
        if h in crops}
    con.close()
    # roll a record held against a strain or forma specialis up to the species node
    out = set(direct)
    for f in direct:
        sp = tax.species_of(f)
        if sp:
            out.add(sp)
    return out


def primary_route(routes: str, on_cereal: bool) -> str:
    """One colour per tip: the strongest evidence, with `target` split by cereal record.

    Splitting inside `target` rather than drawing a separate overlay ring keeps the figure to
    one mark per tip. The ring was the first attempt and it broke that, which is why it went.

    The split loses one thing worth knowing: 92 species that are recorded on a cohort cereal
    arrive by a route other than `target` (88 plant_associate, 2 host_associate, 2 congener),
    so they are fungi recorded on our crops without being known pathogens of anything. They
    keep their own route colour. That is a real category and it stays in assembly_survey.tsv.
    """
    present = routes.split(";")
    if "target" in present:
        return "target_cereal" if on_cereal else "target_other"
    for route in PRECEDENCE:
        if route in present:
            return route
    return "congener"


def load_species(tax, on_cereal: set) -> tuple:
    """Survey rows joined to the taxonomy. Only species that can ENTER db_v4 become tips.

    The first version drew all 3,774 species with an assembly, so the 2,172 with no gene model
    appeared as all-pale tips and the miniprot gap was a visible feature of the ring. Leon's
    call after seeing it: too many tips to be informative. At 3,744 tips the angular pitch is
    0.09 degrees, which is NARROWER THAN THE BARS THEMSELVES near the origin, so the figure
    was not merely crowded, it was wrong -- see draw() for the arithmetic. The excluded count
    stays in the subtitle and the caption, where it is still a number a reader needs.
    """
    out, unplaced, no_model = [], 0, 0
    for r in csv.DictReader(open(SURVEY), delimiter="\t"):
        if not int(r["n_annotated"] or 0):
            no_model += 1
            continue
        tid = tax.by_ncbi.get(int(r["ncbi_taxid"]))
        if tid is None:
            unplaced += 1
            continue
        out.append({
            "tid": tid,
            "name": r["name"],
            "total": int(r["n_assemblies"] or 0),
            "annotated": int(r["n_annotated"] or 0),
            "route": primary_route(r["routes"] or "congener", tid in on_cereal),
        })
    if unplaced:
        print(f"  {unplaced} species are not in holobase's tree and cannot be drawn "
              f"(newer than the snapshot); they stay in the table")
    print(f"  {no_model:,} species have assemblies but no gene model and are NOT tips; "
          f"they are the miniprot worklist")
    return out, no_model, unplaced


def build_tree(tax, species: list) -> tuple:
    """Lineages -> a collapsed cladogram. Returns (children, root, parent)."""
    tips = {s["tid"] for s in species}
    paths = {}
    for tid in tips:
        chain, cur, seen = [], tid, set()
        while cur is not None and cur not in seen:
            seen.add(cur)
            chain.append(cur)
            nxt = tax.parent.get(cur)
            if nxt is None or nxt == cur:
                break
            cur = nxt
        paths[tid] = list(reversed(chain))

    children = defaultdict(set)
    for chain in paths.values():
        for a, b in zip(chain, chain[1:]):
            children[a].add(b)
    root = next(iter(paths.values()))[0]

    # Keep ONLY the standard ranks. Everything else (clade, subkingdom, no rank, subphylum,
    # section, ...) is skipped and its children re-attach to the nearest ranked ancestor. This
    # is what turns the tree from a 20-deep tangle into seven concentric rings, and it also
    # removes the single-child chains that an unranked hierarchy is mostly made of.
    def keep(node):
        return node in tips or tax.rank.get(node) in RANK_LEVEL

    def ranked_kids(node):
        out = []
        for kid in children.get(node, ()):
            if keep(kid):
                out.append(kid)
            else:
                out.extend(ranked_kids(kid))
        return out

    sys.setrecursionlimit(50000)
    rebuilt = {}
    stack = [root]
    while stack:
        node = stack.pop()
        if node in rebuilt:
            continue
        kids = [] if node in tips else ranked_kids(node)
        rebuilt[node] = set(kids)
        stack.extend(kids)
    children = defaultdict(set, rebuilt)
    while root not in tips and len(children.get(root, ())) == 1:
        root = next(iter(children[root]))

    # Prune to what is actually reachable from the collapsed root. `children` was accumulated
    # over the raw lineages, so it still holds every node the collapse skipped over -- including
    # the nodes ABOVE the new root. Leaving them in means draw() meets a node with no radius.
    reachable, stack = {}, [root]
    while stack:
        node = stack.pop()
        if node in reachable:
            continue
        kids = set(children.get(node, ()))
        reachable[node] = kids
        stack.extend(kids)

    parent = {}
    for node, kids in reachable.items():
        for kid in kids:
            parent[kid] = node
    return reachable, root, parent


def layout(tax, children, root, tips: set) -> tuple:
    """Angles for every node, and depth for the radius. Tips in traversal order stay contiguous."""
    order = []

    def walk(node, depth, depths):
        depths[node] = depth
        kids = sorted(children.get(node, ()), key=lambda k: tax.label(k))
        if not kids or node in tips:
            order.append(node)
            return
        for kid in kids:
            walk(kid, depth + 1, depths)

    depths = {}
    walk(root, 0, depths)
    span = math.radians(360.0 - GAP_DEG)
    start = math.radians(90.0 + GAP_DEG / 2.0)
    angle = {t: start + span * (i + 0.5) / len(order) for i, t in enumerate(order)}

    def set_angle(node):
        if node in angle:
            return angle[node]
        kids = children.get(node, ())
        vals = [set_angle(k) for k in kids]
        angle[node] = (min(vals) + max(vals)) / 2.0
        return angle[node]

    set_angle(root)
    # Radius by RANK, spread over the band [R_INNER, R_TREE] rather than from the centre, so
    # branch length is set by R_INNER. A node at an unlisted rank (only the root can be, after
    # the prune) sits at R_INNER; a tip always sits on the rim even if NCBI ranks it below
    # species.
    top = len(RANKS)
    radius = {}
    for node in depths:
        lvl = RANK_LEVEL.get(tax.rank.get(node), 0)
        radius[node] = R_INNER + (R_TREE - R_INNER) * lvl / top
    for t in order:
        radius[t] = R_TREE
    return angle, radius, order


def polar(a, r):
    return r * math.cos(a), r * math.sin(a)


def draw(tax, species, children, root, angle, radius, order, no_model, assembly: bool):
    tip_of = {s["tid"]: s for s in species}
    # Self-contained mode reserves a band at the bottom for the legend, because the circle
    # fills its axes and an overlaid legend collided with the class labels. Assembly mode has
    # no text at all, so the marks can use the whole square canvas.
    fig, ax = plt.subplots(figsize=(15.0, 15.0) if assembly else (15.0, 15.6))
    if assembly:
        fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    else:
        fig.subplots_adjust(left=0.02, right=0.98, bottom=0.10, top=0.930)

    # ── the cladogram: a radial segment per edge, an arc per internal node ────────────────
    # Internal edges and TERMINAL edges are drawn separately. There are 3,744 terminal edges at
    # 0.09 degrees apart, so at the same weight as the rest they merge into a solid ring and
    # bury the structure the figure is for. Lighter and thinner, they read as the species ring.
    segs, tip_segs = [], []
    for node, kids in children.items():
        if not kids:
            continue
        r = radius[node]
        kid_angles = [angle[k] for k in kids]
        lo, hi = min(kid_angles), max(kid_angles)
        steps = max(2, int(math.degrees(hi - lo) * 2))
        segs.append([polar(lo + (hi - lo) * i / steps, r) for i in range(steps + 1)])
        for k in kids:
            (tip_segs if k in tip_of else segs).append(
                [polar(angle[k], r), polar(angle[k], radius[k])])
    ax.add_collection(LineCollection(tip_segs, colors=TIP_BRANCH_COLOUR,
                                     linewidths=TIP_BRANCH_WIDTH, zorder=1))
    ax.add_collection(LineCollection(segs, colors=BRANCH_COLOUR,
                                     linewidths=BRANCH_WIDTH, zorder=2))

    # ── the bars ──────────────────────────────────────────────────────────────────────────
    # WEDGES, not line segments. A line has a width in POINTS, fixed in screen space, while
    # the slot available to a tip has a width in DATA space that grows with radius. At 3,744
    # tips the pitch at the bar origin was 0.30 pt against a 0.45 pt line, so every bar was 34%
    # wider than its slot; because the pale layer is drawn before the solid layer, solid
    # systematically overdrew its neighbours' pale and inflated the one ratio this figure
    # exists to show. A wedge is defined in data space and tiles its slot exactly, so the
    # encoding stays faithful at any tip count.
    half = math.radians(360.0 - GAP_DEG) / len(order) / 2.0
    verts, cols = [], []
    for tid in order:
        s = tip_of[tid]
        a = angle[tid]
        col = ROUTE_COLOUR[s["route"]]
        # The bar is the USABLE assemblies and nothing else. An earlier version drew a pale
        # bar for every assembly that exists with the annotated subset solid over it, but once
        # the tips were restricted to species that carry a gene model, the pale layer stopped
        # answering the figure's question ("what can db_v4 see") and started answering a
        # different one ("how much of NCBI do we use"). That second number is real and is kept
        # in assembly_survey.tsv; it does not belong on this ring. Leon's call 2026-10-16.
        r1 = R_BAR0 + math.log10(s["annotated"] + 1) * BAR_SCALE
        verts.append([polar(a - half, R_BAR0), polar(a + half, R_BAR0),
                      polar(a + half, r1), polar(a - half, r1)])
        cols.append(col)
    ax.add_collection(PolyCollection(verts, facecolors=cols, edgecolors="none",
                                     linewidths=0, zorder=3))

    # ── log gridlines behind the bars, at 1 / 10 / 100 / 1000 assemblies ──────────────────
    top = max(s["annotated"] for s in species)
    for n in [n for n in (1, 10, 100, 1000) if n <= top * 1.5]:   # labelled up the GAP
        r = R_BAR0 + math.log10(n + 1) * BAR_SCALE
        circle = plt.Circle((0, 0), r, fill=False, ec="#b0b0b0", lw=0.5,
                            ls=(0, (2, 3)), zorder=5)
        ax.add_patch(circle)
        # The log scale numerals stay in BOTH modes. Assembly mode strips labels that get
        # re-set in Illustrator, but an axis with no numbers is not a stripped label, it is an
        # unreadable axis -- and the house assembly style keeps axis numerals for that reason.
        ax.text(*polar(math.radians(90.0), r + 0.012), f"{n:,}",
                fontsize=14 if assembly else 10, color="#333333" if assembly else "#555555",
                ha="center", va="bottom", zorder=6)

    # ── class labels outside the ring ─────────────────────────────────────────────────────
    label_widths = []
    by_class = defaultdict(list)
    for tid in order:
        cls = tax.ancestor_at(tid, "class")
        by_class[tax.label(cls) if cls else None].append(angle[tid])
    # Sit the class ring just outside the LONGEST BAR, not at a hardcoded decade. Pinning it
    # to 1,000 when the deepest species holds 294 annotated assemblies left a dead annulus
    # across a third of the radius.
    r_label = R_BAR0 + math.log10(top + 1) * BAR_SCALE + 0.05
    for name, angs in by_class.items():
        if name is None or len(angs) < LABEL_MIN:
            continue
        label_widths.append((name, len(angs)))
        mid = (min(angs) + max(angs)) / 2.0
        deg = math.degrees(mid) % 360
        flip = 90 < deg < 270
        if not assembly:
            ax.text(*polar(mid, r_label), f"{name}  ({len(angs)})",
                    rotation=deg + 180 if flip else deg, rotation_mode="anchor",
                    ha="right" if flip else "left", va="center", fontsize=10.5, zorder=7)
        lo, hi = min(angs), max(angs)
        steps = max(2, int(math.degrees(hi - lo) * 2))
        ax.add_collection(LineCollection(
            [[polar(lo + (hi - lo) * i / steps, r_label - 0.035) for i in range(steps + 1)]],
            colors="#333333", linewidths=1.6, zorder=6))

    n_ann = sum(1 for s in species if s["annotated"])
    # ASSEMBLY MODE strips every piece of text and the legend, leaving the marks for Leon to
    # label in Illustrator. The clade names, counts and angular positions go to the sidecar
    # TSV instead, so nothing has to be eyeballed off the image.
    if not assembly:
        fig.suptitle("What db_v4 can see", fontsize=20, y=0.995)
        fig.text(0.5, 0.967,
                 f"{len(species):,} fungal species carrying a gene model, so usable in a "
                 f"CDS-only index.  Bar: annotated assemblies, log scale",
                 fontsize=13, ha="center", color="#333333")
        fig.text(0.5, 0.950,
                 f"a further {no_model:,} have assemblies but no gene model and are not "
                 f"shown: they are the gene-prediction worklist",
                 fontsize=11.5, ha="center", color="#666666")
        handles = [Patch(facecolor=ROUTE_COLOUR[r], edgecolor="none", label=ROUTE_LABEL[r])
                   for r in PRECEDENCE]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=3,
                   frameon=False, fontsize=12, handlelength=1.5, labelspacing=0.5,
                   columnspacing=2.5)

    # Class labels are drawn OUTSIDE r_label and matplotlib does not grow the limits for text,
    # so the pad has to cover the longest label. "Eurotiomycetes (642)" sits at 180 degrees and
    # was clipped by the frame at +0.30.
    # ~0.030 data units per character at fontsize 10.5 on this canvas. Measured by eye against
    # the longest label rather than computed: matplotlib can only report text extents after a
    # draw, and a renderer round-trip is not worth it for one padding constant.
    longest = max((len(f"{n}  ({c})") for n, c in label_widths), default=20)
    lim = r_label + (0.02 if assembly else 0.030 * longest)

    # What the geometry actually bought, in the units a printer cares about.
    pitch = math.radians(360.0 - GAP_DEG) / len(order)
    axes_in = fig.get_size_inches()[0] * (0.98 if assembly else 0.96)
    width_in = pitch * (axes_in / 2.0) * (R_BAR0 / lim)
    print(f"  bar width at the origin: {width_in * 72:.2f} pt "
          f"({width_in * 300:.1f} px at 300 dpi), "
          f"R_BAR0/lim = {R_BAR0 / lim:.2f}, pitch = {math.degrees(pitch):.3f} deg")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")
    return fig


def clade_table(tax, species, angle, order, path: Path) -> list:
    """Every clade's name, size and angular position, in ring order, at class and order rank.

    This is what makes assembly mode usable: with the labels stripped there is nothing on the
    image saying which wedge is which, so the positions have to leave the script as DATA
    rather than be measured off the PNG by hand. Angles are degrees anticlockwise from east,
    matching matplotlib, with the ring starting and ending at the gap near 90 degrees.
    """
    rows = []
    for rank in ("class", "order"):
        groups = defaultdict(list)
        for tid in order:
            anc = tax.ancestor_at(tid, rank)
            # RAW degrees, no modulo. The ring runs from ~97 to ~443 degrees, so a clade
            # straddling 360 gets min 0.2 / max 360.0 if each tip is wrapped first, which
            # reported Agaricales as spanning the entire circle. Take min and max on the
            # monotonic values, then wrap once for display.
            groups[tax.label(anc) if anc else "(unplaced)"].append(math.degrees(angle[tid]))
        for name, degs in groups.items():
            lo, hi = min(degs), max(degs)
            rows.append({
                "rank": rank, "clade": name, "n_species": len(degs),
                "start_deg": round(lo % 360, 2), "end_deg": round(hi % 360, 2),
                "mid_deg": round(((lo + hi) / 2.0) % 360, 2),
                "span_deg": round(hi - lo, 2),
                "wraps_zero": "yes" if (lo % 360) > (hi % 360) else "",
            })
    rows.sort(key=lambda r: (r["rank"], r["start_deg"]))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def caption(species: list, no_model: int, n_cereal: int, n_cereal_other: int,
            unplaced: int, clades: list, assembly: bool) -> str:
    """The full record of the figure, generated from the data it was drawn from.

    Structured like `01_literature/report/figures/map_crop_origins_caption.txt`: a full
    caption, an explicit WHAT IS NOT SHOWN section, a short caption, and the sources. Nothing
    here is hand-written, because a hand-written caption goes stale silently -- this one
    already caught itself describing pale bars three renders after they were removed.
    """
    n = len(species)
    asm = sum(s["total"] for s in species)
    ann = sum(s["annotated"] for s in species)
    by_route = defaultdict(int)
    for s in species:
        by_route[s["route"]] += 1
    deepest = max(species, key=lambda s: s["annotated"])
    classes = [c for c in clades if c["rank"] == "class"]
    # TRUE ring order: distance travelled from the ring's start, not raw degrees. Sorting on
    # start_deg put Agaricomycetes first because it begins at 350 and wraps through 0, when it
    # actually sits at the END of the ring. The ring starts at 90 + GAP/2.
    ring_start = 90.0 + GAP_DEG / 2.0
    classes.sort(key=lambda c: (c["start_deg"] - ring_start) % 360)
    named = [c for c in classes if c["n_species"] >= LABEL_MIN]
    rest = [c for c in classes if c["n_species"] < LABEL_MIN]

    def row(c):
        return (f"  {c['clade']:<22} {c['n_species']:>5}   "
                f"{c['start_deg']:>6.1f} to {c['end_deg']:>6.1f}   {c['mid_deg']:>6.1f}   "
                f"{c['span_deg']:>6.1f}" + ("   wraps 0" if c["wraps_zero"] else ""))

    return f"""Reference scope of db_v4
Caption for report/figures/db_contents.{{png,pdf}}

GENERATED by report/figures/db_contents.py -- do not hand-edit, the numbers will drift.

Rendered in {'ASSEMBLY' if assembly else 'SELF-CONTAINED'} mode.{
' All text except the log-axis numerals is stripped from the image; the class'
' names, counts and angular positions below are what you place in Illustrator.'
if assembly else ' Title, legend, class labels and scale are drawn on the image.'}


FULL CAPTION

Reference scope of db_v4, the Kraken2 database this chapter classifies against. Circular
cladogram of the {n:,} fungal species that carry a gene model and can therefore enter a
CDS-only index.

TOPOLOGY IS TAXONOMIC, NOT PHYLOGENETIC. The tree is the NCBI taxonomy held in holobase:
parent/child links only, with polytomies wherever NCBI did not resolve a rank. Branch lengths
are not evolutionary distance, and RADIUS ENCODES TAXONOMIC RANK -- kingdom, phylum, class,
order, family, genus, species, as seven concentric rings -- rather than any measure of
divergence. A real phylogeny over {n:,} species would need marker-gene alignment and tree
inference; it is a separate exercise and would not change what this figure is for.

THE BAR is the number of ANNOTATED assemblies the species holds, on a log10 scale, {ann:,}
across the ring. The deepest is {deepest['name']} at {deepest['annotated']:,}. Bars are drawn
as WEDGES that tile their angular slot exactly rather than as fixed-width lines, so no tip
overdraws its neighbour at any radius; with a line of fixed point-width the bars were 34%
wider than their slot near the origin and systematically overdrew each other.

TIP COLOUR is the route by which the species entered scope, taking the strongest evidence
where several apply. Five routes, unioned, with no floors or ceilings:

  {ROUTE_LABEL['target_cereal']:<52} {by_route['target_cereal']:>5}
  {ROUTE_LABEL['target_other']:<52} {by_route['target_other']:>5}
  {ROUTE_LABEL['host_associate']:<52} {by_route['host_associate']:>5}
  {ROUTE_LABEL['plant_associate']:<52} {by_route['plant_associate']:>5}
  {ROUTE_LABEL['congener']:<52} {by_route['congener']:>5}
  {ROUTE_LABEL['host']:<52} {by_route['host']:>5}

Sources per route: targets are FUNGuild's Plant Pathogen guild plus PHI-base pathogenic
associations WHOSE HOST IS A PLANT; congeners are every genus holding a target; plant
associates are USDA Fungus-Host `recorded_on`; host associates come from the phytobiome
abundance literature for the four crops; hosts are the four cohort cereals plus two
outgroups.

CONGENERS ARE GREY BECAUSE THEY ARE NOT REPORTABLE. They are in the index to give Kraken2's
LCA somewhere to park a read, and a detection on a congener-only taxon must not be reported
the way a detection on a target is. The deepest species in the whole database is a congener.

KNOWN PLANT PATHOGENS ARE SPLIT by whether they are recorded on one of the four cohort
cereals ({n_cereal:,}, {n_cereal / n:.0%} of tips) or not ({by_route['target_other']:,}). A
pathogen detected on a crop it is already recorded from is an EXPECTED finding; the same
pathogen on a crop it has never been recorded from is a potential new host association, which
is what this chapter looks for. The split is a LABEL AND NOT A FILTER: scoping the database
by the cohort's own hosts is forbidden, because the study exists to find UNREPORTED
associations and such a filter would make them undetectable by construction. ABSENCE OF A
RECORD IS NOT EVIDENCE that a fungus cannot infect a cereal -- the source holds only what has
been observed and written down -- so the species without one are precisely the population in
which an unreported association could be found.


CLASSES, IN RING ORDER

Degrees anticlockwise from east, matching the drawing. The ring starts and ends at the gap
centred on 90 degrees, so it runs 97 -> 360 -> 83. {len(named)} classes hold at least
{LABEL_MIN} species and are the ones worth labelling.

  class                species        start to end      mid     span
{chr(10).join(row(c) for c in named)}

  below {LABEL_MIN} species: {len(rest)} classes, {sum(c['n_species'] for c in rest):,} species
{chr(10).join(row(c) for c in rest)}

Order-rank spans are in db_contents_clades.tsv alongside these.


WHAT IS NOT SHOWN

{no_model:,} fungal species have at least one assembly but NO GENE MODEL anywhere, and are
absent from the figure entirely. They are not a rejection: a CDS-only index cannot hold them
until gene prediction has run, and they are the worklist in data/silver/miniprot_worklist.tsv.
Including them was the first version of this figure; at {n + no_model:,} tips the angular
pitch fell below the width of the bars themselves and the encoding stopped being faithful.

{asm - ann:,} further assemblies belonging to the species that ARE shown carry no gene model
and are not drawn. The figure therefore answers "what can db_v4 see", not "how much of NCBI
do we use"; the second number is in data/silver/assembly_survey.tsv.

{n_cereal_other:,} species are recorded on a cohort cereal but entered by a route other than
`target`, so they are fungi recorded on our crops WITHOUT being known pathogens of anything.
They keep their own route colour and are not counted in the {n_cereal:,} above.

{unplaced:,} surveyed species are newer than holobase's taxonomy snapshot, cannot be placed
on the tree, and are absent from the figure while remaining in the survey table.

Oomycetes, bacteria, viruses and mammals are all absent by decision rather than oversight.
Fungi-only was settled 2026-10-13; off-target sinks were measured and rejected 2026-10-15,
because bacterial reads go UNCLASSIFIED rather than into fungi (at >20% bacterial load the
unclassified fraction rises to 52% while the fungal call falls), human signal is negligible
in these libraries (16 of 7,486 cohort runs above 1%), and viral genomes hold too few
minimizers for a fractional coverage floor to behave.

Only six plant hosts are present: the four cohort cereals plus Hordeum and Arabidopsis as
monocot and dicot outgroups. db_v3 held 38, and the outgroups are kept because the host call
RANKS rather than thresholds, so a mislabelled non-cereal library stays detectable.


SHORT CAPTION

Reference scope of db_v4: {n:,} fungal species carrying a gene model, on a circular NCBI
taxonomic cladogram where radius is taxonomic rank, not divergence. Bar length is the
species' annotated assembly count on a log10 scale ({ann:,} total). Colour is the route into
scope: known plant pathogens recorded on the four cohort cereals ({by_route['target_cereal']:,})
or not ({by_route['target_other']:,}), host associates ({by_route['host_associate']:,}), plant
associates ({by_route['plant_associate']:,}), congeners present only as classifier competitors
and not reportable ({by_route['congener']:,}), and plant hosts ({by_route['host']:,}). A
further {no_model:,} species have assemblies but no gene model and are not shown.


SOURCES

  figure   02_classify/report/figures/db_contents.py
  clades   02_classify/report/figures/db_contents_clades.tsv
  data     02_classify/data/silver/assembly_survey.tsv
  routes   02_classify/01_refs/routes.py   (membership, five routes)
  survey   02_classify/01_refs/survey.py   (NCBI assembly + annotation counts)
  seeds    02_classify/seeds/  (hosts.tsv, residents.tsv, fungi_usda.tsv, exclusions.tsv)
  taxonomy _refdata/holobase.db via 02_classify/_taxonomy.py
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", default="assembly", choices=("assembly", "self-contained"),
                    help="assembly (default): no text at all, for labelling in Illustrator, "
                         "with every clade's position in db_contents_clades.tsv. "
                         "self-contained: title, legend, class labels and scale kept")
    args = ap.parse_args()
    assembly = args.mode == "assembly"

    tax = _taxonomy.load()
    on_cereal = cereal_recorded(tax)
    species, no_model, unplaced = load_species(tax, on_cereal)
    # Recorded on a cohort cereal but NOT a known pathogen, so not in the target_cereal count.
    n_cereal_other = sum(1 for s in species
                         if s["tid"] in on_cereal and s["route"] != "target_cereal")
    print(f"  {len(species):,} species to draw")
    children, root, parent = build_tree(tax, species)
    tips = {s["tid"] for s in species}
    angle, radius, order = layout(tax, children, root, tips)
    print(f"  cladogram: {len(order):,} tips, {len(children):,} internal nodes, "
          f"root {tax.label(root)}")
    fig = draw(tax, species, children, root, angle, radius, order, no_model, assembly)
    n_cereal = sum(1 for s in species if s["route"] == "target_cereal")
    print(f"  {n_cereal:,} of {len(species):,} drawn species are pathogens recorded on a "
          f"cohort cereal")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    clades = clade_table(tax, species, angle, order,
                         Path(f"{OUT}_clades.tsv"))
    for ext in ("png", "pdf"):
        # No bbox_inches="tight": the reserved legend band is deliberate and tight cropping
        # would undo it.
        fig.savefig(f"{OUT}.{ext}", dpi=300, facecolor="white")
    Path(f"{OUT}_caption.txt").write_text(
        caption(species, no_model, n_cereal, n_cereal_other, unplaced, clades,
                assembly))
    print(f"  -> {OUT}.png / .pdf / _caption.txt / _clades.tsv   (mode: {args.mode})")


if __name__ == "__main__":
    main()
