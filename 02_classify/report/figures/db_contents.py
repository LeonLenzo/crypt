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
    "host":            "#cc79a7",
    "target_cereal":   "#8c2d04",
    "target_other":    "#d55e00",
    "host_associate":  "#009e73",
    "plant_associate": "#0072b2",
    "congener":        "#8c8c8c",
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

R_TREE = 1.00      # tips of the cladogram sit here
R_BAR0 = 1.045     # bars start just outside
BAR_SCALE = 0.30   # radial length of one log10 decade
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
    return out, no_model


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
    # Radius by RANK. A node at an unlisted rank (only the root can be, after the prune) sits at
    # the centre; a tip always sits on the rim even if NCBI ranks it below species.
    top = len(RANKS)
    radius = {}
    for node in depths:
        lvl = RANK_LEVEL.get(tax.rank.get(node), 0)
        radius[node] = R_TREE * lvl / top
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
    fig, ax = plt.subplots(figsize=(13.5, 13.5) if assembly else (13.5, 14.1))
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
    ax.add_collection(LineCollection(tip_segs, colors="#dcdcdc", linewidths=0.12, zorder=1))
    ax.add_collection(LineCollection(segs, colors="#8a8a8a", linewidths=0.3, zorder=2))

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


def caption(species: list, no_model: int, n_cereal: int, clades: list) -> str:
    """Generated from the data. A hand-written caption goes stale silently; this cannot.

    It has already proved the point: an earlier edit dropped the pale half of the bars from the
    figure and the hand-edited sentence describing them survived three renders.
    """
    n = len(species)
    asm = sum(s["total"] for s in species)
    ann = sum(s["annotated"] for s in species)
    by_route = defaultdict(int)
    for s in species:
        by_route[s["route"]] += 1
    routes = ", ".join(f"{ROUTE_LABEL[r].split(':')[0]} {by_route[r]:,}"
                       for r in PRECEDENCE if by_route[r])
    return (
        f"Reference scope of db_v4. Circular cladogram of the {n:,} fungal species that carry "
        f"a gene model and can therefore enter a CDS-only index, drawn on the NCBI taxonomy "
        f"held in holobase: topology is taxonomic, branch lengths are not evolutionary "
        f"distance, and radius encodes taxonomic RANK (kingdom, phylum, class, order, family, "
        f"genus, species) rather than any measure of divergence. "
        f"Tip colour gives the route that brought the species into scope, taking the strongest "
        f"evidence where several apply ({routes}); congeners are grey because they are present "
        f"to give the classifier a competitor and are not reportable as detections. "
        f"The radial bar is the number of ANNOTATED assemblies the species holds on a log10 "
        f"scale, {ann:,} in total; the {asm - ann:,} further assemblies these species hold "
        f"without a gene model are not drawn, because a CDS-only index cannot use them. Bars "
        f"are wedges that tile their angular slot exactly rather than fixed-width lines, so no "
        f"tip overdraws its neighbour at any radius. "
        f"Known plant pathogens are split by whether they are recorded on one of the four "
        f"cohort cereals ({n_cereal:,} of them, {n_cereal / n:.0%} of all tips): a pathogen "
        f"detected on a crop it is already recorded from is an expected finding, while one "
        f"detected on a crop it has never been recorded from is a potential new host "
        f"association. The split is a LABEL and not a filter: scoping the database by the cohort's own hosts is forbidden, because the "
        f"study exists to find UNREPORTED associations and such a filter would make them "
        f"undetectable by construction. Absence of a record is NOT evidence that a fungus "
        f"cannot infect a cereal, because the source holds only what has been observed and "
        f"written down; the species without one are precisely the population in which an "
        f"unreported association could be found. "
        f"A further {no_model:,} species have at least one assembly but no gene model anywhere "
        f"and are not shown, because a CDS-only index cannot hold them until gene prediction "
        f"has run; they are the worklist, not a rejection. Classes with at least {LABEL_MIN} "
        f"species are labelled.\n\n"
        f"ASSEMBLY MODE: all text is stripped from the image. Counts to set in Illustrator, "
        f"by route: "
        + ", ".join(f"{ROUTE_LABEL[r]} {by_route[r]:,}" for r in PRECEDENCE if by_route[r])
        + ".\n\nClades in ring order (degrees anticlockwise from east; the ring starts and "
        f"ends at the gap near 90 deg). Full table with class rank as well in "
        f"db_contents_clades.tsv:\n"
        + "\n".join(
            f"  {c['clade']:<26} {c['n_species']:>5} species   "
            f"{c['start_deg']:>6.1f} to {c['end_deg']:>6.1f} deg  (mid {c['mid_deg']:>6.1f}, "
            f"span {c['span_deg']:>5.1f}{', wraps 0' if c['wraps_zero'] else ''})"
            for c in clades if c["rank"] == "order" and c["n_species"] >= 20)
    )


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
    species, no_model = load_species(tax, on_cereal)
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
        caption(species, no_model, n_cereal, clades) + "\n")
    print(f"  -> {OUT}.png / .pdf / _caption.txt / _clades.tsv   (mode: {args.mode})")


if __name__ == "__main__":
    main()
