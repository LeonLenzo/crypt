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

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _layout import SILVER, FIGURES   # noqa: E402
import _taxonomy                      # noqa: E402

SURVEY = SILVER / "assembly_survey.tsv"
OUT = FIGURES / "db_contents"

# Okabe-Ito, colourblind-safe. congener is grey ON PURPOSE: see the docstring.
ROUTE_COLOUR = {
    "host":            "#cc79a7",
    "target":          "#d55e00",
    "host_associate":  "#009e73",
    "plant_associate": "#0072b2",
    "congener":        "#8c8c8c",
}
PRECEDENCE = ("host", "target", "host_associate", "plant_associate", "congener")
ROUTE_LABEL = {
    "host":            "host (plant)",
    "target":          "target: known plant pathogen",
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


def pale(hex_colour: str, amount: float = 0.74) -> tuple:
    """Blend toward white. The pale/solid pair must stay the same hue to read as one bar."""
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return tuple(c + (1.0 - c) * amount for c in (r, g, b))


def primary_route(routes: str) -> str:
    for route in PRECEDENCE:
        if route in routes.split(";"):
            return route
    return "congener"


def load_species(tax) -> list:
    """Survey rows joined to the taxonomy, dropping anything the tree cannot place."""
    out, unplaced = [], 0
    for r in csv.DictReader(open(SURVEY), delimiter="\t"):
        tid = tax.by_ncbi.get(int(r["ncbi_taxid"]))
        if tid is None:
            unplaced += 1
            continue
        out.append({
            "tid": tid,
            "name": r["name"],
            "total": int(r["n_assemblies"] or 0),
            "annotated": int(r["n_annotated"] or 0),
            "route": primary_route(r["routes"] or "congener"),
        })
    if unplaced:
        print(f"  {unplaced} surveyed species are not in holobase's tree and cannot be drawn "
              f"(newer than the snapshot); they stay in the table")
    return out


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


def draw(tax, species, children, root, angle, radius, order):
    tip_of = {s["tid"]: s for s in species}
    fig, ax = plt.subplots(figsize=(13.5, 14.1))
    # The circle fills its axes, so the legend can only be kept off the class labels by
    # reserving a band for it. Overlaying it collided with Pichiomycetes and Saccharomycetes.
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.10, top=0.945)

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
    half = math.radians(360.0 - GAP_DEG) / len(order) * 0.46
    pale_segs, pale_cols, solid_segs, solid_cols = [], [], [], []
    for tid in order:
        s = tip_of[tid]
        a = angle[tid]
        col = ROUTE_COLOUR[s["route"]]
        r_tot = R_BAR0 + math.log10(s["total"] + 1) * BAR_SCALE
        pale_segs.append([polar(a, R_BAR0), polar(a, r_tot)])
        pale_cols.append(pale(col))
        if s["annotated"]:
            r_ann = R_BAR0 + math.log10(s["annotated"] + 1) * BAR_SCALE
            solid_segs.append([polar(a, R_BAR0), polar(a, r_ann)])
            solid_cols.append(col)
    lw = max(0.45, 1800.0 / len(order) * 0.55)
    ax.add_collection(LineCollection(pale_segs, colors=pale_cols, linewidths=lw,
                                     zorder=3, capstyle="butt"))
    ax.add_collection(LineCollection(solid_segs, colors=solid_cols, linewidths=lw,
                                     zorder=4, capstyle="butt"))

    # ── log gridlines behind the bars, at 1 / 10 / 100 / 1000 assemblies ──────────────────
    for n in (1, 10, 100, 1000):   # labelled up the GAP, so no tick sits over a bar
        r = R_BAR0 + math.log10(n + 1) * BAR_SCALE
        circle = plt.Circle((0, 0), r, fill=False, ec="#b0b0b0", lw=0.5,
                            ls=(0, (2, 3)), zorder=5)
        ax.add_patch(circle)
        ax.text(*polar(math.radians(90.0), r + 0.012), f"{n:,}", fontsize=10, color="#555555",
                ha="center", va="bottom", zorder=6)

    # ── class labels outside the ring ─────────────────────────────────────────────────────
    label_widths = []
    by_class = defaultdict(list)
    for tid in order:
        cls = tax.ancestor_at(tid, "class")
        by_class[tax.label(cls) if cls else None].append(angle[tid])
    r_label = R_BAR0 + math.log10(1000 + 1) * BAR_SCALE + 0.06
    for name, angs in by_class.items():
        if name is None or len(angs) < LABEL_MIN:
            continue
        label_widths.append((name, len(angs)))
        mid = (min(angs) + max(angs)) / 2.0
        deg = math.degrees(mid) % 360
        flip = 90 < deg < 270
        ax.text(*polar(mid, r_label), f"{name}  ({len(angs)})",
                rotation=deg + 180 if flip else deg, rotation_mode="anchor",
                ha="right" if flip else "left", va="center", fontsize=10.5, zorder=7)
        lo, hi = min(angs), max(angs)
        steps = max(2, int(math.degrees(hi - lo) * 2))
        ax.add_collection(LineCollection(
            [[polar(lo + (hi - lo) * i / steps, r_label - 0.035) for i in range(steps + 1)]],
            colors="#333333", linewidths=1.6, zorder=6))

    n_ann = sum(1 for s in species if s["annotated"])
    fig.suptitle("What db_v4 can see", fontsize=20, y=0.988)
    fig.text(0.5, 0.962,
             f"{len(species):,} fungal species with an assembly; "
             f"{n_ann:,} carry a gene model and can enter the database",
             fontsize=13, ha="center", color="#333333")

    handles = [Patch(facecolor=ROUTE_COLOUR[r], edgecolor="none", label=ROUTE_LABEL[r])
               for r in PRECEDENCE]
    handles += [
        Line2D([], [], color="none", label=""),
        Patch(facecolor="#555555", edgecolor="none", label="solid: has a gene model (usable)"),
        Patch(facecolor=pale("#555555"), edgecolor="none",
              label="pale: assembly exists, no gene model"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=2,
               frameon=False, fontsize=12, handlelength=1.5, labelspacing=0.5,
               columnspacing=2.5)

    # Class labels are drawn OUTSIDE r_label and matplotlib does not grow the limits for text,
    # so the pad has to cover the longest label. "Eurotiomycetes (642)" sits at 180 degrees and
    # was clipped by the frame at +0.30.
    # ~0.030 data units per character at fontsize 10.5 on this canvas. Measured by eye against
    # the longest label rather than computed: matplotlib can only report text extents after a
    # draw, and a renderer round-trip is not worth it for one padding constant.
    longest = max((len(f"{n}  ({c})") for n, c in label_widths), default=20)
    lim = r_label + 0.030 * longest
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")
    return fig


def caption(species: list) -> str:
    """Generated from the data. A hand-written caption goes stale silently; this cannot."""
    n = len(species)
    n_ann = sum(1 for s in species if s["annotated"])
    asm = sum(s["total"] for s in species)
    ann = sum(s["annotated"] for s in species)
    by_route = defaultdict(int)
    for s in species:
        by_route[s["route"]] += 1
    routes = ", ".join(f"{ROUTE_LABEL[r].split(':')[0]} {by_route[r]:,}"
                       for r in PRECEDENCE if by_route[r])
    return (
        f"Reference scope of db_v4. Circular cladogram of the {n:,} fungal species for which "
        f"NCBI holds at least one assembly, drawn on the NCBI taxonomy held in holobase: "
        f"topology is taxonomic, branch lengths are not evolutionary distance, and radius "
        f"encodes taxonomic RANK (kingdom, phylum, class, order, family, genus, species) "
        f"rather than any measure of divergence. Tip colour gives the route that brought the "
        f"species "
        f"into scope, taking the strongest evidence where several apply ({routes}); congeners "
        f"are grey because they are present to give the classifier a competitor and are not "
        f"reportable as detections. Each tip carries a radial bar on a log10 scale, pale for "
        f"all {asm:,} assemblies that exist and solid for the {ann:,} that carry a gene model. "
        f"Both bars start at the same origin, so the pale remainder is the log-fold gap "
        f"between what exists and what a CDS-only database can use; the bars are not an "
        f"additive stack. {n_ann:,} species ({n_ann / n:.0%}) have at least one annotated "
        f"assembly and can enter db_v4; the {n - n_ann:,} entirely pale tips have assemblies "
        f"but no gene model and are the gene-prediction worklist. Classes with at least "
        f"{LABEL_MIN} species are labelled."
    )


def main():
    tax = _taxonomy.load()
    species = load_species(tax)
    print(f"  {len(species):,} species to draw")
    children, root, parent = build_tree(tax, species)
    tips = {s["tid"] for s in species}
    angle, radius, order = layout(tax, children, root, tips)
    print(f"  cladogram: {len(order):,} tips, {len(children):,} internal nodes, "
          f"root {tax.label(root)}")
    fig = draw(tax, species, children, root, angle, radius, order)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        # No bbox_inches="tight": the reserved legend band is deliberate and tight cropping
        # would undo it.
        fig.savefig(f"{OUT}.{ext}", dpi=300, facecolor="white")
    Path(f"{OUT}_caption.txt").write_text(caption(species) + "\n")
    print(f"  -> {OUT}.png / .pdf / _caption.txt")


if __name__ == "__main__":
    main()
