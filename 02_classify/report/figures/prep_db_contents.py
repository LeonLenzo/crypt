#!/usr/bin/env python3
"""prep_db_contents.py — the data behind the db_v4 contents figure, and its caption.

    python 02_classify/report/figures/prep_db_contents.py
    Rscript 02_classify/report/figures/db_contents.R

Everything that decides WHAT the figure shows lives here; `db_contents.R` decides how it
looks and nothing else. That split is the house pattern (`01_literature` does prep in Python
and the figure in R) and it is also why this file exists at all: the figure was first written
in matplotlib, and keeping both renderers would have meant two programs writing one PNG and
free to disagree about it. The R version is a quarter the length and is the one that stayed.

Writes four files:

    db_contents_tree.nwk      the cladogram, Newick, ALL BRANCH LENGTHS 1
    db_contents_tips.tsv      one row per tip: label, species, route, assembly counts, class
    db_contents_clades.tsv    class and order spans, for labelling in Illustrator
    db_contents_caption.txt   the full generated caption

**Every branch length is 1 because the tree is a RANK LADDER.** NCBI's uneven nesting is
collapsed to exactly seven ranks -- kingdom, phylum, class, order, family, genus, species --
so every edge is one rank step, unit lengths give evenly spaced concentric rings with no
further work, and radius means taxonomic RANK rather than divergence. That is also why this
is a taxonomic cladogram and must never be captioned as a phylogeny.

The colours are NOT here. They live in `db_contents.R` beside the rest of the styling, so
there is one place to edit them; this file only needs the route NAMES.
"""

import argparse
import csv
import math
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _layout import SILVER, FIGURES, HOLOBASE, HOSTS   # noqa: E402
from _seeds import read_seed_tsv                       # noqa: E402
import _taxonomy                                       # noqa: E402

SURVEY = SILVER / "assembly_survey.tsv"
OUT = FIGURES / "db_contents"

# Route NAMES and their precedence. One colour per tip in the figure, chosen as the strongest
# evidence that applies; `target` is split by whether the pathogen is recorded on one of the
# four cohort cereals. The hex values are in db_contents.R.
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

GAP_DEG = 14.0     # ring left open at the top; must match OPEN_ANGLE in db_contents.R
LABEL_MIN = 40     # a class needs this many species to earn a label

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
    """An angle for every node, in ring order. Tips in traversal order stay contiguous."""
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
    # No radii here. Radial position is ggtree's job: it comes from the unit branch lengths in
    # the Newick, which make the seven ranks evenly spaced rings, and from HOLE in
    # db_contents.R. This function exists to put the tips in ring order and give each one an
    # angle, which is what the clade table reports.
    return angle, order


def clade_table(tax, species, angle, order, path: Path) -> list:
    """Every clade's name, size and angular position, in ring order, at class and order rank.

    This is what makes assembly mode usable: with the labels stripped there is nothing on the
    image saying which wedge is which, so the positions have to leave the script as DATA
    rather than be measured off the PNG by hand. Angles are degrees anticlockwise from east,
    with the ring starting and ending at the gap centred on 90 degrees.
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
            unplaced: int, clades: list) -> str:
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

Drawn by db_contents.R. In its default ASSEMBLY mode all text except the log-axis numerals
is stripped from the image, and the class names, counts and angular positions below are what
you place in Illustrator; `--self-contained` draws the title and legend instead.


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
as WEDGES -- a geom_rect in tree coordinates, which ggtree's fan layout turns into a wedge
that tiles its angular slot exactly. A mark of fixed width in POINTS cannot do that: the slot
grows with radius while the mark does not, and at this tip count a fixed-width bar was 34%
wider than its slot and overdrew its neighbours.

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




def newick(children: dict, root: int, tax) -> str:
    """Newick with unit branch lengths, iteratively.

    Iterative rather than recursive on purpose: the tree is ~2,400 nodes and a recursive
    writer needs the interpreter's limit raised, which is a trap waiting for a deeper
    taxonomy. Labels are `t<internal taxon id>`, unique by construction and the join key to
    the tip table; a scientific name is not safe in Newick, which gives parentheses, commas
    and colons their own meaning.
    """
    out, stack = {}, [(root, False)]
    while stack:
        node, done = stack.pop()
        kids = sorted(children.get(node, ()), key=lambda k: tax.label(k))
        if not kids:
            out[node] = f"t{node}:1"
            continue
        if done:
            out[node] = "(" + ",".join(out[k] for k in kids) + f")t{node}:1"
        else:
            stack.append((node, True))
            stack.extend((k, False) for k in kids)
    return out[root].rsplit(":", 1)[0] + ";"


def main():
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args()
    tax = _taxonomy.load()
    on_cereal = cereal_recorded(tax)
    species, no_model, unplaced = load_species(tax, on_cereal)
    children, root, parent = build_tree(tax, species)
    tips = {s["tid"] for s in species}
    angle, order = layout(tax, children, root, tips)
    internal = sum(1 for n in children if n not in tips)
    print(f"  {len(order):,} tips, {internal:,} internal nodes, root {tax.label(root)}")

    Path(f"{OUT}_tree.nwk").write_text(newick(children, root, tax) + "\n")

    by_tid = {s["tid"]: s for s in species}
    cols = ["label", "ncbi_taxid", "species", "route", "n_annotated", "n_assemblies",
            "class", "order", "on_cohort_cereal"]
    with open(f"{OUT}_tips.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for tid in order:
            s = by_tid[tid]
            cls, ordr = tax.ancestor_at(tid, "class"), tax.ancestor_at(tid, "order")
            w.writerow({
                "label": f"t{tid}", "ncbi_taxid": tax.ncbi_of(tid), "species": s["name"],
                "route": s["route"], "n_annotated": s["annotated"],
                "n_assemblies": s["total"],
                "class": tax.label(cls) if cls else "",
                "order": tax.label(ordr) if ordr else "",
                "on_cohort_cereal": "yes" if tid in on_cereal else "",
            })

    clades = clade_table(tax, species, angle, order, Path(f"{OUT}_clades.tsv"))
    n_cereal = sum(1 for s in species if s["route"] == "target_cereal")
    n_cereal_other = sum(1 for s in species
                         if s["tid"] in on_cereal and s["route"] != "target_cereal")
    Path(f"{OUT}_caption.txt").write_text(
        caption(species, no_model, n_cereal, n_cereal_other, unplaced, clades))
    print(f"  {n_cereal:,} drawn species are pathogens recorded on a cohort cereal; "
          f"{n_cereal_other:,} more are recorded there by another route")
    print(f"  -> {OUT}_tree.nwk / _tips.tsv / _clades.tsv / _caption.txt")
    print(f"\nNext: Rscript 02_classify/report/figures/db_contents.R")


if __name__ == "__main__":
    main()
