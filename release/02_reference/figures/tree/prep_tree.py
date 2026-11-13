#!/usr/bin/env python3
"""prep_tree.py — the data behind the reference-contents tree (panel B of the figure).

    python release/02_reference/figures/tree/prep_tree.py
    Rscript release/02_reference/figures/tree/tree.R

Reads  ../../data/outputs/supp_data_s3_reference.tsv   the published reference table
Writes tree_tree.nwk, tree_tips.tsv, tree_clades.tsv, tree_params.tsv (beside this script)

THE TAXDUMP IS NOT READ HERE. S3 carries the six parent ranks per taxon, so the topology is
rebuilt from those columns and the figure depends only on the published deliverable. Someone
holding the TSV can regenerate the panel.

One tip per OBSERVATION species (USDA-recorded, annotated assembly = db0); congeners are
excluded, since they are competition-only and never reported. The tree is the NCBI TAXONOMY, collapsed
to seven rank rings (kingdom..species), so radius encodes RANK, not divergence -- a taxonomic
cladogram, never a phylogeny.

Tip colour is the taxon's CELL IN PANEL A, so the two panels share one colour vocabulary
exactly. The five categories are the five occupied observation cells of the matrix, which is
why there are five and not four: a plant pathogen USDA records on wheat is a different cell
from one it does not, and the matrix distinguishes them, so the tree must too.

    wheat_pathogen     Wheat x Wheat pathogen     PHI-base names a disease on Triticum
    onwheat_pathogen   Wheat x Plant pathogen     pathogen elsewhere, USDA records it on wheat
    onwheat_resident   Wheat x Non pathogen       on wheat, no pathogenicity record
    plant_pathogen     Plant x Plant pathogen     pathogen elsewhere, never recorded on wheat
    plant_resident     Plant x Non pathogen       the community the mining is for

Congeners have no row here at all: they are competition-only and never reported.

Bar length (drawn in R) is the species' GENOME count on a log scale, so the hyper-sequenced
species stand out as spikes and the "take every genome, no cap" decision is visible.

The algorithm (collapse to ranked rings, rank-gap branch lengths, fan layout) follows the
db_v4 db_contents figure; the taxonomy source and the colour scheme are what differ.
"""

import csv
import math
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]          # release/02_reference
S3 = HERE / "data/outputs/supp_data_s3_reference.tsv"
FIG = Path(__file__).resolve().parent
OUT = FIG / "tree"
csv.field_size_limit(10 ** 7)

RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
RANK_LEVEL = {r: i + 1 for i, r in enumerate(RANKS)}
GAP_DEG = 14.0
LABEL_MIN = 40

# Draw/stack order is set in tree.R; this is the precedence, most specific first.
CATEGORIES = ("wheat_pathogen", "onwheat_pathogen", "onwheat_resident",
              "plant_pathogen", "plant_resident")


def ratio_colour(r):
    """annotation ratio -> hex. r=1 (all assemblies annotated) pale; r->0 dark purple."""
    pale, dark = (242, 239, 250), (84, 39, 143)
    rgb = tuple(round(dark[i] + (pale[i] - dark[i]) * r) for i in range(3))
    return "#%02x%02x%02x" % rgb


class Tax:
    """Taxonomy assembled from S3's lineage columns, with the same interface the tree code
    used to get from the taxdump.

    Node ids are mixed on purpose: a TIP is its integer NCBI taxid, an INTERNAL node is the
    string "rank:name". S3 names the parent ranks without numbering them, and a name is
    unique within a rank, so the pair identifies the node. The rank prefix guards the one
    collision that could otherwise happen, a class and an order sharing a name.
    """

    def __init__(self, rows):
        self.parent, self.rank, self.name = {}, {}, {}
        for r in rows:
            chain = []
            for rk in RANKS[:-1]:                      # kingdom .. genus
                v = r[rk].strip()
                if not v:
                    continue                           # a rank NCBI leaves unplaced
                nid = f"{rk}:{v}"
                self.rank[nid], self.name[nid] = rk, v
                chain.append(nid)
            tid = int(r["ncbi_taxid"])
            self.rank[tid], self.name[tid] = "species", r["taxon"]
            chain.append(tid)
            for a, b in zip(chain, chain[1:]):
                self.parent[b] = a
            self.parent.setdefault(chain[0], None)

    def label(self, t):
        return self.name.get(t, str(t))

    def ancestor_at(self, t, want):
        seen = set()
        while t is not None and t not in seen:
            if self.rank.get(t) == want:
                return t
            seen.add(t)
            t = self.parent.get(t)
        return None


def load_species(rows):
    """The tree's tips: the S3 rows flagged in_tree, i.e. annotated observations."""
    out = []
    for r in rows:
        if r["in_tree"] != "yes":
            continue
        ratio = float(r["annotation_ratio"])
        out.append({"tid": int(r["ncbi_taxid"]), "name": r["taxon"],
                    "n": int(r["n_annotated"]), "n_total": int(r["n_assemblies"]),
                    "ratio": round(ratio, 4), "band": ratio_colour(ratio),
                    "cat": r["category"]})
    bad = {s["cat"] for s in out} - set(CATEGORIES)
    assert not bad, f"S3 carries a category the tree does not know: {bad}"
    return out


def build_tree(tax, species):
    tips = {s["tid"] for s in species}
    children = defaultdict(set)
    root = None
    for tid in tips:
        chain, cur, seen = [], tid, set()
        while cur is not None and cur not in seen:
            seen.add(cur)
            chain.append(cur)
            nxt = tax.parent.get(cur)
            if nxt is None or nxt == cur:
                break
            cur = nxt
        chain = list(reversed(chain))
        root = chain[0] if root is None else root
        for a, b in zip(chain, chain[1:]):
            children[a].add(b)

    def keep(node):
        return node in tips or tax.rank.get(node) in RANK_LEVEL

    def ranked_kids(node):
        out = []
        for kid in children.get(node, ()):
            out.append(kid) if keep(kid) else out.extend(ranked_kids(kid))
        return out

    rebuilt, stack = {}, [root]
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

    reachable, stack = {}, [root]
    while stack:
        node = stack.pop()
        if node in reachable:
            continue
        reachable[node] = set(children.get(node, ()))
        stack.extend(reachable[node])
    return reachable, root


def layout(tax, children, root, tips):
    order = []

    def walk(node):
        kids = sorted(children.get(node, ()), key=tax.label)
        if not kids or node in tips:
            order.append(node)
            return
        for kid in kids:
            walk(kid)

    walk(root)
    span = math.radians(360.0 - GAP_DEG)
    start = math.radians(90.0 + GAP_DEG / 2.0)
    angle = {t: start + span * (i + 0.5) / len(order) for i, t in enumerate(order)}

    def set_angle(node):
        if node in angle:
            return angle[node]
        vals = [set_angle(k) for k in children.get(node, ())]
        angle[node] = (min(vals) + max(vals)) / 2.0
        return angle[node]

    set_angle(root)
    return angle, order


def rank_level(node, tips, tax):
    if node in tips:
        return len(RANKS)
    return RANK_LEVEL.get(tax.rank.get(node), 0)


def newick(children, root, tips, tax):
    # Node LABELS must be newick-safe. A TIP is an int taxid, which already is: keep
    # "t<taxid>", because tree_tips.tsv joins to the R side on exactly that string. An
    # INTERNAL node is now "rank:name", carrying a space and a colon, so number those
    # instead; their names travel in tree_clades.tsv, which is where the R side reads them.
    nid = {}

    def tag(n):
        if isinstance(n, int):
            return f"t{n}"
        return f"n{nid.setdefault(n, len(nid))}"

    level, out, stack = {}, {}, [(root, False)]
    while stack:
        node, done = stack.pop()
        level.setdefault(node, rank_level(node, tips, tax))
        kids = sorted(children.get(node, ()), key=tax.label)
        if not kids:
            continue
        if done:
            parts = [f"{out.get(k, tag(k))}:{max(1, level[k] - level[node])}" for k in kids]
            out[node] = "(" + ",".join(parts) + f"){tag(node)}"
        else:
            stack.append((node, True))
            stack.extend((k, False) for k in kids)
    for node in level:
        out.setdefault(node, tag(node))
    return out[root] + ";"


def clade_table(tax, order, angle, path):
    rows = []
    for rank in ("class", "order"):
        groups = defaultdict(list)
        for tid in order:
            anc = tax.ancestor_at(tid, rank)
            groups[tax.label(anc) if anc else "(unplaced)"].append(math.degrees(angle[tid]))
        for name, degs in groups.items():
            lo, hi = min(degs), max(degs)
            rows.append({"rank": rank, "clade": name, "n_species": len(degs),
                         "start_deg": round(lo % 360, 2), "end_deg": round(hi % 360, 2),
                         "mid_deg": round(((lo + hi) / 2) % 360, 2),
                         "span_deg": round(hi - lo, 2),
                         "wraps_zero": "yes" if (lo % 360) > (hi % 360) else ""})
    rows.sort(key=lambda r: (r["rank"], r["start_deg"]))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    with open(S3) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    tax = Tax([r for r in rows if r["in_tree"] == "yes"])
    species = load_species(rows)
    children, root = build_tree(tax, species)
    tips = {s["tid"] for s in species}
    angle, order = layout(tax, children, root, tips)
    internal = sum(1 for n in children if n not in tips)
    print(f"  {len(order):,} tips, {internal:,} internal nodes, root {tax.label(root)}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{OUT}_tree.nwk").write_text(newick(children, root, tips, tax) + "\n")

    by = {s["tid"]: s for s in species}
    with open(f"{OUT}_tips.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["label", "species", "n_genomes", "n_total", "annotation_ratio",
                    "band_colour", "category", "class"])
        for tid in order:
            s = by[tid]
            cls = tax.ancestor_at(tid, "class")
            w.writerow([f"t{tid}", s["name"], s["n"], s["n_total"], s["ratio"], s["band"],
                        s["cat"], tax.label(cls) if cls else ""])

    clade_table(tax, order, angle, f"{OUT}_clades.tsv")
    with open(f"{OUT}_params.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["param", "value"])
        w.writerow(["open_angle_deg", GAP_DEG])

    cat_n = defaultdict(int)
    for s in species:
        cat_n[s["cat"]] += 1
    print("  categories:", {c: cat_n[c] for c in CATEGORIES})
    print(f"  -> {OUT}_tree.nwk / _tips.tsv / _clades.tsv / _params.tsv")


if __name__ == "__main__":
    main()
