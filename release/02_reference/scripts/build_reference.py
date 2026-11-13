#!/usr/bin/env python3
"""Assemble the reference cohort into one published table, one row per taxon.

    python release/02_reference/scripts/build_reference.py

Reads  data/curation/screen_assemblies.tsv   what is buildable, and by which route
       data/curation/screen_species.tsv      USDA names resolved, and on_wheat
       data/curation/screen_pathogens.tsv    the PHI-base labels
       data/curation/screen_merges.tsv       provisional names folded onto their binomial
       data/raw/raw_taxdmp.zip               nodes.dmp + names.dmp, for the lineage
Writes data/outputs/supp_data_s3_reference.tsv

THE LAST STEP OF THE MODULE, and the only one the figures read. `release/` publishes its
result as a numbered supplementary table and everything downstream reads that, the way
01_cohort's figures read supp_data_s2_cohort.tsv rather than reaching back into curation.
The curation tables stay as the audit trail; this is the deliverable.

Self-sufficiency is the design goal: the table carries the six parent ranks, so panel B's
topology is rebuilt from these columns and nothing in `figures/` needs the taxdump. Someone
holding only this TSV can regenerate both panels.

    matrix   group by usda_row x phibase_col, count rows, sum n_assemblies and n_annotated
    tree     the rows with in_tree=yes; topology from kingdom..genus; bar = n_annotated;
             ring = annotation_ratio; arcs = class

One row per taxon IN THE REFERENCE, so a taxon that carries a label but has no in-scope
assembly is absent by construction (one wheat pathogen is in that position). The pathogen
labels are a property of the taxon and travel with it; they are never a filter on scope.
"""

import argparse
import csv
import io
import zipfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
CUR = HERE / "data/curation"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
OUT = HERE / "data/outputs/supp_data_s3_reference.tsv"

RANKS = ("kingdom", "phylum", "class", "order", "family", "genus")
COLS = (["ncbi_taxid", "taxon"] + list(RANKS) +
        ["route", "usda_row", "on_wheat",
         "phibase_col", "is_plant_pathogen", "is_wheat_pathogen", "wheat_pathogen_source",
         "plant_pathogen_evidence", "n_assemblies", "n_annotated", "annotation_ratio",
         "category", "in_tree", "merged_from"])
csv.field_size_limit(10 ** 7)


def rd(p):
    with open(p) as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                   delimiter="\t"))


def load_taxonomy():
    parent, rank, name = {}, {}, {}
    with zipfile.ZipFile(TAXDMP) as zf:
        with zf.open("nodes.dmp") as raw:
            for ln in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
                p = [f.strip() for f in ln.split("|")]
                parent[int(p[0])], rank[int(p[0])] = int(p[1]), p[2]
        with zf.open("names.dmp") as raw:
            for ln in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
                p = [f.strip() for f in ln.split("|")]
                if len(p) >= 4 and p[3] == "scientific name":
                    name[int(p[0])] = p[1]
    return parent, rank, name


def lineage(tid, parent, rank, name):
    """-> {rank: name} for the six parent ranks. Cycle-guarded: root is its own parent."""
    out, seen, cur = {}, set(), tid
    while cur and cur not in seen:
        seen.add(cur)
        r = rank.get(cur)
        if r in RANKS and r not in out:
            out[r] = name.get(cur, "")
        nxt = parent.get(cur)
        if nxt is None or nxt == cur:
            break
        cur = nxt
    return out


def species_of(tid, parent, rank):
    seen, cur = set(), tid
    while cur and cur not in seen:
        if rank.get(cur) == "species":
            return cur
        seen.add(cur)
        nxt = parent.get(cur)
        if nxt is None or nxt == cur:
            return cur
        cur = nxt
    return cur


def main():
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args()
    parent, rank, name = load_taxonomy()

    n_asm, n_ann, routes, organism = (defaultdict(int), defaultdict(int),
                                      defaultdict(set), {})
    for r in rd(CUR / "screen_assemblies.tsv"):
        t = int(r["species_taxid"])
        n_asm[t] += 1
        routes[t].add(r["route"])
        organism.setdefault(t, r["organism_name"])
        if r["annotated"] == "yes":
            n_ann[t] += 1
    scope = set(n_asm)

    on_wheat = {species_of(int(r["ncbi_taxid"]), parent, rank)
                for r in rd(CUR / "screen_species.tsv")
                if r["on_wheat"] == "yes" and r["ncbi_taxid"]}

    pth = {int(r["ncbi_taxid"]): r for r in rd(CUR / "screen_pathogens.tsv")}

    merged = defaultdict(list)
    for r in rd(CUR / "screen_merges.tsv"):
        if r["action"] == "merge" and r["to_taxid"].isdigit():
            merged[int(r["to_taxid"])].append(r["from_taxid"])

    rows = []
    for t in sorted(scope, key=lambda x: name.get(x, organism.get(x, str(x)))):
        obs = bool(routes[t] & {"usda_species", "exception"})
        onw = obs and t in on_wheat
        p = pth.get(t, {})
        is_plant = p.get("is_plant_pathogen") == "yes"
        is_wheat = p.get("is_wheat_pathogen") == "yes"

        usda_row = "wheat" if onw else ("plant" if obs else "congener")
        phibase_col = ("wheat_pathogen" if is_wheat else
                       "plant_pathogen" if is_plant else "none")
        # The five tree categories ARE the five occupied observation cells of the matrix.
        category = "" if not obs else (
            "wheat_pathogen" if is_wheat else
            ("onwheat_pathogen" if onw else "plant_pathogen") if is_plant else
            ("onwheat_resident" if onw else "plant_resident"))

        lin = lineage(t, parent, rank, name)
        row = {"ncbi_taxid": t, "taxon": name.get(t, organism.get(t, "")),
               "route": "+".join(sorted(routes[t])), "usda_row": usda_row,
               "on_wheat": "yes" if onw else "no",
               "phibase_col": phibase_col,
               "is_plant_pathogen": "yes" if is_plant else "no",
               "is_wheat_pathogen": "yes" if is_wheat else "no",
               "wheat_pathogen_source": p.get("wheat_pathogen_source", ""),
               "plant_pathogen_evidence": p.get("plant_pathogen_evidence", ""),
               "n_assemblies": n_asm[t], "n_annotated": n_ann[t],
               "annotation_ratio": round(n_ann[t] / n_asm[t], 4),
               "category": category,
               "in_tree": "yes" if (obs and n_ann[t]) else "no",
               "merged_from": ";".join(merged.get(t, ()))}
        row.update({r: lin.get(r, "") for r in RANKS})
        rows.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    cells = defaultdict(lambda: [0, 0, 0])
    for r in rows:
        c = cells[(r["usda_row"], r["phibase_col"])]
        c[0] += 1
        c[1] += r["n_assemblies"]
        c[2] += r["n_annotated"]
    print(f"{len(rows):,} taxa, {sum(r['n_assemblies'] for r in rows):,} assemblies, "
          f"{sum(r['n_annotated'] for r in rows):,} with a gene model")
    for k in sorted(cells):
        t, a, c = cells[k]
        print(f"  {k[0]:<9} x {k[1]:<15} {t:>6,} taxa {a:>7,} asm {c:>6,} cds")
    print(f"  in_tree {sum(1 for r in rows if r['in_tree'] == 'yes'):,}")
    missing = [r["taxon"] for r in rows if not r["genus"]]
    if missing:
        print(f"  note: {len(missing)} taxa have no genus in the lineage: "
              + ", ".join(missing[:5]))
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
