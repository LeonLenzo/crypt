#!/usr/bin/env python3
"""Assemble the reference cohort into one published table, one row per taxon.

    python release/02_reference/scripts/make_supp_tables.py

Reads  data/curation/screen_assemblies.tsv   what is buildable, and by which route
       data/curation/screen_species.tsv      USDA names resolved, and on_wheat
       data/curation/screen_pathogens.tsv    the PHI-base labels
       data/curation/screen_merges.tsv       provisional names folded onto their binomial
       data/curation/screen_hosts.tsv        the plant hosts, the database's sink
       data/curation/extra_cds_db2.tsv       the gene models we predicted ourselves
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

## THE HOSTS ARE ROWS HERE, AND THAT CHANGES WHAT THIS TABLE CLAIMS TO BE

The database contains plant CDS as a SINK: without it a plant read is absorbed by the nearest
fungal CDS and reported as a fungal hit, and these are plant RNA-seq libraries. A table that
called itself the reference while omitting wheat would misdescribe what was built, so the
hosts are published here rather than only in build_plan.tsv, which is not a supp_data_*
table and so would leave them out of the supplement entirely.

The consequence is that S3 is "THE DATABASE CONTENTS", not "the fungal scope". Host rows
carry `route=host`, `kingdom=Viridiplantae` and EMPTY fungal columns, because no clause of
the scope rule or the PHI-base label means anything for a wheat genome. **Anything counting
the scope funnel or drawing a panel must filter `kingdom == "Fungi"`** -- the figures do, and
the matrix cell report below does. `in_tree` is `no` for every host, so panel B is unaffected
either way.
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
         "plant_pathogen_evidence", "n_assemblies", "n_annotated", "n_predicted",
         "annotation_ratio",
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


    # Gene models PREDICTED here rather than downloaded (research/predict, vipsania), which
    # enter the db2 build through --extra-cds. They are published as their own column and
    # never folded into n_annotated: n_annotated is what NCBI carries, and a reader has to
    # be able to see which is which. Every one of these is a taxon with NO annotated
    # assembly, so n_predicted > 0 with n_annotated = 0 means the taxon became buildable
    # here and was uncallable before.
    n_pred = defaultdict(int)
    for r in rd(CUR / "extra_cds_db2.tsv"):
        n_pred[int(r["taxid"])] += 1

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

        # TWO VALUES, and the RANK OF THE RECORD IS NOT ONE OF THEM. Leon's call,
        # 2026-11-24: the matrix describes taxa, and a taxon is either named on wheat or it
        # is not. Whether it entered the reference because USDA named the species or because
        # USDA named its genus is how it got here, not what it is, and that is already
        # recorded in `route` for anyone who needs it. Three earlier attempts at this axis
        # (an observed/congener split, a genus_only row, a host x rank crossing) all tried
        # to draw the route and all made the figure assert something about the organism that
        # the route cannot support.
        usda_row = "wheat_resident" if onw else "plant_resident"
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
               "n_predicted": n_pred.get(t, 0),
               "annotation_ratio": round(n_ann[t] / n_asm[t], 4),
               "category": category,
               # A PREDICTED GENE MODEL PUTS A TAXON IN THE TREE. The rule used to be
               # "observation with an NCBI annotation", which excluded all 135 taxa
               # research/predict was run on, every one of which had no annotation at all:
               # the figure could not show the arm whose whole purpose is to close that
               # gap. The tree draws what the db2 database can be built from, so a
               # predicted model counts exactly as a downloaded one does for membership,
               # and `n_predicted` keeps them distinguishable once drawn.
               "in_tree": "yes" if (obs and (n_ann[t] or n_pred.get(t, 0))) else "no",
               "merged_from": ";".join(merged.get(t, ()))}
        row.update({r: lin.get(r, "") for r in RANKS})
        rows.append(row)

    # The host sink. Kept rows only; screen_hosts.tsv carries the dropped hosts and their
    # reasons as the audit trail, which is not a database input.
    hosts = CUR / "screen_hosts.tsv"
    n_host = 0
    if hosts.exists():
        for h in rd(hosts):
            if h["decision"] != "keep" or not h["ensembl_species"]:
                continue
            t = int(h["ncbi_taxid"])
            n = int(h["n_genomes"])
            lin = lineage(t, parent, rank, name)
            row = {c: "" for c in COLS}
            row.update({
                "ncbi_taxid": t, "taxon": name.get(t, h["host"]),
                "route": "host", "usda_row": "", "on_wheat": "",
                "phibase_col": "", "is_plant_pathogen": "", "is_wheat_pathogen": "",
                # Every host genome carries a gene model by construction: these come from
                # Ensembl Plants, where a genome is only published with a gene build, so
                # n_annotated == n_assemblies and the ratio is 1 for a host by definition.
                "n_assemblies": n, "n_annotated": n, "n_predicted": 0,
                "annotation_ratio": 1.0,
                "category": f"host_{h['role']}", "in_tree": "no", "merged_from": "",
            })
            row.update({r: lin.get(r, "") for r in RANKS})
            rows.append(row)
            n_host += 1
    else:
        print("  WARNING: no screen_hosts.tsv, so S3 omits the host sink. "
              "Run screen_hosts.py first.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    fungi = [r for r in rows if r["kingdom"] == "Fungi"]
    cells = defaultdict(lambda: [0, 0, 0])
    for r in fungi:
        c = cells[(r["usda_row"], r["phibase_col"])]
        c[0] += 1
        c[1] += r["n_assemblies"]
        c[2] += r["n_annotated"]
    print(f"{len(fungi):,} fungal taxa, {sum(r['n_assemblies'] for r in fungi):,} "
          f"assemblies, {sum(r['n_annotated'] for r in fungi):,} with a gene model")
    for k in sorted(cells):
        t, a, c = cells[k]
        print(f"  {k[0]:<9} x {k[1]:<15} {t:>6,} taxa {a:>7,} asm {c:>6,} cds")
    print(f"  in_tree {sum(1 for r in rows if r['in_tree'] == 'yes'):,}")
    if n_host:
        print(f"\n  + {n_host} HOST taxa (the sink, kingdom Viridiplantae, excluded from "
              f"every count above):")
        for r in rows[-n_host:]:
            print(f"      {r['category']:<18}{r['taxon']:<24}{r['genus']}")
    print(f"  {len(rows):,} rows total in the table")
    missing = [r["taxon"] for r in fungi if not r["genus"]]
    if missing:
        print(f"  note: {len(missing)} taxa have no genus in the lineage: "
              + ", ".join(missing[:5]))
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
