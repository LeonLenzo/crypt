#!/usr/bin/env python3
"""Select the assemblies inside the scope, and record what each one can contribute.

    python release/02_reference/scripts/screen_assemblies.py

Reads  data/inputs/scope_exceptions.tsv     the genera admitted outside the master list
       data/curation/screen_species.tsv     USDA names resolved to NCBI taxids
       data/curation/screen_merges.tsv      provisional NCBI names folded onto their binomial
       data/raw/raw_taxdmp.zip              nodes.dmp, for the walk to genus
       data/raw/raw_assembly_summary.txt    every GenBank fungal assembly
Writes data/curation/screen_assemblies.tsv  one row per in-scope assembly

The scope is read at GENUS: an assembly is in if its species sits inside a genus the master
list reached, or inside one of the exception genera. `route` records which, because a taxon
present only because its genus qualified must not be reportable the way a named one is.

NOTHING IS FILTERED OUT HERE. Unannotated assemblies stay in the table with annotated=no,
because they are the gene-prediction worklist and dropping them would hide it. The columns
that decide that worklist travel with each row: `species_annotated` says whether the species
already has a gene model from anywhere, and `genus_annotated` counts the annotated
assemblies in its genus, which is what separates a genuine gap from redundancy. Predicting
into a genus that already holds several annotated assemblies mostly adds competitors.

Annotation comes from NCBI's own columns, never inferred. `annotation_provider` is `na` when
absent, and an assembly without one has no CDS to extract.

`screen_merges.tsv` is applied FIRST, before the genus walk and before `route`. NCBI files
some genomes under a provisional node (`Cercospora cf. flagellaris`) whose binomial USDA
names separately, which would otherwise split one fungus into an observation with no genome
and a congener with no host record. A merged assembly carries the TARGET taxid and name, so
every count downstream sees one taxon. Written by `screen_merges.py`; absent file = no
merging, and the previous behaviour.
"""

import argparse
import csv
import io
import zipfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
EXC = HERE / "data/inputs/scope_exceptions.tsv"
SPECIES = HERE / "data/curation/screen_species.tsv"
MERGES = HERE / "data/curation/screen_merges.tsv"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
SUMMARY = HERE / "data/raw/raw_assembly_summary.txt"
OUT = HERE / "data/curation/screen_assemblies.tsv"
LEVEL_ORDER = {"Complete Genome": 3, "Chromosome": 2, "Scaffold": 1, "Contig": 0}

csv.field_size_limit(10 ** 7)


def read_tsv(path):
    with open(path) as fh:
        yield from csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                  delimiter="\t")


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


def ancestor_at(tid, want, parent, rank):
    """Cycle-guarded: NCBI's root is its own parent."""
    seen = set()
    while tid and tid not in seen:
        if rank.get(tid) == want:
            return tid
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt == tid:
            return None
        tid = nxt
    return None


def load_merges():
    """-> ({from_taxid: to_taxid}, {taxids promoted to observation in place})."""
    merge, reroute = {}, set()
    if not MERGES.exists():
        print(f"  note: {MERGES.name} not found; no provisional names will be folded. "
              f"Run screen_merges.py first.")
        return merge, reroute
    for r in read_tsv(MERGES):
        if r["action"] == "merge" and r["to_taxid"].isdigit():
            merge[int(r["from_taxid"])] = int(r["to_taxid"])
        elif r["action"] == "reroute":
            reroute.add(int(r["from_taxid"]))
    return merge, reroute


def assembly_rows():
    with open(SUMMARY, encoding="utf-8", errors="replace") as fh:
        hdr = None
        for ln in fh:
            if ln.startswith("#assembly_accession"):
                hdr = ln[1:].rstrip("\n").split("\t")
                continue
            if ln.startswith("#") or hdr is None:
                continue
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def main():
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args()
    parent, rank, name = load_taxonomy()
    at = lambda t, w: ancestor_at(t, w, parent, rank)

    usda_sp, wheat_sp = set(), set()
    for r in read_tsv(SPECIES):
        if r["ncbi_taxid"]:
            usda_sp.add(int(r["ncbi_taxid"]))
            if r["on_wheat"] == "yes":
                wheat_sp.add(int(r["ncbi_taxid"]))
    usda_gn = {g for g in (at(t, "genus") for t in usda_sp) if g}
    wheat_gn = {g for g in (at(t, "genus") for t in wheat_sp) if g}
    exc_gn = {int(r["ncbi_taxid"]) for r in read_tsv(EXC) if r["status"] == "in"}
    merge, reroute = load_merges()
    print(f"scope: {len(usda_gn):,} USDA genera + {len(exc_gn)} exception genera"
          f" | {len(merge)} provisional taxa merged, {len(reroute)} rerouted")

    annotated = lambda d: d.get("annotation_provider", "").strip() not in ("", "na")
    kept, n_merged = [], [0]
    for d in assembly_rows():
        if not d["species_taxid"].isdigit():
            continue
        sp = int(d["species_taxid"])
        was = sp
        sp = merge.get(sp, sp)
        gn = at(sp, "genus")
        if gn in exc_gn:
            route = "exception"
        elif sp in usda_sp or sp in reroute:
            route = "usda_species"
        elif gn in usda_gn:
            route = "usda_genus"
        else:
            continue
        if sp != was:
            n_merged[0] += 1
            d["organism_name"] = name.get(sp, d["organism_name"])
        d["_sp"], d["_gn"], d["_route"] = sp, gn, route
        kept.append(d)

    sp_ann, sp_n, gn_ann = defaultdict(int), defaultdict(int), defaultdict(int)
    for d in kept:
        sp_n[d["_sp"]] += 1
        if annotated(d):
            sp_ann[d["_sp"]] += 1
            gn_ann[d["_gn"]] += 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["accession", "species_taxid", "organism_name", "genus_taxid", "genus", "route",
            "on_wheat", "assembly_level", "annotated", "annotation_provider",
            "protein_coding_genes", "genome_size_mb", "contig_count", "scaffold_count",
            "seq_rel_date", "n_asm_species", "species_annotated", "genus_annotated"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for d in sorted(kept, key=lambda x: (x["organism_name"], x["assembly_accession"])):
            size = d.get("genome_size", "")
            w.writerow({
                "accession": d["assembly_accession"], "species_taxid": d["_sp"],
                "organism_name": d["organism_name"], "genus_taxid": d["_gn"] or "",
                "genus": name.get(d["_gn"], ""), "route": d["_route"],
                "on_wheat": "yes" if d["_gn"] in wheat_gn else "no",
                "assembly_level": d["assembly_level"],
                "annotated": "yes" if annotated(d) else "no",
                "annotation_provider": d.get("annotation_provider", ""),
                "protein_coding_genes": d.get("protein_coding_gene_count", ""),
                "genome_size_mb": f"{int(size) / 1e6:.1f}" if size.isdigit() else "",
                "contig_count": d.get("contig_count", ""),
                "scaffold_count": d.get("scaffold_count", ""),
                "seq_rel_date": d.get("seq_rel_date", ""),
                "n_asm_species": sp_n[d["_sp"]],
                "species_annotated": "yes" if sp_ann[d["_sp"]] else "no",
                "genus_annotated": gn_ann[d["_gn"]],
            })

    n_ann = sum(1 for d in kept if annotated(d))
    species = set(sp_n)
    need = {s for s in species if not sp_ann[s]}
    orphan = {s for s in need if gn_ann[at(s, 'genus')] == 0}
    print(f"\n  {len(kept):,} assemblies in scope over {len(species):,} taxa"
          f"  ({n_merged[0]} assemblies moved onto a merged taxon)")
    print(f"    annotated                      {n_ann:,} assemblies")
    print(f"    taxa with a gene model         {len(species) - len(need):,}")
    print(f"    taxa with none                 {len(need):,}")
    print(f"      genus has no annotation      {len(orphan):,}  <- the prediction gap")
    print(f"      genus already has some       {len(need) - len(orphan):,}  "
          f"(predicting here mostly adds competitors)")
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
