#!/usr/bin/env python3
"""Turn the in-scope assemblies into a tiered build plan for the three databases.

    python release/02_reference/scripts/plan_build.py
    python release/02_reference/scripts/plan_build.py --per-species 3   # cap per species

Reads  data/curation/screen_assemblies.tsv  every in-scope assembly, annotated or not
       data/raw/raw_taxdmp.zip             nodes.dmp, for the genus walk
Writes data/curation/build_plan.tsv         one row per SELECTED assembly, with its db tier

It sits in `curation/` and not `outputs/` because it is a DECISION, the fifth screen of the
module, and `outputs/` is the published supplement that the figures read. The plan is read by
fetch_cds.py and build_db.py, never by a figure.

The chapter builds three databases as an incremental ladder (kraken2-build --add-to-library
reuses the previous library), and each tier isolates one effect:

    host  the sink    plant CDS, in EVERY database (screen_hosts.py). Not a tier of the
                      ladder: it is present in all three, so it never enters a subtraction.
    db0  baseline    observations: route in {usda_species, exception}, annotated
    db1  +congeners  db0 + route usda_genus, annotated      (a sibling of a recorded genus)
    db2  +predicted  db1 + gene-predicted CDS, 135 genomes (research/predict/worklist.py)

THE HOST TIER IS NOT PART OF THE LADDER AND MUST NOT BE SUBTRACTED. db0/db1 were built
without any plant CDS at all (`Fungi 100%` in db_inspect.tsv), which is a regression: without
a host sink a plant read is absorbed by the nearest fungal CDS and reported as a fungal hit,
and these are plant RNA-seq libraries. Hosts go into every database so that db1-db0 and
db2-db1 still isolate one effect each; a comparison where one arm has the sink and the other
does not measures the sink, not the tier.

The tiers touch disjoint genera, which is what makes the subtraction clean. A congener is an
annotated species, so its genus is never an orphan; a predicted species is in a genus with no
annotation at all. So db1-db0 measures competition among already-represented genera and
db2-db1 measures reach into genera nothing could represent before. db2's predicted set is NOT
emitted here: it is genomic FASTA fetched and predicted in research/predict, then added to the
library, and predicted CDS is a COMPETITOR only, never a positive call. This plan covers the
CDS that exist to download, tiers db0 and db1.

ONE ASSEMBLY PER SPECIES by default. db_v4 took several per species for pangenome breadth, but
the saturation work found open pangenomes with marginal return ~gamma/n (median gamma 0.162),
so a second assembly adds ~16% k-mer content for 2x the fetch. --per-species raises the cap;
selection is by assembly level then fewest contigs.

Nothing about Acacia is consulted here: this is the scope the build NEEDS. What is already
downloaded is a fetch optimisation, applied by the fetch step with --have.
"""

import argparse
import csv
import io
import zipfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
SCREEN = HERE / "data/curation/screen_assemblies.tsv"
HOSTS = HERE / "data/curation/screen_hosts.tsv"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
OUT = HERE / "data/curation/build_plan.tsv"
LEVEL = {"Complete Genome": 3, "Chromosome": 2, "Scaffold": 1, "Contig": 0}
# CDS as a fraction of genome, used only to SIZE the download. Fungal genomes are gene-dense
# and the db_v4 ratio is ~0.50. PLANTS ARE NOT: a hexaploid wheat genome is mostly repeat, and
# its CDS is ~0.2 GB against ~15 GB genomic (CLAUDE.md), so applying the fungal ratio to a
# host overstates it by ~40x and made the host sink look like 9.9 GB when it is nearer 0.3.
CDS_RATIO = {"fungus": 0.50, "host": 0.013}

csv.field_size_limit(10 ** 7)


def genus_of(tid, parent, rank):
    seen = set()
    while tid and tid not in seen:
        if rank.get(tid) == "genus":
            return tid
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt == tid:
            return None
        tid = nxt
    return None


def load_nodes():
    parent, rank = {}, {}
    with zipfile.ZipFile(TAXDMP) as zf, zf.open("nodes.dmp") as raw:
        for ln in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
            p = [f.strip() for f in ln.split("|")]
            parent[int(p[0])], rank[int(p[0])] = int(p[1]), p[2]
    return parent, rank


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--per-species", default="all",
                    help="max assemblies kept per species: an integer, or 'all' (default). "
                         "db_v4 took all; the per-species ceiling is a BUILD EXPERIMENT, not "
                         "a settled default, because the pangenomes are open (every assembly "
                         "adds k-mer content) but whether that content changes DETECTIONS is "
                         "untested. Compare builds at different ceilings to find out")
    args = ap.parse_args()
    cap = None if str(args.per_species).lower() == "all" else int(args.per_species)
    parent, rank = load_nodes()

    # Annotated in-scope assemblies, grouped by species. The route is per-assembly but
    # constant within a species for the tier decision, so take it from any row.
    by_sp = defaultdict(list)
    for r in csv.DictReader((ln for ln in open(SCREEN) if not ln.startswith("#")),
                            delimiter="\t"):
        if r["annotated"] == "yes":
            by_sp[int(r["species_taxid"])].append(r)

    def tier(route):
        return "db1_congener" if route == "usda_genus" else "db0_baseline"

    rows = []
    for sp, cands in by_sp.items():
        cands.sort(key=lambda r: (LEVEL.get(r["assembly_level"], 0),
                                  -int(r["contig_count"] or 10 ** 9)), reverse=True)
        for r in (cands if cap is None else cands[:cap]):
            routes = {c["route"] for c in cands}
            route = ("usda_species" if "usda_species" in routes else
                     "exception" if "exception" in routes else "usda_genus")
            rows.append((tier(route), r))

    # The host sink. Kept rows only; screen_hosts.tsv also carries every dropped host with
    # its reason, which is the audit trail and not a build input.
    n_host = 0
    if HOSTS.exists():
        for h in csv.DictReader((ln for ln in open(HOSTS) if not ln.startswith("#")),
                                delimiter="\t"):
            if h["decision"] != "keep" or not h["ensembl_species"]:
                continue
            # EVERY cultivar for the analysed host, the reference alone for the others
            # (screen_hosts.py decides which). Taking all 19 wheat cultivars is safe in a
            # way one-species-per-genus is not, for the same reason: they all share taxid
            # 4565, so they consolidate on ONE node and no genus-level LCA bleed is
            # possible, and Kraken2 stores each minimizer once so the redundancy adds no
            # weight. What it buys is k-mer coverage of the plant transcriptome, which is
            # what decides whether a plant read is absorbed or leaks onto a fungus.
            #
            # HOST CDS COMES FROM ENSEMBL PLANTS, NOT NCBI, so the row carries the Ensembl
            # species key as its `cds_ref` and the NCBI GCA only for traceability. NCBI has
            # the wheat pangenome as assemblies without gene models (8 of 12 annotated
            # T. aestivum entries carry 0 protein-coding genes), so it cannot supply this.
            sps = [x for x in h["ensembl_species"].split(";") if x]
            accs = h["assembly_accessions"].split(";")
            for i, sp in enumerate(sps):
                rows.append(("host", {
                    "accession": accs[i] if i < len(accs) else "",
                    "species_taxid": h["ncbi_taxid"],
                    "organism_name": h["host"], "genus": h["genus"], "route": h["role"],
                    "on_wheat": "", "assembly_level": "", "contig_count": "",
                    "genome_size_mb": "", "protein_coding_genes": "",
                    "cds_source": "ensembl_plants", "cds_ref": sp,
                }))
                n_host += 1
    else:
        print("  WARNING: no screen_hosts.tsv; the build plan has NO HOST SINK. "
              "Run screen_hosts.py first.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    # cds_source / cds_ref exist because the host CDS does not come from NCBI. Without them
    # build_plan would describe only the half of the database `datasets` can fetch, and the
    # whole point of carrying hosts here is that ONE table says what was built.
    cols = ["db_tier", "accession", "species_taxid", "organism_name", "genus", "route",
            "on_wheat", "assembly_level", "contig_count", "genome_size_mb",
            "protein_coding_genes", "cds_source", "cds_ref"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for t, r in sorted(rows, key=lambda tr: (tr[0], tr[1]["organism_name"])):
            w.writerow([t, r["accession"], r["species_taxid"], r["organism_name"],
                        r["genus"], r["route"], r["on_wheat"], r["assembly_level"],
                        r["contig_count"], r["genome_size_mb"], r["protein_coding_genes"],
                        r.get("cds_source", "ncbi"), r.get("cds_ref", r["accession"])])

    # report
    n = defaultdict(int)
    sp = defaultdict(set)
    gb = defaultdict(float)
    for t, r in rows:
        n[t] += 1
        sp[t].add(r["species_taxid"])
        gb[t] += (float(r["genome_size_mb"] or 0) * 1e6
                  * CDS_RATIO["host" if t == "host" else "fungus"])
    print(f"build plan: {len(rows):,} assemblies, "
          f"{'all' if cap is None else cap} per species\n")
    if n_host:
        print(f"  {'host sink (in EVERY database)':<32}{n['host']:>5} cds  "
              f"{len(sp['host']):>5} sp  (Ensembl Plants; size from "
              f"raw/host_cds_manifest.tsv)")
    cum_sp = cum_gb = cum_n = 0
    for t, label in (("db0_baseline", "db0 baseline (observations)"),
                     ("db1_congener", "db1 adds congeners")):
        cum_sp += len(sp[t]); cum_gb += gb[t]; cum_n += n[t]
        print(f"  {label:<32}{n[t]:>5} asm  {len(sp[t]):>5} sp  "
              f"+{gb[t] / 1e9:4.1f} GB CDS   (cumulative {cum_n:,} asm, "
              f"~{cum_gb / 1e9:.1f} GB)")
    print(f"\n  db2 adds the 135 gene-predicted genomes in research/predict (worklist.py)")
    print(f"  CDS sized as {CDS_RATIO['fungus']:.2f} x genome for fungi "
          f"(the db_v4 ratio), {CDS_RATIO['host']:.3f} x for hosts (plants are mostly "
          f"repeat: wheat CDS is ~0.2 GB against ~15 GB genomic)")
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
