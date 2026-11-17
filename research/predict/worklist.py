#!/usr/bin/env python3
"""worklist.py — which assemblies to gene-predict, and why each one is on the list.

    python research/predict/worklist.py
    python research/predict/worklist.py --tiers A,B        # drop the competitor tier
    python research/predict/worklist.py --all-dark         # tier C: every dark genus

Reads  release/02_reference/data/outputs/supp_data_s3_reference.tsv   the published scope
       release/02_reference/data/curation/screen_assemblies.tsv       to pick the accession
Writes research/predict/data/worklist.tsv

A SELECTION, NEVER A COPIED FILE, so it stays correct when `assembly_summary` refreshes.
S3 decides who is on the list; screen_assemblies.tsv only supplies which accession to send,
because S3 is one row per taxon and carries no accessions.

ONE TOOL: VIPSANIA, FOR EVERY TARGET. miniprot was kept for targets with an annotated
congener, and that split is dropped deliberately. The reason is comparability rather than
convenience: miniprot biases toward genes conserved with its guide and vipsania over-calls
in repeats, so a database built from both has two bias profiles and the minimizer
denominators of two taxa differ for reasons that are not biological. Coverage is already
not comparable ACROSS database versions; mixing methods would make it not comparable
BETWEEN TAXA WITHIN one, with nothing to flag it. The benchmark backs the choice: vipsania
scores recall 0.963 / precision 0.945 by k-mer containment over 60 annotated assemblies,
which is the metric a classifier reference is judged on (locus F1 0.660 answers a question
nobody here is asking).

## The tiers, and why these and not the 932

The old worklist was "species with no gene model whose genus has none either" = 932 species.
It is the wrong target for two reasons. The unit is the GENUS, since one prediction gives
every species in it a k-mer neighbour, which turns 932 species into 437 genera. And ranked
by species count those genera are Leccinum, Inocybe, Xerocomus and lichens: ectomycorrhizal
and lichenised taxa that will not appear in a wheat leaf.

    A  SILENT REPORTABLE   in the wheat row or a pathogen column, and no gene model, so it
                           is in the reference and cannot be called. Led by Fusarium
                           asiaticum, a wheat pathogen with 259 assemblies and zero
                           annotated, and Puccinia hordei, the documented case of a read
                           being misattributed to P. triticina.
    B  ON-WHEAT DARK GENUS one assembly for each genus wholly absent from the database that
                           USDA records on wheat. 13 of the 437 dark genera qualify.
    C2 CONGENER COMPETITOR a callable wheat-relevant taxon with no ANNOTATED congener has no
                           within-genus competitor, so Kraken2's LCA can never decline and
                           it gets a confident call whether or not it is present. That is
                           the db_v2 failure that put soybean rust in every wheat sample.
                           Where an unannotated congener exists, predicting it supplies the
                           competitor. Parastagonospora nodorum is the worst case: 155
                           annotated assemblies and nothing to compete with them, fixed by
                           predicting P. avenae.

    --all-dark adds tier C, one assembly per dark genus for all 437. Do this only after the
    purpose test below passes on A+B+C2.

## What cannot be fixed, and must be stated rather than solved

A genus with no congener anywhere in GenBank cannot be given one. Oculimacula (1 assembly
in total) and Rhizoctonia (51, all one species) are in that position, and both already have
an annotated FAMILY neighbour in scope (Cadophora, Ceratobasidium), which is weaker
competition because only the more conserved k-mers are shared, but is real. Phakopsora
pachyrhizi and Oidiodendron maius have no neighbour at genus OR family level, so no
prediction helps them and their calls rest on read-level validation alone. Across the whole
reference 259 callable taxa are in that position.

## How a prediction gets checked

predict_counts compares predicted CDS with the true annotation by FracMinHash containment.
It CANNOT be used here: the worklist is by definition the set with no annotation, so there
is no truth to compare against. Annotation-free checks instead, in order of strength:

    1. BUSCO on the predicted proteins. A real fungal gene set recovers most of its
       lineage's complement, and this needs no reference.
    2. Cross-assembly reproducibility, where a species has more than one assembly: predict
       both independently and measure containment between the two predictions.
    3. `predict_risk` below, as a PRIOR rather than a check. Benchmark deviation tracks
       genome size at r = -0.53 and Blumeria graminis over-called 3.9x at precision 0.513,
       so repeat-rich targets inflate the minimizer denominator and their coverage is not
       comparable with a downloaded taxon's.
    4. The purpose test, which is the only one that matters: classify a sample set before
       and after adding the predicted CDS and confirm reads move DOWN, from genus and order
       nodes to species. 31% of fungal reads currently resolve no finer than genus.

Predicted CDS stays COMPETITOR-ONLY and never supports a positive call about the species it
came from. That policy is unchanged by the tool choice.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent                      # phd/01-review
REF = REPO / "release/02_reference"
S3 = REF / "data/outputs/supp_data_s3_reference.tsv"
ASM = REF / "data/curation/screen_assemblies.tsv"
OUT = HERE / "data/worklist.tsv"

# mirrors plan_build.py: better level first, then fewest contigs
LEVEL = {"Complete Genome": 3, "Chromosome": 2, "Scaffold": 1, "Contig": 0}
# ...but SIZE OUTRANKS LEVEL, because NCBI's assembly_level is a submitter label and nothing
# validates it. GCA_986280995.1 is 528 bp on one contig, registered as "Complete Genome",
# and a level-first sort put it ahead of a real 40 Mb scaffold assembly for Nigrospora
# sphaerica. Sending that to a gene finder returns nothing and the task looks like a crash.
# The floor is deliberately low: the smallest real fungal genomes (microsporidia) are ~2.3 Mb.
MIN_GENOME_MB = 5.0
SU_PER_GENOME = 12.0   # MEASURED on GPU, pilot 53221324 (120 SU / 10 genomes). The 51.1 in
                       # the README is the CPU route, which the GPU trial superseded.
BIG_MB, MANY_CONTIGS = 100.0, 5000             # repeat/fragmentation risk flags
csv.field_size_limit(10 ** 7)


def rd(p):
    with open(p) as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                   delimiter="\t"))


def best_assembly(rows):
    """The one assembly to send: plausible size first, then level, then fewest contigs.

    An undersized record is ranked last rather than dropped, so a species whose only
    assembly is a fragment still appears and is visible as such in `predict_risk`.
    """
    return sorted(rows, key=lambda r: (float(r["genome_size_mb"] or 0) >= MIN_GENOME_MB,
                                       LEVEL.get(r["assembly_level"], 0),
                                       -int(r["contig_count"] or 10 ** 9)),
                  reverse=True)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tiers", default="A,B,C2", help="comma-separated subset to emit")
    ap.add_argument("--all-dark", action="store_true",
                    help="tier C: one assembly for EVERY dark genus, not just on-wheat")
    args = ap.parse_args()
    want = {t.strip() for t in args.tiers.split(",") if t.strip()}
    if args.all_dark:
        want.add("C")

    ref = rd(S3)
    for r in ref:
        r["n_ann"] = int(r["n_annotated"])
        r["n_asm"] = int(r["n_assemblies"])
    by_gen, gen_ann = defaultdict(list), defaultdict(int)
    for r in ref:
        by_gen[r["genus"]].append(r)
        gen_ann[r["genus"]] += r["n_ann"]
    reportable = lambda r: r["usda_row"] == "wheat" or r["phibase_col"] != "none"

    picked = {}                                 # taxid -> (tier, reason)

    def add(r, tier, reason):
        if r["ncbi_taxid"] not in picked:
            picked[r["ncbi_taxid"]] = (tier, reason)

    # A: in the reference, reportable or pathogen-labelled, and silent
    if "A" in want:
        for r in ref:
            if r["n_ann"] == 0 and reportable(r):
                add(r, "A", f"silent {r['usda_row']} x {r['phibase_col']}, "
                            f"{r['n_asm']} assemblies, no gene model")

    # B / C: one assembly per dark genus
    dark = [g for g in by_gen if g and gen_ann[g] == 0]
    for g in dark:
        on_wheat = any(x["usda_row"] == "wheat" for x in by_gen[g])
        if on_wheat and "B" in want:
            add(max(by_gen[g], key=lambda r: r["n_asm"]), "B",
                f"dark genus recorded on wheat ({len(by_gen[g])} species in scope)")
        elif "C" in want:
            add(max(by_gen[g], key=lambda r: r["n_asm"]), "C",
                f"dark genus ({len(by_gen[g])} species in scope)")

    # C2: a competitor for every callable wheat-relevant taxon that has no annotated congener
    if "C2" in want:
        for r in ref:
            if not r["n_ann"] or not reportable(r):
                continue
            sibs = [s for s in by_gen[r["genus"]] if s["ncbi_taxid"] != r["ncbi_taxid"]]
            if not sibs or any(s["n_ann"] > 0 for s in sibs):
                continue
            add(max(sibs, key=lambda s: s["n_asm"]), "C2",
                f"competitor for {r['taxon']} ({r['n_ann']} annotated assemblies, "
                f"no annotated congener)")

    asm_by_sp = defaultdict(list)
    for a in rd(ASM):
        asm_by_sp[a["species_taxid"]].append(a)

    cols = ["tier", "ncbi_taxid", "taxon", "genus", "family", "usda_row", "phibase_col",
            "n_assemblies", "accession", "assembly_level", "contig_count", "genome_size_mb",
            "predict_risk", "reason"]
    info = {r["ncbi_taxid"]: r for r in ref}
    rows = []
    for tid, (tier, reason) in picked.items():
        r = info[tid]
        cands = asm_by_sp.get(tid, [])
        a = best_assembly(cands) if cands else {}
        mb = float(a.get("genome_size_mb") or 0)
        nc = int(a.get("contig_count") or 0)
        flags = []
        if mb and mb < MIN_GENOME_MB:
            flags.append("UNDERSIZED")      # no real genome available for this taxon
        if mb >= BIG_MB:
            flags.append("large_genome")
        if nc >= MANY_CONTIGS:
            flags.append("fragmented")
        risk = ";".join(flags)
        rows.append({
            "tier": tier, "ncbi_taxid": tid, "taxon": r["taxon"], "genus": r["genus"],
            "family": r["family"], "usda_row": r["usda_row"],
            "phibase_col": r["phibase_col"], "n_assemblies": r["n_asm"],
            "accession": a.get("accession", ""), "assembly_level": a.get("assembly_level", ""),
            "contig_count": nc or "", "genome_size_mb": a.get("genome_size_mb", ""),
            "predict_risk": risk, "reason": reason})
    rows.sort(key=lambda x: (x["tier"], x["taxon"]))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    per = defaultdict(int)
    for x in rows:
        per[x["tier"]] += 1
    print(f"{len(rows)} genomes, {len(rows) * SU_PER_GENOME:,.0f} SU at "
          f"{SU_PER_GENOME} SU/genome")
    for t in sorted(per):
        print(f"  tier {t:<3}{per[t]:>5}")
    flagged = [x for x in rows if x["predict_risk"]]
    print(f"  flagged repeat/fragmentation risk: {len(flagged)}")
    missing = [x for x in rows if not x["accession"]]
    if missing:
        print(f"  WARNING no assembly resolved for {len(missing)}: "
              + ", ".join(x["taxon"] for x in missing[:5]))
    print(f"\n  -> {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
