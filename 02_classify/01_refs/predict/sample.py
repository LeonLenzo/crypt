#!/usr/bin/env python3
"""sample.py — choose the assemblies for the gene-prediction benchmark.

Step 1 of the predict evaluation. Writes `data/silver/predict_sample.tsv` and nothing else.

## The question this benchmark answers

Vipsania (unsupervised deep ab-initio, Fungi locus-level F1 0.660) could give CDS to the
**2,138 species that have assemblies but no gene model**, which would roughly double db_v4's
species count. Before spending 2,138 genome inferences we need to know whether predicted CDS
is worth having, and `01_refs` already made the opposite call once: `kraken_annotate.py` chose
homology projection (miniprot) over ab initio on the grounds that *"in a competitive index an
invented gene corrupts every species' apportionment, so precision matters more than recall."*
That argument was made without measuring an ab-initio predictor on fungi. This measures it.

## Two arms, measuring two different things, and conflating them would be the mistake

**Arm A — ANNOTATED species. Measures ACCURACY, because ground truth exists.** These species
already have NCBI CDS on disk, so predicting on the same assembly and comparing against its
own `cds_from_genomic.fna` is a clean test. They are already in db_v4, so prediction adds
nothing to the database here; this arm exists purely to calibrate how much to trust a
predictor.

**Arm B — UNANNOTATED species. Measures VALUE ADDED, because no ground truth exists.** Nothing
can be scored for correctness. What can be measured is how much CDS k-mer space a predictor
yields, and **how much of that space is NOVEL against db_v4's existing k-mers.** That novelty
fraction is the number that decides a db_v5: a predicted gene set 95% contained in species
db_v4 already holds adds no competitor value regardless of how accurate it is.

## Why locus-level F1 is the wrong metric for us, and what replaces it

F1 asks "is this gene model correct". A Kraken2 CDS reference does not need correct gene
models; it needs CDS-like sequence from the right organism, because a misplaced exon boundary
still yields real sequence from that species and its k-mers still point at the right taxon.
So the benchmark measures k-mer containment in BOTH directions, which maps onto the two
failure modes that actually cost us:

  * `containment(truth in predicted)` — RECALL of real CDS k-mer space. Too low and the
    species cannot be detected even when present.
  * `containment(predicted in truth)` — PRECISION. Invented sequence inflates the taxon's
    minimizer count, which is the DENOMINATOR every coverage figure in `05_detect` divides by,
    so over-prediction deflates that species' coverage and weakens it as a competitor.

Note this is a narrower claim than `kraken_annotate.py`'s. Invented k-mers do not corrupt
*other* species' apportionment — k-mers shared between taxa resolve to their LCA either way.
What they corrupt is the predicted species' own coverage denominator.

## Stratification, and why assembly level is in it

Contiguity drives ab-initio quality: a gene split across two contigs is a gene the predictor
cannot see whole. The 2,138 unannotated species are mostly fragmented — 1,320 Scaffold, 681
Contig, 75 Chromosome, 62 Complete — so Arm A is stratified on assembly level TOO, over the
same range. Calibrating only on chromosome-level assemblies and applying the result to 681
contig-level ones would be the obvious way to get a confidently wrong answer.

Sampling is seeded and deterministic, so the sample is reproducible from this file alone.

    python 02_classify/01_refs/predict/sample.py --per-arm 60
"""

import argparse
import csv
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

MODULE = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import PREDICT_SAMPLE, REF_CANDIDATES, SILVER  # noqa: E402
from _util import atomic_write  # noqa: E402

ASSEMBLIES = SILVER / "assemblies.tsv"
WORKLIST = SILVER / "miniprot_worklist.tsv"

LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

# Must appear in Arm A whatever the sampler draws. Calibration has to cover the organisms the
# chapter actually makes claims about, and a stratified draw over 1,547 species will usually
# miss all of them. P. nodorum, Z. tritici and the two Pyrenophora are the thesis's wheat and
# barley pathogens; R. solani and F. oxysporum are the open-pangenome extremes (ratio 14.3 and
# 5.2), where a predictor has the most accessory sequence to miss.
FOCAL = [
    "Parastagonospora nodorum", "Zymoseptoria tritici", "Pyrenophora tritici-repentis",
    "Pyrenophora teres", "Puccinia striiformis", "Puccinia triticina",
    "Fusarium graminearum", "Fusarium oxysporum", "Rhizoctonia solani",
    "Pyricularia oryzae", "Bipolaris maydis", "Cercospora beticola",
]

# Residents flagged in CLAUDE.md as present in cereal phyllosphere at real abundance and
# ABSENT from db_v3, i.e. exactly the taxa a db_v5 would exist to add. Over-weighted in Arm B
# so the benchmark speaks to the case that motivates it, not to an average unannotated fungus.
RESIDENT_GENERA = {
    "Vishniacozyma", "Cystofilobasidium", "Nigrospora", "Phoma", "Leptosphaerulina",
    "Acremonium", "Bullera", "Udeniomyces", "Dioszegia", "Cladosporium",
    "Sporobolomyces", "Papiliotrema", "Cryptococcus", "Aureobasidium",
}

# Named ABSENT TARGETS: known plant pathogens with assemblies and ZERO gene models, so they
# are missing from db_v4 entirely. Measured 2026-10-16: **1,159 of 1,974 target species (59%)
# are in this state**, holding 2,093 genomic assemblies. That is the majority of the category
# this chapter is about, which makes Arm B's novelty question larger than "is a db_v5 worth
# it". These are forced into the sample because each already has a documented cost:
#   * F. asiaticum — 257 assemblies, none annotated. In the F. graminearum species complex and
#     a head-blight pathogen, so its reads can only go to F. graminearum.
#   * P. avenae — P. nodorum's ONLY congener. crypt-assembly-saturation records that capping
#     P. nodorum "buys no specificity" because its genus is otherwise empty; the genus is not
#     empty, its sibling just has no gene model.
#   * P. hordei — barley leaf rust, the documented 04_kallisto blind spot, currently
#     misattributed to P. triticina. This is why kraken_annotate.py was written.
ABSENT_TARGETS = [
    "Fusarium asiaticum", "Parastagonospora avenae", "Puccinia hordei",
    "Didymella pinodes", "Bipolaris gigantea", "Penicillium commune",
    "Epichloe typhina", "Microbotryum violaceum",
]


def read_tsv(path: Path) -> list:
    if not path.exists():
        sys.exit(f"missing {path} — run 01_refs/survey.py and select_refs.py first")
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(r, key) -> float:
    try:
        return float(r.get(key) or 0)
    except ValueError:
        return 0.0


def best_assembly(rows: list, fasta_type: str):
    """Most contiguous assembly of the requested type, deterministically."""
    cands = [r for r in rows if r.get("fasta_type") == fasta_type]
    if not cands:
        return None
    return max(cands, key=lambda r: (LEVEL_RANK.get((r.get("assembly_level") or "").strip(), 0),
                                     num(r, "scaffold_n50_kb"), r["accession"]))


def stratified(pool: dict, n: int, rng) -> list:
    """Take n items spread as evenly as possible over the strata, largest-remainder.

    A flat quota per stratum loses whole strata when one is small — the mistake the kallisto
    stratification made, which dropped Zea mays entirely. Here a stratum smaller than its
    quota gives up what it cannot fill and the remainder is redistributed.
    """
    picked, leftovers = [], []
    keys = sorted(pool)
    if not keys:
        return []
    per = max(1, n // len(keys))
    for k in keys:
        items = sorted(pool[k], key=lambda r: r["accession"])
        rng.shuffle(items)
        picked += items[:per]
        leftovers += items[per:]
    rng.shuffle(leftovers)
    return (picked + leftovers)[:n]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--per-arm", type=int, default=60,
                    help="assemblies per arm (default 60). 60 is sized for a containment "
                         "DISTRIBUTION, not a rare event, and keeps the genomic download "
                         "near 7 GB so the benchmark is cheap to redo.")
    ap.add_argument("--seed", type=int, default=20261016,
                    help="sampling seed; the sample is reproducible from it")
    ap.add_argument("--out", default=None, help=f"default {PREDICT_SAMPLE}")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    cands = read_tsv(REF_CANDIDATES)
    asm = read_tsv(ASSEMBLIES)
    work = read_tsv(WORKLIST)

    by_name = defaultdict(list)
    for r in asm:
        by_name[r["organism_name"]].append(r)

    rows = []

    # ── Arm A: annotated, ground truth is its own NCBI CDS ──────────────────────────────────
    # Keyed off ref_candidates rather than assemblies.tsv so every Arm A species is one db_v4
    # actually holds, which is what makes its CDS the truth we compare against.
    ann_by_name = defaultdict(list)
    for r in cands:
        ann_by_name[r["organism_name"]].append(r)

    pool_a = defaultdict(list)
    for name, rs in ann_by_name.items():
        best = best_assembly(rs, "cds")
        if best:
            pool_a[(best.get("assembly_level") or "").strip() or "unknown"].append(best)

    chosen_a = {}
    for name in FOCAL:                                  # forced first, before any quota
        best = best_assembly(ann_by_name.get(name, []), "cds")
        if best:
            chosen_a[best["accession"]] = best
    for r in stratified(pool_a, max(0, args.per_arm - len(chosen_a)), rng):
        chosen_a.setdefault(r["accession"], r)
    chosen_a = dict(list(chosen_a.items())[:max(args.per_arm, len(FOCAL))])

    for acc, r in chosen_a.items():
        rows.append({
            "arm": "annotated", "accession": acc, "taxid": r["taxid"],
            "organism_name": r["organism_name"], "routes": r.get("routes", ""),
            "assembly_level": r.get("assembly_level", ""),
            "scaffold_n50_kb": r.get("scaffold_n50_kb", ""),
            "protein_coding_genes": r.get("protein_coding_genes", ""),
            "has_truth": "yes",
            "stratum": (r.get("assembly_level") or "").strip(),
            "focal": "yes" if r["organism_name"] in FOCAL else "no",
        })

    # ── Arm B: unannotated, no ground truth; measures novelty against db_v4 ─────────────────
    work_names = {r["name"] for r in work}
    routes_of = {r["name"]: r.get("routes", "") for r in work}
    pool_b = defaultdict(list)
    resident_hits = []
    for name in sorted(work_names):
        best = best_assembly(by_name.get(name, []), "genome")
        if not best:
            continue
        best = dict(best)
        best["routes"] = routes_of.get(name, best.get("routes", ""))
        lvl = (best.get("assembly_level") or "").strip() or "unknown"
        if name.split()[0] in RESIDENT_GENERA:
            resident_hits.append(best)
        else:
            pool_b[lvl].append(best)

    chosen_b = {}
    for name in ABSENT_TARGETS:                          # forced, before any quota
        best = best_assembly(by_name.get(name, []), "genome")
        if best:
            best = dict(best)
            best["routes"] = routes_of.get(name, best.get("routes", ""))
            chosen_b[best["accession"]] = best
    rng.shuffle(resident_hits)
    for r in resident_hits[:args.per_arm // 3]:          # ~a third reserved for residents
        chosen_b[r["accession"]] = r
    for r in stratified(pool_b, max(0, args.per_arm - len(chosen_b)), rng):
        chosen_b.setdefault(r["accession"], r)
    chosen_b = dict(list(chosen_b.items())[:args.per_arm])

    for acc, r in chosen_b.items():
        rows.append({
            "arm": "unannotated", "accession": acc, "taxid": r["taxid"],
            "organism_name": r["organism_name"], "routes": r.get("routes", ""),
            "assembly_level": r.get("assembly_level", ""),
            "scaffold_n50_kb": r.get("scaffold_n50_kb", ""),
            "protein_coding_genes": "",
            "has_truth": "no",
            "stratum": (r.get("assembly_level") or "").strip(),
            "focal": ("absent_target" if r["organism_name"] in ABSENT_TARGETS
                      else "resident" if r["organism_name"].split()[0] in RESIDENT_GENERA
                      else "no"),
        })

    out = Path(args.out) if args.out else PREDICT_SAMPLE
    cols = ["arm", "accession", "taxid", "organism_name", "routes", "assembly_level",
            "scaffold_n50_kb", "protein_coding_genes", "has_truth", "stratum", "focal"]
    with atomic_write(out) as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"sample → {out}")
    for arm in ("annotated", "unannotated"):
        sub = [r for r in rows if r["arm"] == arm]
        print(f"\n  {arm}: {len(sub)} assemblies, "
              f"{len({r['organism_name'] for r in sub})} species")
        for lvl, n in sorted(Counter(r["stratum"] for r in sub).items(),
                             key=lambda kv: -LEVEL_RANK.get(kv[0], 0)):
            print(f"      {lvl or 'unknown':16s} {n:3d}")
        for tag in ("yes", "absent_target", "resident"):
            k = len([r for r in sub if r["focal"] == tag])
            if k:
                label = "(focal)" if tag == "yes" else f"({tag})"
                print(f"      {label:16s} {k:3d}")
    print(f"\nseed {args.seed}; rerunning with the same seed reproduces this sample exactly.")
    print("Next: 01_refs/predict/fetch_genomes.py  (GENOMIC fasta; db_v4 holds only CDS)")


if __name__ == "__main__":
    main()
