#!/usr/bin/env python3
"""Screen the unassessed BioProjects for ones that cannot answer the question at all.

PRJEB58762 is why this exists. It sat near the top of the worklist at 473 runs, and it is
building-surface shotgun metagenomics from a UK study: `scientific_name` is `indoor
metagenome`, `library_strategy` is WGS, and one sample is named `Blank4`. It entered through a
Europe PMC accession mention and nothing afterwards asked whether the project was even the
right assay. Reading a paper for it would have been pure waste.

This writes a report. **It does not change runs.tsv and it does not remove anything.** Leon's
standing instruction is that candidates are not silently dropped, and the two screens below
differ sharply in how much they can be trusted, so each project gets a recommendation and a
reason for a human to act on.

## The assay screen is decisive

The cohort is bulk polyA mRNA of plant tissue. A WGS, AMPLICON, Hi-C, WGA, WXS or POOLCLONE
library cannot be that, whatever the organism. This screen is mechanical and safe.

`ncRNA-Seq` and `miRNA-Seq` are separated out rather than lumped in: they are RNA, but
size-selected for small RNA, so they will not carry fungal mRNA at useful depth. Flagged for
review, not rejected, because a few are arms of otherwise usable projects.

## The organism screen is NOT decisive, and must not be applied blindly

`organism` records what the submitter considered the subject, not what is in the tube. An
in-planta rust sample is routinely registered under the rust rather than the wheat, so a
project reading `Puccinia striiformis` may be exactly the dual RNA-seq this study wants, or it
may be axenic spores. The same applies to every fungal and oomycete name here. CLAUDE.md
already records the related trap: `txid{n}[Host]` does not work in SRA esearch because host is
free text.

So fungal, bacterial, viral and metagenome organisms are reported with
`recommend = review`, never `reject`, except for the two classes where the organism alone
settles it: a human sample, and a built-environment metagenome.

Run:  python 01a_Literature/offtarget_screen.py
      python 01a_Literature/offtarget_screen.py --all   (include already-assessed projects)
"""

from __future__ import annotations

import argparse, collections, csv, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUNS = HERE / "data" / "runs.tsv"
BPS = HERE / "data" / "bioprojects.tsv"
OUT = HERE / "data" / "offtarget.tsv"

# Assays that cannot be bulk polyA mRNA of plant tissue.
WRONG_ASSAY = {"WGS", "WGA", "WXS", "AMPLICON", "Hi-C", "POOLCLONE", "CLONE", "ChIP-Seq",
               "Bisulfite-Seq", "ATAC-seq", "RAD-Seq"}
# RNA, but size-selected away from mRNA.
SMALL_RNA = {"ncRNA-Seq", "miRNA-Seq", "ssRNA-seq"}

# Organism classification. `organism` is the submitter's idea of the subject, so this sorts
# projects for review; it decides nothing on its own except the two cases noted in the
# docstring. Matched as lowercase substrings, first hit wins, so order matters.
ORG_CLASS = [
    ("animal",     ["homo sapiens", "mus musculus"]),
    ("built-env",  ["indoor metagenome"]),
    ("metagenome", ["metagenome", "mixed sample", "organismal metagenomes"]),
    ("virus",      ["virus", "viroid", "riboviria", "virga", "alfamovirus", "orthotospovirus",
                    "dichorhavirus", "carlavirus", "torradovirus", "carmovirus", "potyvirus"]),
    ("bacterium",  ["erwinia", "ralstonia", "xanthomonas", "pseudomonas", "agrobacterium",
                    "clavibacter", "streptomyces", "bacillus", "escherichia"]),
    ("oomycete",   ["phytophthora", "plasmopara", "peronospora", "pythium", "bremia",
                    "albugo", "hyaloperonospora"]),
    ("fungus",     ["puccinia", "phakopsora", "fusarium", "colletotrichum", "pyricularia",
                    "magnaporthe", "sclerotinia", "botrytis", "rhizoctonia", "blumeria",
                    "penicillium", "aspergillus", "alternaria", "ustilaginoidea", "epichloe",
                    "moniliophthora", "plenodomus", "pyrenophora", "fulvia", "cytospora",
                    "zymoseptoria", "verticillium", "cercospora", "diaporthe", "monilinia",
                    "venturia", "cladosporium", "septoria", "melampsora", "tilletia",
                    "ustilago", "claviceps", "thielaviopsis", "macrophomina"]),
]


def classify_organism(org: str) -> str:
    o = (org or "").strip().lower()
    if not o:
        return "unknown"
    for label, keys in ORG_CLASS:
        if any(k in o for k in keys):
            return label
    return "plant"      # the default, since the cohort is plant-anchored by construction


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true",
                    help="screen every project, not only the unassessed ones")
    args = ap.parse_args()

    runs = list(csv.DictReader(RUNS.open(), delimiter="\t"))
    bps = {r["BioProject"]: r for r in csv.DictReader(BPS.open(), delimiter="\t")}

    by_bp = collections.defaultdict(list)
    for r in runs:
        if args.all or not r["setting"]:
            by_bp[r["BioProject"]].append(r)

    rows = []
    for bp, rs in by_bp.items():
        strat = collections.Counter(r["LibraryStrategy"] or "(blank)" for r in rs)
        n = len(rs)
        n_rna = strat.get("RNA-Seq", 0)
        n_bad = sum(v for k, v in strat.items() if k in WRONG_ASSAY)
        n_small = sum(v for k, v in strat.items() if k in SMALL_RNA)
        org = (bps.get(bp, {}).get("organism") or "").strip()
        oc = classify_organism(org)

        # The assay test first, because it is the one that decides.
        if n_rna == 0 and n_bad == n:
            rec, why = "reject", f"no RNA-Seq at all; {n} runs are {'/'.join(sorted(k for k in strat if k in WRONG_ASSAY))}"
        elif oc == "built-env":
            rec, why = "reject", "built-environment metagenome, not a plant sample"
        elif oc == "animal":
            rec, why = "reject", f"organism is {org}, not a plant"
        elif n_rna == 0 and n_small == n:
            rec, why = "review", f"all {n} runs are small-RNA ({'/'.join(sorted(k for k in strat if k in SMALL_RNA))}); size-selected away from mRNA"
        elif n_rna < n:
            rec, why = "review", f"mixed assay: {n_rna}/{n} RNA-Seq, rest {'/'.join(f'{k}x{v}' for k, v in strat.most_common() if k != 'RNA-Seq')}"
        elif oc in ("fungus", "oomycete"):
            rec, why = "review", f"organism is the PATHOGEN ({org}); may be in-planta dual RNA-seq or may be axenic culture - the archive cannot tell you which"
        elif oc == "metagenome":
            rec, why = "review", f"organism is '{org}'; check whether the tissue is plant and the library polyA"
        elif oc in ("virus", "bacterium"):
            rec, why = "review", f"organism is a {oc} ({org}); out of scope as a SUBJECT per CLAUDE.md, but the library may still be infected plant tissue"
        else:
            rec, why = "keep", "plant organism, all RNA-Seq"

        rows.append(dict(BioProject=bp, runs=n, rna_seq=n_rna, organism=org,
                         organism_class=oc, assays="; ".join(f"{k}:{v}" for k, v in strat.most_common()),
                         recommend=rec, reason=why,
                         primary_paper=bps.get(bp, {}).get("primary_paper", "")))

    rows.sort(key=lambda r: (r["recommend"] != "reject", -r["runs"]))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    tot = collections.Counter()
    for r in rows:
        tot[r["recommend"]] += 1
        tot[r["recommend"] + "_runs"] += r["runs"]
    print(f"screened {len(rows)} projects, {sum(r['runs'] for r in rows):,} runs\n", file=sys.stderr)
    for rec in ("reject", "review", "keep"):
        print(f"  {rec:<8} {tot[rec]:>4} projects  {tot[rec + '_runs']:>6} runs", file=sys.stderr)
    print(f"\nwrote {OUT.relative_to(ROOT)}", file=sys.stderr)
    print("\nlargest rejects:", file=sys.stderr)
    for r in [x for x in rows if x["recommend"] == "reject"][:8]:
        print(f"  {r['runs']:>5}  {r['BioProject']:<14} {r['reason'][:70]}", file=sys.stderr)


if __name__ == "__main__":
    main()
