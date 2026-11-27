#!/usr/bin/env python3
"""Decide which taxa and which assembly go into the alignment reference.

    python research/03_align/01_index/plan_index.py --limit-family Pleosporaceae
    python research/03_align/01_index/plan_index.py --scope genus
    python research/03_align/01_index/plan_index.py --scope family --check-local

Reads  release/03_classify/data/outputs/supp_data_s5_detections.tsv   what db2 called
       release/02_reference/data/outputs/supp_data_s3_reference.tsv   taxon -> genus, family
       release/02_reference/data/curation/screen_assemblies.tsv       every annotated assembly
       release/02_reference/data/curation/build_plan.tsv              the assembly ALREADY chosen
       release/02_reference/data/curation/extra_cds_db2.tsv           CDS we predicted
       release/02_reference/data/raw/host_cds_manifest.tsv            the host CDS file
Writes data/index_plan_<label>.tsv    one row per taxon in the index, with the reason
       data/cds_files_<label>.txt     the CDS paths, one per line, for the fetch step

Research, not release: this runs to find out what the data looks like. The clean version is
written once the first index has been aligned against.

## WHO IS IN THE INDEX

Three roles, and the distinction is load-bearing rather than cosmetic:

    confirm      db2 called it. This is what job 1 is testing.
    competitor   in scope, not called. It is here so that `species-unique` means unique
                 against the plausible alternatives rather than unique by their absence.
    host         Chinese Spring, as the read sink.

Scope is the competitor question and the only real choice. `--scope genus` takes every
reference taxon in a genus holding a call; `--scope family` widens it to the family.
`--limit-family` restricts either to one clade, which is how the first index is built.

**The k-mer measurement says genus and does not get the last word.** On the crossover matrix,
no cross-genus pair of the 254 wheat-row taxa shares more than 4% of the smaller taxon's
minimizers and only 4 pairs of ~10,800 pass 1%. But that is exact k=35 sharing and bowtie2
scores mismatches, and the artefact this module exists to adjudicate was cross-genus inside a
family: Pleosporaceae reads landing on *Stemphylium lycopersici* for want of a nearer
reference. So the exact-k-mer number is a floor on how much competitors matter, not a ceiling,
and the boundary gets MEASURED by quantifying the same reads against the genus-only subset of
a family-scoped index. Hence both scopes come out of one script.

## WHICH ASSEMBLY, AND WHY THE CHOICE IS OURS

**`build_plan.tsv` is not a per-species pick, which is the first thing checked here rather
than assumed.** It carries all 5,302 annotated assemblies in scope, 34 of them for
*A. alternata* alone, because the classifier wanted every isolate's k-mer space. So it is used
as a MEMBERSHIP FILTER, not as a choice: an assembly must be one the cohort was classified
against, since aligning reads to sequence the classifier never held would confirm a detection
against something that could not have produced it.

**One assembly per species, for a different reason than the classifier had.** Alignment needs
one gene model per locus: near-duplicate transcripts do not hurt species-level uniqueness,
because the tags are species, but they split a gene's reads across redundant copies and inflate
gene breadth, which is one of the three presence facts job 1 reads.

The rule is assembly level first, as `plan_build.py` has it, because a fragmented assembly
fragments its gene models. Then **nearest the species median gene count**, and that is the part
that was measured rather than chosen on principle:

- *A. alternata*: level-then-fewest-contigs picks GCA_016097525.1, a 24,347-gene annotation,
  over a 13,466-gene one on a ONE CONTIG difference. The species' 34 annotations sit at a
  median of 11,928, so the fewest-contigs rule would have put two copies of most loci in the
  index, in the one genus the whole adjudication is about.
- *A. melanogenum*: its 48 annotations are BIMODAL, 7,480 to 10,428 genes at ~26 Mb and 17,274
  to 27,149 at ~49 Mb. A median-band reject does not help here, it just picks whichever cluster
  is more numerous. Proximity to the median does not pretend to know which cluster is the real
  organism; it goes with the bulk of the evidence and FLAGS the species as bimodal.

**Two different tests, and conflating them is what the first pass got wrong.** A HARD FLOOR
rejects an assembly with no gene model at all: *A. clavatophorus*' only annotated assembly
declares five protein-coding genes, and a rule that merely flagged outliers would have admitted
it as a competitor that absorbs nothing. A FLAG, not a rejection, records a gene model that
disagrees with its relatives, because metadata cannot tell a hybrid from a duplicated
annotation from a genuinely larger genome. **The sequence decides that, not the table**: the
flagged picks get a within-species redundancy check once the CDS is on disk, which is a
measurement this script cannot make.

A called taxon is never dropped. If its only annotation is below the floor it goes in with
`no_usable_annotation` set, because testing it is the job and a silent omission would read as
an absence.

## METADATA PROPOSES, SEQUENCE DISPOSES, AND IT IS A LOOP

`--sequence-check` reads `redundancy.py`'s verdicts and rejects a pick the SEQUENCE fails,
taking the next by the same ranking rule. Two stages rather than one because this script is
metadata-only on purpose: it reads four small tables, runs in a second, and is safe on a login
node, while measuring 1.2 Gbp of CDS is a job.

So the order is: plan, measure, re-plan with `--sequence-check`, measure the new picks, and
stop when no pick moves. A substituted pick has never been measured, which is why the run
prints how many picks are still unmeasured instead of implying the gate has seen everything.
Measured 2026-11-26, the first pass moved nothing in the Pleosporaceae plan, because proximity
to the species median had already declined the one over-annotated assembly for its own reasons.

**What the sequence gate catches that the metadata cannot** is in `redundancy.py`'s docstring,
and it is not what this script originally guessed: *A. alternata*'s 24,347-gene annotation
holds no duplicated loci at all (0.002), it holds 24,347 SHORT ones, median 576 bp against
1,149. The flag was right and the stated reason was wrong.

Predicted CDS keeps its provenance in `cds_source`. For a taxon with no NCBI annotation it is
the only CDS there is, and it enters as a competitor; a POSITIVE resting on it is a presence
claim about our own gene prediction and has to be labelled that way wherever it is quoted.
"""

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common as C  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "data"
LEVEL_ORDER = ("Complete Genome", "Chromosome", "Scaffold", "Contig")
HARD_MIN_GENES = 500             # below this an "annotated" assembly has no gene model at all
HARD_MIN_FRACTION = 0.1          # ... or less than a tenth of what its relatives carry
OUTLIER_FRACTION = 0.5           # |genes - median| / median past which a pick is FLAGGED
BIMODAL_RATIO = 1.8              # max/min gene count at which a species' annotations disagree
MIN_FOR_SPECIES_BAND = 3         # annotated assemblies needed before a species sets its own median
COLS = ("taxid", "taxon", "genus", "family", "role", "is_wheat_pathogen", "db2_calls",
        "cds_source", "accession", "cds_path", "assembly_level", "contig_count",
        "protein_coding_genes", "median_genes", "median_from", "flags", "pick_reason",
        "runners_up")


def rd(path):
    with open(path, newline="", encoding="utf-8") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def genes(row):
    """Protein-coding gene count, or None when NCBI gives no usable number."""
    value = (row.get("protein_coding_genes") or "").strip()
    return int(value) if value.isdigit() and int(value) > 0 else None


def rank(row, median=None):
    """Sort key: assembly level, then NEAREST THE SPECIES MEDIAN gene count, then contiguity.

    Level first, as `plan_build.py` has it, because a fragmented assembly fragments its gene
    models. Then proximity to the median rather than fewest contigs, which is the correction
    that matters: on *A. alternata* the fewest-contigs rule chose a 24,347-gene annotation over
    a 13,466-gene one on a one-contig difference, and a doubled gene set puts two copies of
    every locus in the index. Proximity is robust to an outlier in either direction and does
    not require deciding which of a bimodal species is correct. Contigs, date and accession
    break the remaining ties; the accession is last so the result is deterministic.
    """
    level = LEVEL_ORDER.index(row["assembly_level"]) if row["assembly_level"] in LEVEL_ORDER else 9
    contigs = int(row["contig_count"]) if row["contig_count"].isdigit() else 10 ** 9
    distance = abs((genes(row) or 0) - median) if median else 0
    newest_first = tuple(-int(part) for part in row["seq_rel_date"].split("-")
                         if part.isdigit()) or (0,)
    return (level, distance, contigs, newest_first, row["accession"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scope", choices=("genus", "family"), default="genus",
                    help="how far out the competitors go (default: genus)")
    ap.add_argument("--limit-family", metavar="NAME",
                    help="restrict to one family, e.g. Pleosporaceae for the first index")
    ap.add_argument("--check-local", action="store_true",
                    help="report how many CDS files are present on this machine")
    ap.add_argument("--sequence-check", metavar="TSV", nargs="?", const="auto",
                    help="apply redundancy.py's verdicts: reject a pick the SEQUENCE fails "
                         "and take the next by the same rule (default file: "
                         "data/redundancy_<label>.tsv)")
    args = ap.parse_args()

    # The sequence verdicts, when there are any. METADATA PROPOSES, SEQUENCE DISPOSES: this
    # file is written by redundancy.py against the CDS itself, and it is the only thing here
    # that can tell an over-annotation from a larger genome. It is optional because the first
    # pass has nothing to read: plan, measure, re-plan, measure again until no pick moves.
    verdicts = {}
    if args.sequence_check:
        label = (args.limit_family or args.scope).lower()
        path = (PLANS / f"redundancy_{label}.tsv" if args.sequence_check == "auto"
                else Path(args.sequence_check))
        if not path.exists():
            sys.exit(f"no sequence measurements at {path}; run redundancy.py --label {label} "
                     f"first, or drop --sequence-check for the metadata-only pass")
        verdicts = {r["accession"]: r for r in rd(path)}

    calls = defaultdict(int)
    wheat_pathogen = {}
    for r in rd(C.DETECTIONS):
        if r["arm"] != C.ARM or r["taxon_class"] == "host":
            continue
        calls[r["taxid"]] += 1
        wheat_pathogen[r["taxid"]] = r["is_wheat_pathogen"]

    reference = {r["ncbi_taxid"]: r for r in rd(C.REFERENCE)}
    missing = [t for t in calls if t not in reference]
    if missing:
        sys.exit(f"{len(missing)} called taxa are absent from the reference table: {missing[:5]}")

    genera = {reference[t]["genus"] for t in calls if reference[t]["genus"]}
    families = {reference[t]["family"] for t in calls if reference[t]["family"]}
    unplaced = {reference[t]["genus"] for t in calls if not reference[t]["family"]}

    # TWO CALLED TAXA HAVE NO FAMILY IN NCBI, Tilletiopsis washingtonensis and Trichothecium
    # roseum, and neither genus is placed either, so there is nothing to infer from their
    # congeners. At family scope their competitor set falls back one rank to their own genus
    # rather than being dropped: a scope rule that silently tested 88 of 90 would report the
    # scope as the reason for a gap it created.
    def in_scope(ref):
        if args.scope == "genus":
            return ref["genus"] in genera
        return (ref["family"] in families if ref["family"] else False) or ref["genus"] in unplaced

    selected = [r for r in reference.values() if in_scope(r)]
    if args.limit_family:
        selected = [r for r in selected if r["family"] == args.limit_family]
        if not selected:
            sys.exit(f"no taxon in scope for family {args.limit_family!r}")

    annotated = defaultdict(list)
    genus_genes = defaultdict(list)
    species_genes = defaultdict(list)
    for a in rd(C.ASSEMBLIES):
        if a["annotated"] != "yes":
            continue
        annotated[a["species_taxid"]].append(a)
        if genes(a):
            genus_genes[a["genus"]].append(genes(a))
            species_genes[a["species_taxid"]].append(genes(a))
    genus_median = {g: statistics.median(v) for g, v in genus_genes.items()}
    species_median = {s: statistics.median(v) for s, v in species_genes.items()
                      if len(v) >= MIN_FOR_SPECIES_BAND}

    # build_plan is the set of assemblies the BUILD used, not a per-species pick: it carries
    # all 5,302 annotated assemblies, 34 of them for A. alternata alone. So it is a membership
    # filter here, never a choice, and the choice below is this module's own.
    in_build = defaultdict(set)
    for p in rd(C.BUILD_PLAN):
        if p["db_tier"].startswith("db") and p["accession"]:
            in_build[p["species_taxid"]].add(p["accession"])

    predicted = {}
    for p in rd(C.EXTRA_CDS):
        predicted.setdefault(p["taxid"], p["path"])

    rows, dropped, flagged = [], [], []
    for ref in sorted(selected, key=lambda r: (r["family"], r["genus"], r["taxon"])):
        taxid = ref["ncbi_taxid"]
        role = "confirm" if taxid in calls else "competitor"
        counts = sorted(species_genes.get(taxid, []))
        median = (statistics.median(counts) if len(counts) >= MIN_FOR_SPECIES_BAND
                  else genus_median.get(ref["genus"]))
        band_from = ("species" if len(counts) >= MIN_FOR_SPECIES_BAND
                     else "genus" if median else "none")
        bimodal = bool(len(counts) >= MIN_FOR_SPECIES_BAND
                       and counts[-1] / max(counts[0], 1) >= BIMODAL_RATIO)

        # Eligible = annotated, carrying a real gene model, and in the build the cohort was
        # classified against. The last condition matters: aligning against sequence the
        # classifier never held would confirm a detection against something that could not
        # have produced it.
        #
        # THE HARD FLOOR IS NOT THE OUTLIER FLAG. A. clavatophorus' only annotated assembly
        # declares FIVE protein-coding genes; a rule that only flags outliers would have put
        # it in the index as a competitor that absorbs nothing and reports as a reference.
        # Below the floor there is no gene model, which is a different statement from a gene
        # model that disagrees with its relatives.
        hard_floor = max(HARD_MIN_GENES, (median or 0) * HARD_MIN_FRACTION)
        usable = [a for a in annotated.get(taxid, []) if (genes(a) or 0) >= hard_floor]
        unusable = [a for a in annotated.get(taxid, []) if a not in usable]
        candidates = usable
        in_this_build = [a for a in candidates if a["accession"] in in_build.get(taxid, ())]
        candidates = in_this_build or candidates

        def sequence_ok(a):
            row = verdicts.get(a["accession"])
            return row is None or row["verdict"] == "ok"

        pick, reason = None, ""
        if not candidates and role == "confirm" and unusable:
            # A called taxon is what the job tests, so it never silently vanishes; it goes in
            # with the problem named instead.
            candidates = unusable
        rejected_by_sequence = []
        if candidates:
            ordered = sorted(candidates, key=lambda a: rank(a, median))
            passing = [a for a in ordered if sequence_ok(a)]
            rejected_by_sequence = [a for a in ordered[:ordered.index(passing[0])]
                                    if not sequence_ok(a)] if passing else []
            pick = (passing or ordered)[0]
            reason = (f"best of {len(candidates)} by level, then nearest the {band_from} "
                      f"median of {median:.0f} genes" if median and len(candidates) > 1
                      else "the only annotated assembly with a gene model")
            for a in rejected_by_sequence:
                reason += (f"; the sequence gate rejected {a['accession']} "
                           f"({verdicts[a['accession']]['verdict']})")
            if not passing and any(a["accession"] in verdicts for a in ordered):
                reason += "; EVERY measured candidate fails the sequence gate"
            if not in_this_build:
                reason += "; NOT in the build the cohort was classified against"
            if (genes(pick) or 0) < hard_floor:
                reason += (f"; the only annotation available declares {genes(pick)} genes, "
                           f"below the floor of {hard_floor:.0f}")

        off = (median and pick and abs(genes(pick) - median) / median > OUTLIER_FRACTION)
        unmeasured = bool(pick and verdicts and pick["accession"] not in verdicts)
        flags = ";".join(f for f, on in (
            ("no_usable_annotation", bool(pick and (genes(pick) or 0) < hard_floor)),
            ("sequence_substituted", bool(rejected_by_sequence)),
            ("sequence_unmeasured", unmeasured),
            ("gene_count_outlier", bool(off)),
            ("species_bimodal", bimodal),
            ("not_in_build", bool(candidates and not in_this_build))) if on)
        if flags and pick:
            flagged.append((ref["taxon"], pick["accession"], genes(pick), median, flags))

        if pick:
            source, accession = "ncbi", pick["accession"]
            path = f"{C.CDS_DIR.relative_to(C.ROOT)}/{accession}/cds_from_genomic.fna"
            level, contigs, n_genes = pick["assembly_level"], pick["contig_count"], genes(pick)
        elif taxid in predicted:
            source, accession = "predicted", ""
            path, level, contigs, n_genes = predicted[taxid], "", "", None
            reason = "no NCBI annotation; CDS predicted in research/predict"
        else:
            dropped.append(ref["taxon"])
            continue

        others = [a for a in candidates if a["accession"] != accession][:3]
        rows.append({
            "taxid": taxid, "taxon": ref["taxon"], "genus": ref["genus"],
            "family": ref["family"], "role": role,
            "is_wheat_pathogen": ref["is_wheat_pathogen"],
            "db2_calls": calls.get(taxid, 0), "cds_source": source, "accession": accession,
            "cds_path": path, "assembly_level": level, "contig_count": contigs,
            "protein_coding_genes": n_genes if n_genes is not None else "",
            "median_genes": f"{median:.0f}" if median else "",
            "median_from": band_from, "flags": flags,
            "pick_reason": reason,
            "runners_up": ";".join(f"{a['accession']}({a['assembly_level']},"
                                   f"{genes(a) or '?'}g)" for a in others),
        })

    host = [h for h in rd(C.HOST_CDS_MANIFEST) if h["ensembl_species"] == C.HOST_ENSEMBL]
    if len(host) != 1:
        sys.exit(f"expected one host row for {C.HOST_ENSEMBL}, found {len(host)}")
    rows.append({
        "taxid": C.HOST_TAXID, "taxon": "Triticum aestivum", "genus": "Triticum",
        "family": "Poaceae", "role": "host", "is_wheat_pathogen": "no", "db2_calls": 0,
        "cds_source": "ensembl_plants", "accession": C.HOST_ACCESSION,
        "cds_path": f"{C.HOST_CDS_DIR.relative_to(C.ROOT)}/{C.HOST_ENSEMBL}/{host[0]['file']}",
        "assembly_level": "Chromosome", "contig_count": "", "protein_coding_genes": "",
        "median_genes": "", "median_from": "", "flags": "", "pick_reason": "the sink; Chinese Spring alone, decided "
                                                 "2026-11-25, NOT db2's 19 cultivars",
        "runners_up": "",
    })

    label = (args.limit_family or args.scope).lower()
    OUT.mkdir(parents=True, exist_ok=True)
    plan_path = OUT / f"index_plan_{label}.tsv"
    with open(plan_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    list_path = OUT / f"cds_files_{label}.txt"
    list_path.write_text("".join(f"{r['cds_path']}\n" for r in rows), encoding="utf-8")

    fungal = [r for r in rows if r["role"] != "host"]
    confirm = [r for r in fungal if r["role"] == "confirm"]
    pred = [r for r in fungal if r["cds_source"] == "predicted"]
    volume = sum(int(r["protein_coding_genes"]) if r["protein_coding_genes"] else 11_000
                 for r in fungal) * 1.6e3
    print(f"scope {args.scope}" + (f", family {args.limit_family}" if args.limit_family else ""))
    print(f"  {len(fungal)} fungal taxa: {len(confirm)} confirm, "
          f"{len(fungal) - len(confirm)} competitor")
    print(f"  CDS: {len(fungal) - len(pred)} NCBI-annotated, {len(pred)} predicted here")
    print(f"  wheat pathogens among the confirm set: "
          f"{sum(1 for r in confirm if r['is_wheat_pathogen'] == 'yes')}")
    print(f"  ~{volume / 1e9:.2f} Gbp of CDS, before the host")
    print(f"  {len(dropped)} in-scope taxa have no CDS at all and cannot compete")
    if flagged:
        print(f"  {len(flagged)} picks FLAGGED for the sequence-level check, not rejected here:")
        for taxon, accession, n, median, flags in flagged[:8]:
            print(f"      {taxon}: {accession} at {n}g against a median of {median:.0f} [{flags}]")
        if len(flagged) > 8:
            print(f"      ... {len(flagged) - 8} more, all in the plan's flags column")
    else:
        print("  no pick flagged: measured on this scope, not assumed")
    if verdicts:
        moved = [r for r in rows if "sequence_substituted" in r["flags"]]
        unknown = [r for r in rows if "sequence_unmeasured" in r["flags"]]
        print(f"  sequence gate: {len(verdicts)} assemblies measured, "
              f"{len(moved)} picks moved, {len(unknown)} picks not yet measured")
        for r in moved:
            print(f"      {r['taxon']}: now {r['accession']} -- {r['pick_reason']}")
        if unknown:
            print("      re-run redundancy.py to measure the new picks, then this again; "
                  "the loop ends when no pick moves")
    if args.check_local:
        present = sum(1 for r in rows if (C.ROOT / r["cds_path"]).exists())
        print(f"  CDS present on this machine: {present}/{len(rows)}")
    print(f"  wrote {plan_path.relative_to(C.ROOT)} and {list_path.relative_to(C.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
