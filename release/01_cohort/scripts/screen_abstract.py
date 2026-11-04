#!/usr/bin/env python3
"""Classify each paper's study type from its title and abstract, before full text is chased.

    python release/01_cohort/scripts/screen_abstract.py

Reads  data/raw/abstracts.tsv          title and abstract for every DOI in the corpus
Writes data/curation/screen_abstract.tsv   one row per paper: data_type, signal, recommendation

WHY. The corpus is a semantic-search result, and semantic search matches TOPIC, not study
TYPE: it returns reviews, methods chapters and marker-assay papers that are about wheat
pathogenomics without generating the shotgun RNA-seq this cohort is built from. The abstract
is enough to tell those apart, and doing it here means a reviewer does not hand-download a
book chapter's PDF.

IT CLASSIFIES AND RECOMMENDS; IT DOES NOT DROP. Excluding on the absence of a keyword is the
same mistake as reading a host off the archive `organism` field: `pathogenomic` and
`high-throughput sequencing` are RNA-seq by another name, and a rule keyed on the literal
`RNA-seq` flagged a core cohort paper. So `data_type` is a best guess, `signal` is the phrase
that drove it, and `recommend` is keep / review / drop with `drop` reserved for the cases a
human would call the same way every time (a review, a pure marker assay). Everything uncertain
is `review`, for a person to settle.

data_type:
    rna_seq      shotgun RNA-seq / transcriptome / HTS of plant or infected tissue
    genome_only  a single virus genome cloned or assembled, with no field-survey sequencing
    marker       microsatellite, KASP, SNP array or GBS genotyping, no shotgun RNA-seq
    review       review, perspective, methods or book chapter: no primary data
    unclear      abstract too thin to call, or no abstract retrieved
"""

import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
DATA = MODULE / "data"
SRC = DATA / "raw" / "abstracts.tsv"
OUT = DATA / "curation" / "screen_abstract.tsv"

# Deliberately broad: every way a paper says "we shotgun-sequenced RNA".
RNASEQ = re.compile(
    r"rna[- ]?seq|transcriptom|metatranscriptom|\bmrna\b|total rna|poly[- ]?a|"
    r"ribo[- ]?deplet|small[- ]?rna|high[- ]?throughput sequenc|deep sequenc|"
    r"next[- ]generation sequenc|\bHTS\b|shotgun sequenc|illumina|nanopore|hiseq|"
    r"novaseq|miseq|pathogenomic|field pathogenomic|rna[- ]?tag[- ]?seq|tag[- ]?seq|"
    r"3'[- ]?(?:rna[- ]?seq|mrna)|quant[- ]?seq|\bcDNA librar|sequenc\w* (?:librar|read|"
    r"the (?:transcriptom|rna|genom))", re.I)
REVIEW = re.compile(
    r"\bthis review\b|\bwe review\b|\ba review of\b|\breview[,.]|perspective on|"
    r"book chapter|this chapter|this protocol describes|an advanced tool|an overview|"
    r"methods in molecular biology|springer protocols|state of the art", re.I)
MARKER = re.compile(
    r"microsatellite|\bSSR\b|\bKASP\b|kompetitive allele|\bSNP array\b|\bDArT\b|"
    r"cleaved amplified|\bCAPS\b marker|\bSCAR\b marker|genotyping[- ]by[- ]sequencing|"
    r"\bGBS\b", re.I)
GENOME_ONLY = re.compile(
    r"complete (?:genome|nucleotide) sequence|full[- ]length (?:genome|sequence)|"
    r"complete genome of|genome organization of|molecular characteriz\w+ of a novel", re.I)
# First Report notes: one-paragraph diagnostic records, rarely a survey
FIRST_REPORT = re.compile(r"first report of|first record of|detected in .* for the first time", re.I)


def classify(title, abstract):
    t = f"{title} {abstract}"
    seq = bool(RNASEQ.search(t))
    if not abstract:
        return "unclear", "no abstract retrieved", "review"
    if REVIEW.search(t):
        return "review", REVIEW.search(t).group(0), "drop"
    if MARKER.search(t) and not seq:
        return "marker", MARKER.search(t).group(0), "drop"
    if MARKER.search(t) and seq:
        # Derives markers from sequence, or is a marker paper citing transcriptomes: a person
        # must read it. The India pathotype papers sit here.
        return "rna_seq", f"{MARKER.search(t).group(0)} + {RNASEQ.search(t).group(0)}", "review"
    if GENOME_ONLY.search(t) and not re.search(
            r"survey|field|lineage|population|diversity|co-?infect|multiple", t, re.I):
        # A single assembled/cloned virus genome, with no survey language. These deposit a
        # consensus genome to GenBank, not raw reads to SRA, so there is nothing for a
        # read-based cohort to classify. Leon, drop. The survey filter above is what keeps a
        # genuine field virome study (which says "survey"/"lineage"/"diversity") out of here.
        return "genome_only", GENOME_ONLY.search(t).group(0), "drop"
    if seq:
        return "rna_seq", RNASEQ.search(t).group(0), "keep"
    if FIRST_REPORT.search(t):
        return "unclear", FIRST_REPORT.search(t).group(0), "review"
    return "unclear", "no sequencing language found", "review"


def main():
    rows = []
    with open(SRC, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            dt, sig, rec = classify(r.get("title", ""), r.get("abstract", ""))
            rows.append({"doi": r["doi"], "data_type": dt, "recommend": rec,
                         "signal": sig, "has_abstract": int(bool(r.get("abstract"))),
                         "open_access": r.get("open_access", ""), "title": r.get("title", "")})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["doi", "data_type", "recommend", "signal",
                                           "has_abstract", "open_access", "title"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    import collections
    print(f"-> {OUT.name}  {len(rows)} papers")
    print(f"   data_type: {dict(collections.Counter(r['data_type'] for r in rows).most_common())}")
    print(f"   recommend: {dict(collections.Counter(r['recommend'] for r in rows).most_common())}")
    for rec in ("drop", "review"):
        papers = [r for r in rows if r["recommend"] == rec]
        print(f"\n   {rec} ({len(papers)}):")
        for r in papers:
            print(f"     {r['doi']:<32}{r['data_type']:<12}{r['signal'][:26]:<28}{r['title'][:40]}")


if __name__ == "__main__":
    main()
