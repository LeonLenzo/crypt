#!/usr/bin/env python3
"""Concatenate the planned CDS into one tagged reference, and write the transcript table.

    python research/03_align/01_index/build_reference.py --label pleosporaceae
    python research/03_align/01_index/build_reference.py --label genus --dry-run

Reads  data/index_plan_<label>.tsv        what plan_index.py decided
       the CDS files it names             NCBI, predicted, and the host
Writes $ALIGN_WORK/<label>/reference.fa   the tagged reference, one record per transcript
       $ALIGN_WORK/<label>/transcripts.tsv  tag -> taxid, species, accession, gene, length
       $ALIGN_WORK/<label>/reference_summary.tsv  per taxon: transcripts, bases

Both products are OUTSIDE the repo: a 1.2 GB FASTA and a transcript table of ~750,000 rows are
not code and not a small figure. `reference_summary.tsv` is small and is the one worth reading.

## THE TAG IS THE WHOLE POINT

    {species_slug}|{accession}|{transcript_id}

bowtie2 truncates a reference name at the first whitespace, so everything that has to survive
into the BAM must be in that one token and nothing else may be. Three fields, carried from the
retired kallisto module where the scheme was designed, and it is what makes the two statistics
possible from one alignment: the first field aggregates to a species, the third keeps the gene
identity a gene-set test needs. The accession is in the middle because provenance has to be
readable off a BAM record without a join.

**Uniqueness is checked, not assumed.** Predicted gene models are named `g1`, `g2`, ... per
genome, so `g1_t1` appears in all 21 predicted files; the accession field is what separates
them. A collision would silently merge two organisms' reads, so the run fails on one.

## THE THREE HEADER DIALECTS, each with its own gene identifier

    NCBI       >kraken:taxid|5599|lcl|KV441469.1_cds_OAG25636.1_1 [locus_tag=CC77DRAFT_923880] ...
    predicted  >g1|g1_t1|JAVSBN010000001.1:12768-13755(+)
    Ensembl    >TraesCS7D03G0420900.1 cds chromosome:... gene:TraesCS7D03G0420900 ...

The gene is `locus_tag` for NCBI (falling back to `gene`, then `protein_id`), the first pipe
field for a predicted model, and the `gene:` token for Ensembl. **The gene is not the
transcript**, and for the host that distinction is the whole of question 3's quantification
decision: counts aggregate to the gene, and homoeologues collapse on top of that from the BAM.

The `kraken:taxid|<taxid>|` prefix on the NCBI files is 02_reference's doing, added so
kraken2-build could read the taxonomy off the sequence names. It is stripped here: the taxid
travels in the transcript table, where a join can use it, rather than in a name that has to
stay short.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common as C  # noqa: E402

HERE = Path(__file__).resolve().parent
PLANS = HERE / "data"
KRAKEN_PREFIX = re.compile(r"^kraken:taxid\|\d+\|")
BRACKET = re.compile(r"\[(\w+)=([^\]]*)\]")


def slug(name):
    """A species name as one whitespace-free token, reversibly enough to read."""
    return re.sub(r"_+", "_", re.sub(r"[^0-9a-z]+", "_", name.lower())).strip("_")


def open_maybe_gzip(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" \
        else open(path, "r", encoding="utf-8")


def parse_header(line, source):
    """Return (transcript_id, gene_id) for one FASTA header, by dialect."""
    head = line[1:].rstrip("\n")
    if source == "ncbi":
        head = KRAKEN_PREFIX.sub("", head)
        name = head.split(None, 1)[0]
        fields = dict(BRACKET.findall(head))
        gene = fields.get("locus_tag") or fields.get("gene") or fields.get("protein_id") or name
        return name, gene
    if source == "predicted":
        parts = head.split(None, 1)[0].split("|")
        gene = parts[0]
        transcript = parts[1] if len(parts) > 1 else parts[0]
        return transcript, gene
    if source == "ensembl_plants":
        name = head.split(None, 1)[0]
        gene = next((tok.split(":", 1)[1] for tok in head.split() if tok.startswith("gene:")),
                    name)
        return name, gene
    raise ValueError(f"unknown cds_source {source!r}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--label", default="pleosporaceae", help="which index plan to build")
    ap.add_argument("--dry-run", action="store_true",
                    help="check every input and the tagging, write nothing")
    args = ap.parse_args()

    plan_path = PLANS / f"index_plan_{args.label}.tsv"
    if not plan_path.exists():
        sys.exit(f"no plan at {plan_path}; run plan_index.py first")
    with open(plan_path, newline="", encoding="utf-8") as fh:
        plan = list(csv.DictReader(fh, delimiter="\t"))

    slugs = defaultdict(list)
    for row in plan:
        slugs[slug(row["taxon"])].append(row["taxon"])
    clashes = {s: v for s, v in slugs.items() if len(set(v)) > 1}
    if clashes:
        sys.exit(f"species slugs collide, which would merge two organisms: {clashes}")

    missing = [r["cds_path"] for r in plan if not (C.ROOT / r["cds_path"]).exists()]
    if missing:
        sys.exit(f"{len(missing)} CDS files are not on this machine, first: {missing[0]}\n"
                 f"pull them from Acacia (team bucket, single-file cp) before building")

    out_dir = C.ALIGN_WORK / args.label
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / "reference.fa"
    table = out_dir / "transcripts.tsv"
    summary = out_dir / "reference_summary.tsv"

    seen = set()
    per_taxon = []
    total_bases = total_transcripts = 0
    fa = None if args.dry_run else open(fasta, "w", encoding="utf-8")
    tb = None if args.dry_run else open(table, "w", newline="", encoding="utf-8")
    writer = None
    if tb:
        writer = csv.writer(tb, delimiter="\t", lineterminator="\n")
        writer.writerow(("tag", "taxid", "species", "accession", "cds_source", "gene",
                         "transcript", "role", "length"))
    try:
        for row in plan:
            species = slug(row["taxon"])
            accession = row["accession"] or Path(row["cds_path"]).parent.name
            n_transcripts = n_bases = 0
            tag = None
            length = 0
            with open_maybe_gzip(C.ROOT / row["cds_path"]) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if tag and writer:
                            writer.writerow((tag, row["taxid"], row["taxon"], accession,
                                             row["cds_source"], gene, transcript, row["role"],
                                             length))
                        transcript, gene = parse_header(line, row["cds_source"])
                        tag = f"{species}|{accession}|{transcript}"
                        if tag in seen:
                            sys.exit(f"duplicate tag {tag}: two records would merge")
                        seen.add(tag)
                        n_transcripts += 1
                        length = 0
                        if fa:
                            fa.write(f">{tag}\n")
                    else:
                        length += len(line.strip())
                        n_bases += len(line.strip())
                        if fa:
                            fa.write(line)
                if tag and writer:
                    writer.writerow((tag, row["taxid"], row["taxon"], accession,
                                     row["cds_source"], gene, transcript, row["role"], length))
            if not n_transcripts:
                sys.exit(f"{row['taxon']} ({row['cds_path']}) holds no records")
            per_taxon.append((row["taxon"], row["role"], row["cds_source"], accession,
                              n_transcripts, n_bases))
            total_transcripts += n_transcripts
            total_bases += n_bases
    finally:
        for handle in (fa, tb):
            if handle:
                handle.close()

    if not args.dry_run:
        with open(summary, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(("taxon", "role", "cds_source", "accession", "transcripts", "bases"))
            w.writerows(per_taxon)

    host = [t for t in per_taxon if t[1] == "host"]
    fungal = [t for t in per_taxon if t[1] != "host"]
    print(f"plan {args.label}: {len(plan)} entries" + (" (DRY RUN)" if args.dry_run else ""))
    print(f"  {total_transcripts:,} transcripts, {total_bases / 1e9:.3f} Gbp")
    print(f"  fungal: {sum(t[4] for t in fungal):,} transcripts, "
          f"{sum(t[5] for t in fungal) / 1e9:.3f} Gbp over {len(fungal)} taxa")
    if host:
        print(f"  host:   {host[0][4]:,} transcripts, {host[0][5] / 1e9:.3f} Gbp "
              f"({host[0][3]})")
    widest = sorted(fungal, key=lambda t: -t[5])[:3]
    print("  largest fungal contributors: " +
          ", ".join(f"{t[0]} {t[5] / 1e6:.0f} Mbp" for t in widest))
    if not args.dry_run:
        print(f"  wrote {fasta}")
        print(f"        {table}")
        print(f"        {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
