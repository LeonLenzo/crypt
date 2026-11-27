#!/usr/bin/env python3
"""Measure within-species transcript redundancy, so the sequence decides what metadata cannot.

    python research/03_align/01_index/redundancy.py --label pleosporaceae
    python research/03_align/01_index/redundancy.py --cds GCA_016097525.1 GCA_001642055.1

Reads  data/index_plan_<label>.tsv   the picks, or --cds for a head-to-head on accessions
       the CDS files themselves
Writes data/redundancy_<label>.tsv   one row per assembly measured (small, tracked)

`plan_index.py` FLAGS a pick whose gene count disagrees with its species and refuses to decide
why, because the metadata cannot: a doubled gene count is a hybrid, a duplicated annotation, or
a genuinely larger genome, and NCBI's table says the same thing in all three cases. **This is
the measurement that separates them.** A duplicated annotation carries a near-identical twin
for most of its loci; a larger genome carries more DISTINCT genes.

## TWO STATISTICS, BECAUSE THE FIRST ONE ANSWERED NO AND THE QUESTION WAS STILL OPEN

**`duplicated`**: for each transcript, sample every 16th 31-mer, then ask whether another
transcript of the same assembly shares at least 90% of that sample. The share of transcripts
for which one does.

**`median_length` and `short_share`**: the median CDS length, and the share of CDS under 300 bp.

The second pair exists because the first answered the *A. alternata* question with a flat no and
the assembly was still wrong, measured 2026-11-26. `GCA_016097525.1` declares 24,347 genes
against the species median of 11,928, and `duplicated` is **0.002**: it carries no duplicated
loci at all. What it carries is an OVER-ANNOTATION. Its median CDS is 576 bp with 21.2% under
300 bp, against 1,149 bp and 3.9% for `GCA_001642055.1`, and the two hold 21.1 and 18.5 Mbp of
coding sequence. Same genome, same coding content, nearly twice as many records, each about
half as long. Every other pick in the Pleosporaceae plan sits between 1,059 and 1,347 bp.

**That distinction decides a presence fact.** Gene breadth, the count of distinct genes carrying
species-unique reads, is one of the three things job 1 reads, and 5,153 sub-300 bp ORFs inflate
it while adding exactly the short targets cross-mapping concentrates on. So the length
distribution is not a tidiness check; it is the statistic that protects the measurement.

`duplicated` stays, as a null result worth having: over the 62 Pleosporaceae picks its median is
0.000 and its maximum 0.104 (*Stemphylium vesicarium*), which is gene families rather than
duplication. Locus duplication is not a problem in this index, and now that is known rather
than assumed.

Sampling rather than every k-mer because the answer does not need the precision and the cost is
16x: a 13,000-transcript assembly reduces to about a million hashes. **It is a redundancy
measure and not an identity measure.** Two transcripts sharing 90% of sampled k-mers are
near-duplicates at nucleotide level, which is what matters for an index, but the number is not
a percent identity and should not be quoted as one.

A real gene family raises this a little and a duplicated assembly raises it a lot, so the run
prints the whole plan's distribution rather than a threshold: the comparison is between
assemblies of the SAME species, where the only difference is the annotation.
"""

import argparse
import csv
import gzip
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common as C  # noqa: E402

HERE = Path(__file__).resolve().parent
PLANS = HERE / "data"
K = 31
MIN_MEDIAN_CDS = 800             # bp; the Pleosporaceae picks run 1,059-1,347 and the one
                                 # rejected annotation sits at 576
MAX_SHORT_SHARE = 0.15           # of CDS under SHORT_CDS; normal is ~0.04, the bad one 0.21
MAX_DUPLICATED = 0.30            # the observed maximum from real gene families is 0.104
STRIDE = 16
SHARE = 0.90
SHORT_CDS = 300                  # bp; below this a CDS is a short ORF, and an over-annotation
                                 # is mostly these. Not a quality judgement on one gene.
KRAKEN_PREFIX = re.compile(r"^kraken:taxid\|\d+\|")


def records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    name, chunks = None, []
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name = KRAKEN_PREFIX.sub("", line[1:].split(None, 1)[0])
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        yield name, "".join(chunks)


def sampled_kmers(sequence):
    up = sequence.upper()
    return {hash(up[i:i + K]) for i in range(0, max(len(up) - K + 1, 0), STRIDE)}


def measure(path):
    """(n, median length, duplicated fraction, worst partner count, short share, coding bp)."""
    kmers, lengths = {}, []
    for name, sequence in records(path):
        if len(sequence) < K:
            continue
        kmers[name] = sampled_kmers(sequence)
        lengths.append(len(sequence))
    index = defaultdict(list)
    for name, ks in kmers.items():
        for k in ks:
            index[k].append(name)
    duplicated = 0
    partners = []
    for name, ks in kmers.items():
        if not ks:
            continue
        shared = defaultdict(int)
        for k in ks:
            for other in index[k]:
                if other != name:
                    shared[other] += 1
        hits = [o for o, n in shared.items() if n >= SHARE * len(ks)]
        if hits:
            duplicated += 1
        partners.append(len(hits))
    n = len(kmers)
    short = sum(1 for length in lengths if length < SHORT_CDS) / len(lengths) if lengths else 0.0
    return (n, int(statistics.median(lengths)) if lengths else 0,
            duplicated / n if n else 0.0, max(partners) if partners else 0,
            short, sum(lengths))


def cds_path(accession):
    for candidate in (C.CDS_DIR / accession / "cds_from_genomic.fna",
                      C.ROOT / "research/predict/data/predicted" / accession / "coding.fa"):
        if candidate.exists():
            return candidate
    sys.exit(f"no CDS on this machine for {accession}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--label", default="pleosporaceae", help="measure every pick in this plan")
    ap.add_argument("--cds", nargs="+", metavar="ACCESSION",
                    help="measure these accessions instead, for a head-to-head")
    ap.add_argument("--flagged-only", action="store_true",
                    help="measure only the picks plan_index.py flagged, and their runners-up")
    args = ap.parse_args()

    targets = []
    if args.cds:
        targets = [(a, a, "", "") for a in args.cds]
    else:
        with open(PLANS / f"index_plan_{args.label}.tsv", newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["role"] == "host":
                    continue
                accession = row["accession"] or Path(row["cds_path"]).parent.name
                if args.flagged_only and not row["flags"]:
                    continue
                targets.append((accession, row["taxon"], row["flags"], row["runners_up"]))
                if args.flagged_only:
                    for other in filter(None, row["runners_up"].split(";")):
                        targets.append((other.split("(")[0], row["taxon"] + " (runner-up)",
                                        "", ""))

    rows = []
    for accession, taxon, flags, _ in targets:
        n, median_len, duplicated, worst, short, coding = measure(cds_path(accession))
        verdict = ";".join(v for v, on in (
            ("over_annotated", median_len and median_len < MIN_MEDIAN_CDS),
            ("short_orf_heavy", short > MAX_SHORT_SHARE),
            ("duplicated_loci", duplicated > MAX_DUPLICATED)) if on) or "ok"
        rows.append({"accession": accession, "taxon": taxon, "flags": flags,
                     "transcripts": n, "median_length": median_len,
                     "short_share": f"{short:.3f}", "coding_mbp": f"{coding / 1e6:.1f}",
                     "duplicated": f"{duplicated:.3f}", "max_partners": worst,
                     "verdict": verdict})
        print(f"  {taxon:<34} {accession:<18} {n:>7,} tx  median {median_len:>5} bp  "
              f"short {short:>5.1%}  dup {duplicated:.3f}  {verdict}"
              + (f"  [{flags}]" if flags else ""))

    if not args.cds:
        out = PLANS / f"redundancy_{args.label}.tsv"
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        values = sorted(float(r["duplicated"]) for r in rows)
        lengths = sorted(r["median_length"] for r in rows)
        failing = [r for r in rows if r["verdict"] != "ok"]
        print(f"\n{len(rows)} assemblies measured")
        print(f"  duplicated fraction: median {statistics.median(values):.3f}, "
              f"range {values[0]:.3f} to {values[-1]:.3f}")
        print(f"  median CDS length:   median {statistics.median(lengths):.0f} bp, "
              f"range {lengths[0]} to {lengths[-1]}")
        if failing:
            print(f"  {len(failing)} picks FAIL the sequence gate and should be substituted:")
            for r in failing:
                print(f"      {r['taxon']} {r['accession']}: {r['verdict']} "
                      f"({r['median_length']} bp, {r['short_share']} short, "
                      f"{r['duplicated']} duplicated)")
        else:
            print("  every pick passes the sequence gate")
        print(f"wrote {out.relative_to(C.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
