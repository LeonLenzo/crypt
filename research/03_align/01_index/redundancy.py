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

## THE STATISTIC

For each transcript, sample every 16th 31-mer, then ask whether another transcript of the same
assembly shares at least 90% of that sample. `duplicated` is the share of transcripts for which
one does.

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
STRIDE = 16
SHARE = 0.90
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
    """(n_transcripts, median length, duplicated fraction, worst partner count)."""
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
    return (n, int(statistics.median(lengths)) if lengths else 0,
            duplicated / n if n else 0.0, max(partners) if partners else 0)


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
        n, median_len, duplicated, worst = measure(cds_path(accession))
        rows.append({"accession": accession, "taxon": taxon, "flags": flags,
                     "transcripts": n, "median_length": median_len,
                     "duplicated": f"{duplicated:.3f}", "max_partners": worst})
        print(f"  {taxon:<34} {accession:<18} {n:>7,} transcripts  "
              f"duplicated {duplicated:.3f}  max partners {worst}"
              + (f"  [{flags}]" if flags else ""))

    if not args.cds:
        out = PLANS / f"redundancy_{args.label}.tsv"
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        values = sorted(float(r["duplicated"]) for r in rows)
        print(f"\n{len(rows)} assemblies measured, duplicated fraction "
              f"median {statistics.median(values):.3f}, "
              f"range {values[0]:.3f} to {values[-1]:.3f}")
        print(f"wrote {out.relative_to(C.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
