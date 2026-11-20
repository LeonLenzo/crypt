#!/usr/bin/env python3
"""The purpose test: do the 135 predicted references move reads down to species.

    python release/03_classify/scripts/compare_arms.py
    python release/03_classify/scripts/compare_arms.py --runs 200

Reads  data/raw/reports/db1p/<BioProject>/<Run>.report.tsv   the baseline arm
       data/raw/reports/db2/<BioProject>/<Run>.report.tsv     the same runs, plus predictions
       data/raw/candidates_db2.tsv                            what the competitors absorbed
Writes data/outputs/arm_resolution.tsv   one row per run: where its reads came to rest, per arm

db2 is db1p plus 135 gene-predicted genomes and NOTHING else, so a difference between the arms
is attributable to those references. The question is not whether more reads are classified; it
is whether reads that previously came to rest at GENUS or above now rest at SPECIES. A read in
a genus with no annotated member cannot be resolved below genus however abundant it is, and
that unresolved fraction is what the prediction arm was built to reduce.

**THE PROFILE IS MEASURED OVER FUNGAL READS, and the all-reads profile is reported beside it
as context.** The host is most of a plant RNA-seq library and resolves at species, so it sets
the all-reads denominator and buries any fungal movement below the second decimal place. The
fungal restriction is a lineage test on the Fungi node, not a name match.

**Only runs classified in BOTH arms are compared**, and the count is printed. A per-run paired
difference is the right statistic here because run depth and library composition vary by orders
of magnitude across the cohort, so a pooled total would be dominated by the deepest projects.

**The predicted taxa are competitors, never calls.** What they absorbed is reported here as
reads and distinct minimizers, which is the measurement the purpose test rests on, and it is
the one place those numbers are allowed to appear as a result.
"""

import argparse
import collections
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _reports import parse_report  # noqa: E402

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
REPORTS = MODULE / "data" / "raw" / "reports"
OUT = MODULE / "data" / "outputs" / "arm_resolution.tsv"

ARMS = ("db1p", "db2")
BUCKETS = ("unclassified", "species", "genus", "above_genus")
FUNGI_TAXID = 4751
# The fungal buckets are the ones the purpose test is about. The host swamps the all-reads
# denominator -- wheat alone is about 79% of classified reads and resolves at species -- so a
# shift of a few hundred thousand fungal reads is invisible there. Measured both ways because
# the all-reads profile is still what says whether the library as a whole moved.
FBUCKETS = ("f_reads", "f_species", "f_genus", "f_above_genus")
OUT_COLS = (["run", "BioProject"]
            + [f"{a}_{b}" for a in ARMS for b in ("reads",) + BUCKETS + FBUCKETS]
            + ["delta_species", "delta_genus", "delta_above_genus", "delta_unclassified",
               "delta_f_species", "delta_f_genus", "delta_f_above_genus"])


def bucket(rank: str) -> str:
    if rank == "U":
        return "unclassified"
    if rank.startswith("S"):
        return "species"
    if rank.startswith("G"):
        return "genus"
    return "above_genus"


def profile(report: Path) -> dict:
    """Reads resting AT each rank bucket, from reads_direct, plus the run's total.

    `reads_direct` and not `reads_clade`: a clade count includes every descendant, so summing
    it across ranks counts the same read many times. The direct column is where a read actually
    came to rest, which is the whole question.
    """
    out = dict.fromkeys(BUCKETS + FBUCKETS, 0)
    total = 0
    for rec in parse_report(report):
        b = bucket(rec["rank"])
        out[b] += rec["reads_direct"]
        total += rec["reads_direct"]
        if rec["taxid"] == FUNGI_TAXID or any(t == FUNGI_TAXID for _, t, _ in rec["lineage"]):
            out["f_reads"] += rec["reads_direct"]
            if b != "unclassified":
                out[f"f_{b}"] += rec["reads_direct"]
    out["reads"] = total
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=int, help="stop after this many paired runs (a smoke test)")
    args = ap.parse_args()

    found = {}
    for arm in ARMS:
        root = REPORTS / arm
        if not root.is_dir():
            sys.exit(f"no reports for arm {arm} at {root}")
        found[arm] = {r.name.split(".")[0]: r
                      for proj in sorted(p for p in root.iterdir() if p.is_dir())
                      for r in sorted(proj.glob("*.report.tsv"))}
    paired = sorted(set(found["db1p"]) & set(found["db2"]))
    only = {a: len(set(found[a]) - set(paired)) for a in ARMS}
    if not paired:
        sys.exit("no run is classified in both arms yet")
    if args.runs:
        paired = paired[:args.runs]

    rows = []
    for run in paired:
        prof = {a: profile(found[a][run]) for a in ARMS}
        row = {"run": run, "BioProject": found["db2"][run].parent.name}
        for a in ARMS:
            for key in ("reads",) + BUCKETS + FBUCKETS:
                row[f"{a}_{key}"] = prof[a][key]
        for b in BUCKETS + ("f_species", "f_genus", "f_above_genus"):
            row[f"delta_{b}"] = prof["db2"][b] - prof["db1p"][b]
        rows.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    def share(row, arm, b, denom):
        d = row[f"{arm}_{denom}"]
        if denom == "reads":
            d -= row[f"{arm}_unclassified"]
        return row[f"{arm}_{b}"] / d if d else None

    print(f"{len(rows):,} runs classified in both arms "
          f"(db1p only {only['db1p']:,}, db2 only {only['db2']:,})")

    for label, buckets, denom in (
            ("FUNGAL reads only -- this is the purpose test",
             ("f_species", "f_genus", "f_above_genus"), "f_reads"),
            ("all classified reads, which the host dominates",
             ("species", "genus", "above_genus"), "reads")):
        print(f"\n{label}")
        print(f"  {'bucket':14s} {'db1p':>9s} {'db2':>9s} {'shift':>9s}   median share per run")
        for b in buckets:
            a1 = [s for r in rows if (s := share(r, "db1p", b, denom)) is not None]
            a2 = [s for r in rows if (s := share(r, "db2", b, denom)) is not None]
            if not a1:
                continue
            m1, m2 = statistics.median(a1), statistics.median(a2)
            print(f"  {b:14s} {m1:8.2%} {m2:8.2%} {m2 - m1:+8.2%}")

    for b in ("f_species", "f_genus", "species", "genus"):
        d = [r[f"delta_{b}"] for r in rows]
        up = sum(1 for x in d if x > 0)
        print(f"\n{b}: more in db2 for {up:,} of {len(rows):,} runs, "
              f"median per-run change {statistics.median(d):+,.0f}, "
              f"total {sum(d):+,.0f}")

    cand = MODULE / "data" / "raw" / "candidates_db2.tsv"
    if cand.is_file():
        comp = [c for c in csv.DictReader(open(cand, newline="", encoding="utf-8"),
                                         delimiter="\t")
                if c["call"] == "competitor_only"]
        if comp:
            reads = sum(int(c["reads"]) for c in comp)
            taxa = collections.Counter(c["taxon"] for c in comp)
            print(f"\npredicted competitors absorbed {reads:,} reads over {len(comp):,} rows "
                  f"and {len(taxa)} species, none of them callable")
            print("  most absorbed: "
                  + ", ".join(f"{t} ({n} runs)" for t, n in taxa.most_common(5)))
    else:
        print(f"\nnote: {cand.name} not present, so what the competitors absorbed is not "
              "reported. Run detect.py --arm db2 first.")
    print(f"\n-> {OUT.relative_to(MODULE)}  {len(rows):,} paired runs")


if __name__ == "__main__":
    main()
