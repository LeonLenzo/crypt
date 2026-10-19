#!/usr/bin/env python3
"""prep_saturation_cloud.py — per-SAMPLE accumulation cloud, with one universal null curve.

The 03_kraken charts plotted reads against distinct k-mers for ONE SPECIES at a time, faceted
by host, and the signal was visual: real detections sit on an accumulation line while
crash-outs fall off it to a floor. That works per species because every point shares one taxon
ceiling, so a single accumulation curve applies to the whole panel.

Per SAMPLE it does not, because every point is a different taxon with a different ceiling M.
Normalising both axes by M fixes it and makes the null universal:

    x = reads * kmers_per_read / M      sequencing redundancy relative to THAT taxon
    y = distinct / M                    which is exactly coverage
    null: y = 1 - exp(-x)               the same curve for every taxon, any genome size

A taxon whose reads are spread over its k-mer space like a real organism's sits ON that curve.
One whose reads pile onto a few loci falls FAR below it. `saturation` is the vertical distance,
y / (1 - exp(-x)), and is what _reports.saturation() returns.

Measured on PRJEB33109 (4 runs, 63 taxa): calls clearing 1% coverage sit at 0.070-0.5595, and
leakage with >=1,000 reads sits at 0.0001-0.0457. No overlap. 99,509 reads assigned to *Oryza
sativa* in a wheat sample score 0.0019 -- if they were real rice they would have produced
~6.6 M distinct minimizers rather than 24,388.

REAL CALLS SIT BELOW THE NULL, not on it, typically at 0.1-0.6. That is not a failure of the
model: the null assumes reads distribute uniformly over the k-mer space, and RNA-seq expression
is wildly non-uniform, so reads concentrate on highly expressed genes. The curve is a CEILING
on what uniform sampling could achieve, and the useful fact is that artefacts sit two further
orders of magnitude beneath it.

    python 02_classify/report/figures/prep_saturation_cloud.py --run ERR3383632
    python 02_classify/report/figures/prep_saturation_cloud.py --project PRJEB33109
"""

import argparse
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import DB_INSPECT, REPORTS, SILVER  # noqa: E402
from _reports import load_inspect, parse_report, saturation  # noqa: E402
from _util import atomic_write  # noqa: E402

TSV = HERE / "saturation_cloud.tsv"
CAPTION = HERE / "saturation_cloud_caption.txt"
K = 35
COV_FLOOR = 0.01
MIN_DISTINCT = 10_000

COLS = ["run", "BioProject", "taxid", "species", "rank", "reads", "distinct", "db_clade",
        "redundancy", "coverage", "saturation", "clears"]


def kmers_per_spot(read_len: float, layout: str) -> float:
    """k-mers contributed by ONE SPOT: (L-K+1) summed over its mates.

    ENA's read_len is bases per spot, so a PAIRED library's per-read length is half and the
    spot contributes twice that. Getting this wrong moves the null by a factor of 7 between a
    51 nt and a 150 nt library, which would make short-read samples look artefactual wholesale.
    """
    L = read_len / 2 if layout == "PAIRED" else read_len
    if L <= K:
        return 0.0
    per = L - K + 1
    return per * 2 if layout == "PAIRED" else per


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", default=[], help="run accession; repeatable")
    ap.add_argument("--project", help="every run with a report in this BioProject")
    ap.add_argument("--all", action="store_true", help="every run with a report, all projects")
    ap.add_argument("--reports", default=str(REPORTS), help="reports tree")
    ap.add_argument("--inspect", default=str(DB_INSPECT))
    ap.add_argument("--min-reads", type=int, default=1,
                    help="drop taxa below this many reads (default 1: keep everything, so the "
                         "cloud shows its own floor rather than a trimmed version of it)")
    args = ap.parse_args()

    insp = load_inspect(Path(args.inspect))
    depth = {r["Run"]: r for r in csv.DictReader(open(SILVER / "run_depth.tsv"),
                                                 delimiter="\t")}
    reports_root = Path(args.reports)

    wanted = []
    if args.all:
        for pdir in sorted(d for d in reports_root.iterdir() if d.is_dir()):
            for f in sorted(pdir.glob("*.report.tsv")):
                wanted.append((f.name.split(".")[0], pdir.name, f))
    if args.project:
        for p in sorted((reports_root / args.project).glob("*.report.tsv")):
            wanted.append((p.name.split(".")[0], args.project, p))
    for run in args.run:
        hits = list(reports_root.glob(f"*/{run}.report.tsv"))
        if not hits:
            sys.exit(f"no report for {run} under {reports_root}")
        wanted.append((run, hits[0].parent.name, hits[0]))
    if not wanted:
        sys.exit("give --run, --project or --all")

    rows = []
    for run, proj, path in wanted:
        d = depth.get(run)
        if not d:
            sys.exit(f"{run} has no depth record; run 03_cohort/fetch_depth.py")
        kps = kmers_per_spot(float(d["read_len"] or 0), d["layout"])
        if not kps:
            print(f"  {run}: read length {d['read_len']} <= k={K}, skipped")
            continue
        for rec in parse_report(path):
            if rec["rank"] != "S" or not rec["reads_clade"]:
                continue
            if rec["reads_clade"] < args.min_reads or rec["distinct"] is None:
                continue
            e = insp.get(rec["taxid"])
            if not e or not e["clade"]:
                continue
            m = e["clade"]
            sat = saturation(rec["distinct"], rec["reads_clade"], rec["taxid"], insp, kps)
            cov = rec["distinct"] / m
            rows.append({
                "run": run, "BioProject": proj, "taxid": rec["taxid"],
                "species": rec["name"], "rank": rec["rank"],
                "reads": rec["reads_clade"], "distinct": rec["distinct"], "db_clade": m,
                # x and y of the plot, both dimensionless and comparable across taxa.
                "redundancy": f"{rec['reads_clade'] * kps / m:.6g}",
                "coverage": f"{cov:.6g}",
                "saturation": f"{sat:.6g}" if sat is not None else "",
                "clears": "yes" if (cov >= COV_FLOOR
                                    and rec["distinct"] >= MIN_DISTINCT) else "no",
            })

    with atomic_write(TSV) as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    runs = sorted({r["run"] for r in rows})
    clears = [r for r in rows if r["clears"] == "yes"]
    print(f"-> {TSV.name}: {len(rows)} taxa over {len(runs)} runs")
    print(f"   {len(clears)} clear both floors")

    L, A = [], None
    A = L.append
    A("Per-sample accumulation cloud")
    A("Caption for report/figures/saturation_cloud.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_saturation_cloud.py -- do not hand-edit.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A("Every species-rank taxon Kraken2 assigned a read to in each sample, placed by how much")
    A("of that taxon's k-mer space the reads could have covered against how much they did.")
    A("Both axes are normalised by the taxon's own total in db_v4, so taxa with genomes")
    A("differing by four orders of magnitude are directly comparable on one panel:")
    A("")
    A("  x  redundancy = reads x k-mers per spot / the taxon's db_v4 minimizers")
    A("  y  coverage   = distinct minimizers observed / the same total")
    A("")
    A("THE CURVE IS THE NULL, y = 1 - exp(-x): where a taxon would sit if the reads assigned")
    A("to it were spread uniformly across its k-mer space. It is the same curve for every")
    A("taxon in every sample, which is what normalising buys. Saturation, the statistic in")
    A("_reports.py, is the vertical distance from it.")
    A("")
    A("REAL CALLS SIT BELOW THE CURVE, typically at 0.1-0.6 of it, and that is expected rather")
    A("than a failure: the null assumes uniform sampling and RNA-seq expression is not uniform,")
    A("so reads concentrate on highly expressed genes. The curve is a CEILING on what uniform")
    A("sampling could reach. What matters is that leakage sits two further orders of magnitude")
    A("beneath it.")
    A("")
    A("Measured on PRJEB33109, 4 runs, 63 taxa:")
    A("  clears 1% coverage         n= 9   median saturation 0.2795   range 0.0700 - 0.5595")
    A("  leakage (>=1k reads, <1%)  n=54   median saturation 0.0017   range 0.0001 - 0.0457")
    A("No overlap, and the separation holds per run (241x, 359x, 230x, 162x) BECAUSE the null")
    A("scales with read count. A raw reads-per-distinct ratio does NOT work and was tried: its")
    A("separation decayed 20.8x to 1.0x across those same four runs as depth rose, since a")
    A("deeply sequenced real organism accumulates redundant reads too.")
    A("")
    A("WHAT THE FALLEN POINTS ARE. Not contamination and not a classifier error, but REFERENCE")
    A("INCOMPLETENESS. A read from a host transcript absent from db_v4's CDS for that host,")
    A("whose ortholog in a related species IS present, has only the relative's k-mers to match")
    A("and takes the best home available; were both present the LCA would resolve it to the")
    A("shared ancestor instead. The clearest case: 99,509 reads assigned to Oryza sativa in a")
    A("wheat sample, saturation 0.0019 -- real rice would have produced ~6.6 M distinct")
    A("minimizers, not 24,388.")
    A("")
    A("WHAT THIS DOES NOT SAY. Saturation is measured, not yet calibrated: the separation above")
    A("comes from four wheat runs of a single Puccinia-dominated project, and no threshold is")
    A("set on it here. The floors drawn are the published operating point (1% coverage, and an")
    A(f"absolute {MIN_DISTINCT:,} distinct minimizers alongside it, since 1% of db_v4's smallest")
    A("species is 29 k-mers and reachable by cross-mapping alone).")
    A("")
    A("SHORT CAPTION")
    A("")
    A("Species-rank taxa per sample, positioned by k-mer redundancy against observed coverage,")
    A("both normalised by the taxon's db_v4 total so all taxa share one scale. The curve is the")
    A("null for reads spread uniformly over a taxon's k-mer space. Points near it behave like a")
    A("real organism; points far below it are reads piled onto a few loci, the signature of a")
    A("gap in the reference for the true source taxon.")
    A("")
    A("SOURCES")
    A(f"  reports  {reports_root}")
    A(f"  ceilings {Path(args.inspect).name}")
    A(f"  depths   run_depth.tsv (read length and layout set k-mers per spot)")
    A("")
    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))
    print(f"-> {CAPTION.name}")


if __name__ == "__main__":
    main()
