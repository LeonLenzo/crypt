#!/usr/bin/env python3
"""Tidy tables for the prediction-validation figure (panels A-D).

Four panels, four claims, each from data already on disk:

    A  fragmentation does not degrade recovery   n=60 benchmark, by assembly level
    B  predictions are as complete as downloads  n=135 BUSCO vs assembly quality
    C  congeners do not cover a species          n=59 novelty against db_v4
    D  large genomes are over-predicted          n=53 over-prediction vs genome size

Panels A, C and D read the n=60/n=59 benchmark; panel B reads the production run.
Nothing here is a release output yet: when supp_data_s4_predictions.tsv exists, panel B
should read that instead of worklist.tsv + all135_busco_results.tsv.
"""
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REVIEW = HERE.parents[3]                      # phd/01-review
PREDICT = HERE.parents[1]                     # research/predict

BENCHMARK = REVIEW / "research/02_classify/data/gold/predict_benchmark.tsv"
BUSCO135 = PREDICT / "data/all135_busco_results.tsv"
WORKLIST = PREDICT / "data/worklist.tsv"
SCREEN = REVIEW / "release/02_reference/data/curation/screen_assemblies.tsv"
DBV4_BUSCO = REVIEW / "research/superseded/03_kraken/utilities/busco/data/busco_scores.tsv"

LEVELS = ["Complete Genome", "Chromosome", "Scaffold", "Contig"]
MIN_GENUS_N = 3          # assemblies needed before a genus median is a usable expectation


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(name, fieldnames, rows):
    path = HERE / name
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{path.name:28s} {len(rows):5d} rows")
    return rows


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def genus_of(name):
    """First token, with the wastebasket brackets stripped: '[Candida] sp.' -> 'Candida'."""
    return name.split()[0].strip("[]")


# ─── panel A ──────────────────────────────────────────────────────────────────
# Arm A of the benchmark: predict a genome that HAS an annotation, score the
# prediction's k-mers against it. Long format so one geom draws both metrics.
def panel_a(bench):
    arm = [r for r in bench if r["arm"] == "annotated"]
    rows = []
    for r in arm:
        for metric in ("recall", "precision"):
            rows.append({
                "accession": r["accession"],
                "taxon": r["organism_name"],
                "assembly_level": r["assembly_level"],
                "metric": metric,
                "value": f"{float(r[metric]):.4f}",
            })
    write("panel_a_benchmark.tsv",
          ["accession", "taxon", "assembly_level", "metric", "value"], rows)
    return arm


# ─── panel B ──────────────────────────────────────────────────────────────────
# The production run has no annotation to score against, so completeness is BUSCO
# on the predicted proteins. Plotted against assembly quality, because an ab initio
# finder can only recover what the assembly contains: size deviation from the genus
# median catches partial and haplotig-inflated assemblies, contig count catches
# fragmentation. The comparator is db_v4's DOWNLOADED annotations, same metric.
def panel_b():
    sizes = defaultdict(list)
    for r in read(SCREEN):
        mb = num(r["genome_size_mb"])
        if mb and r["genus"]:
            sizes[r["genus"]].append(mb)

    worklist = {r["accession"]: r for r in read(WORKLIST)}
    rows = []
    for r in read(BUSCO135):
        w = worklist[r["accession"]]
        mb = num(w["genome_size_mb"])
        pool = sizes.get(w["genus"], [])
        expected = st.median(pool) if len(pool) >= MIN_GENUS_N else None
        rows.append({
            "accession": r["accession"],
            "taxon": w["taxon"],
            "tier": w["tier"],
            "assembly_level": w["assembly_level"],
            "contig_count": w["contig_count"],
            "genome_size_mb": w["genome_size_mb"],
            "genus_median_mb": f"{expected:.1f}" if expected else "",
            "size_dev_log2": f"{math.log2(mb / expected):.3f}" if expected and mb else "",
            "busco_complete": r["complete_pct"],
            "busco_lineage": r["lineage"],
            "n_markers": r["n_markers"],
        })
    rows.sort(key=lambda d: float(d["busco_complete"]))
    write("panel_b_production.tsv", list(rows[0]), rows)

    # the reference distribution: every assembly db_v4 actually took, with a BUSCO score
    ref = [num(r["complete_pct"]) for r in read(DBV4_BUSCO)
           if r["selected"] == "True" and num(r["complete_pct"]) is not None]
    q1, med, q3 = st.quantiles(ref, n=4)
    write("panel_b_reference.tsv", ["source", "n", "q1", "median", "q3", "frac_below_90"],
          [{"source": "downloaded_dbv4", "n": len(ref), "q1": f"{q1:.1f}",
            "median": f"{med:.1f}", "q3": f"{q3:.1f}",
            "frac_below_90": f"{sum(x < 90 for x in ref) / len(ref):.4f}"}])
    return rows, ref


# ─── panel C ──────────────────────────────────────────────────────────────────
# Arm B: predict a genome db_v4 does NOT contain, and ask how much of the prediction
# is k-mers db_v4 has never seen. Split on whether the GENUS is in db_v4, because the
# standing objection to prediction is that a congener already covers the species.
def panel_c(bench):
    in_db = {genus_of(r["organism_name"]) for r in read(DBV4_BUSCO) if r["selected"] == "True"}
    rows = []
    for r in bench:
        if r["arm"] != "unannotated":
            continue
        genus = genus_of(r["organism_name"])
        rows.append({
            "accession": r["accession"],
            "taxon": r["organism_name"],
            "genus": genus,
            "focal": r["focal"],
            "genus_in_dbv4": "yes" if genus in in_db else "no",
            "novelty": r["novelty"],
        })
    write("panel_c_novelty.tsv", list(rows[0]), rows)
    return rows


# ─── panel D ──────────────────────────────────────────────────────────────────
# The limitation. Over-prediction is predicted CDS bases over true CDS bases, which
# is the quantity that inflates a taxon's minimizer denominator and so depresses its
# coverage. Genome size is taken from the assembly screen, falling back to db_v4's
# own record; 7 of 60 are in neither and are dropped, named in the caption.
def panel_d(arm):
    size = {}
    for r in read(SCREEN):
        if num(r["genome_size_mb"]):
            size[r["accession"]] = num(r["genome_size_mb"])
    for r in read(DBV4_BUSCO):
        if r["accession"] not in size and num(r.get("total_length_mb")):
            size[r["accession"]] = num(r["total_length_mb"])

    rows, dropped = [], []
    for r in arm:
        mb = size.get(r["accession"])
        if mb is None:
            dropped.append(r["organism_name"])
            continue
        rows.append({
            "accession": r["accession"],
            "taxon": r["organism_name"],
            "assembly_level": r["assembly_level"],
            "genome_size_mb": f"{mb:.1f}",
            "pred_bases": r["pred_bases"],
            "true_bases": r["true_bases"],
            "overprediction": f"{float(r['pred_bases']) / float(r['true_bases']):.4f}",
            "precision": r["precision"],
            "recall": r["recall"],
        })
    rows.sort(key=lambda d: -float(d["overprediction"]))
    write("panel_d_overprediction.tsv", list(rows[0]), rows)
    if dropped:
        print(f"  dropped (no genome size on record): {len(dropped)} -> {', '.join(dropped)}")
    return rows, dropped


def pearson(xs, ys):
    mx, my = st.mean(xs), st.mean(ys)
    num_ = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = math.sqrt(sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys))
    return num_ / den


# ─── the numbers the caption quotes ───────────────────────────────────────────
def numbers(arm, b_rows, ref, c_rows, d_rows, dropped):
    out = []

    def put(key, value, note=""):
        out.append({"key": key, "value": value, "note": note})

    put("a_n", len(arm))
    put("a_recall_median", f"{st.median(float(r['recall']) for r in arm):.3f}")
    put("a_precision_median", f"{st.median(float(r['precision']) for r in arm):.3f}")
    for lvl in LEVELS:
        g = [r for r in arm if r["assembly_level"] == lvl]
        key = lvl.split()[0].lower()
        put(f"a_{key}_n", len(g))
        put(f"a_{key}_recall", f"{st.median(float(r['recall']) for r in g):.3f}")

    busco = [float(r["busco_complete"]) for r in b_rows]
    put("b_n", len(b_rows))
    put("b_median", f"{st.median(busco):.1f}")
    put("b_q1", f"{st.quantiles(busco)[0]:.1f}")
    put("b_below90", sum(x < 90 for x in busco))
    put("b_ref_n", len(ref))
    put("b_ref_median", f"{st.median(ref):.1f}")
    put("b_ref_below90_pct", f"{100 * sum(x < 90 for x in ref) / len(ref):.0f}")
    put("b_no_expectation", sum(1 for r in b_rows if not r["size_dev_log2"]),
        f"whose genus holds fewer than {MIN_GENUS_N} assemblies")
    put("b_outliers", "; ".join(
        f"{r['taxon']} {float(r['busco_complete']):.0f}%" for r in b_rows
        if float(r["busco_complete"]) < 90))

    for flag in ("yes", "no"):
        g = [float(r["novelty"]) for r in c_rows if r["genus_in_dbv4"] == flag]
        put(f"c_genus_{flag}_n", len(g))
        put(f"c_genus_{flag}_median", f"{st.median(g):.3f}")
    put("c_n", len(c_rows))

    ratios = [float(r["overprediction"]) for r in d_rows]
    mbs = [float(r["genome_size_mb"]) for r in d_rows]
    put("d_n", len(d_rows))
    put("d_dropped", len(dropped), "; ".join(dropped))
    put("d_median_ratio", f"{st.median(ratios):.2f}")
    put("d_r_size_overprediction", f"{pearson(mbs, ratios):+.2f}")
    put("d_r_size_precision", f"{pearson(mbs, [float(r['precision']) for r in d_rows]):+.2f}")
    put("d_r_size_recall", f"{pearson(mbs, [float(r['recall']) for r in d_rows]):+.2f}")
    for r in d_rows[:2]:
        tag = r["taxon"].split()[1][:4]
        put(f"d_{tag}_ratio", f"{float(r['overprediction']):.2f}")
        put(f"d_{tag}_precision", f"{float(r['precision']):.2f}")

    write("validation_numbers.tsv", ["key", "value", "note"], out)


def main():
    bench = read(BENCHMARK)
    arm = panel_a(bench)
    b_rows, ref = panel_b()
    c_rows = panel_c(bench)
    d_rows, dropped = panel_d(arm)
    numbers(arm, b_rows, ref, c_rows, d_rows, dropped)


if __name__ == "__main__":
    main()
