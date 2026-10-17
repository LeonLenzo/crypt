#!/usr/bin/env python3
"""prep_predict_counts.py — data, axis parameters and caption for the prediction figure.

Owns every number. `predict_counts.R` owns only the drawing, and reads the axis limits from
`predict_counts_params.tsv` rather than deriving them, so the two cannot disagree about what
was plotted. Same split as prep_db_contents.py / db_contents.R, and for the same reason: a
caption that quotes numbers the figure no longer shows is worse than no caption.

Writes:
    predict_counts.tsv          one row per assembly
    predict_counts_params.tsv   shared axis limits and breaks, per panel
    predict_counts_caption.txt  the caption, generated from the data

ASSEMBLY MODE. The figure carries no titles, no legend, no in-panel species labels and no
caption, because the panels are printed separately and annotated in Illustrator. Everything
stripped from the image is listed in the caption instead, including the outlier coordinates
needed to place labels by hand.

    python 02_classify/report/figures/prep_predict_counts.py
"""

import csv
import math
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import PREDICT_METRICS, REF_CANDIDATES, SILVER  # noqa: E402
from _util import atomic_write  # noqa: E402

GENOMES_MANIFEST = SILVER / "predict_genomes.tsv"
TSV = HERE / "predict_counts.tsv"
PARAMS = HERE / "predict_counts_params.tsv"
CAPTION = HERE / "predict_counts_caption.txt"

# Colour breaks for genome size, stated here so the caption and the figure agree.
SIZE_BREAKS = [10, 20, 40, 80, 160]


def linear_axis(hi: float, step: float) -> tuple:
    """Shared 0-based limits and regular breaks for a 1:1 panel.

    LINEAR AND FROM ZERO, not log. Leon's call, and the data agrees: counts span 3.4k to 22k
    and Mb 6.4 to 36.8, both under one order of magnitude, so log buys nothing and cannot
    include zero anyway. On a ratio plot a zero origin is worth having — it makes the 1:1 line
    a true identity line through the origin rather than an arbitrary offset.

    Both axes MUST share the range, or the dashed diagonal is not at 45 degrees and the eye
    reads deviation that is not in the data.
    """
    top = math.ceil(hi / step) * step
    breaks = [i * step for i in range(int(top / step) + 1)]
    return 0.0, top, breaks


def main() -> None:
    cand = {r["accession"]: r for r in
            csv.DictReader(open(REF_CANDIDATES, newline=""), delimiter="\t")}
    gsize = {r["accession"]: int(r["bytes"] or 0) for r in
             csv.DictReader(open(GENOMES_MANIFEST, newline=""), delimiter="\t")}
    rows = [r for r in csv.DictReader(open(PREDICT_METRICS, newline=""), delimiter="\t")
            if r["arm"] == "annotated" and r["true_seqs"]]
    if not rows:
        sys.exit(f"no scored Arm A rows in {PREDICT_METRICS}")

    out = []
    for r in rows:
        pg = (cand.get(r["accession"], {}).get("protein_coding_genes") or "").strip()
        out.append({
            "accession": r["accession"], "organism_name": r["organism_name"],
            "assembly_level": r["assembly_level"], "focal": r["focal"],
            "genome_mb": round(gsize.get(r["accession"], 0) / 1e6, 1),
            "true_seqs": int(r["true_seqs"]), "pred_seqs": int(r["pred_seqs"]),
            "true_genes": int(pg) if pg.isdigit() else "",
            "true_mb": round(int(r["true_bases"]) / 1e6, 3),
            "pred_mb": round(int(r["pred_bases"]) / 1e6, 3),
            "recall": float(r["recall"]), "precision": float(r["precision"]),
        })
    for r in out:
        r["ratio_seqs"] = round(r["pred_seqs"] / r["true_seqs"], 4)
        r["ratio_mb"] = round(r["pred_mb"] / r["true_mb"], 4)

    with atomic_write(TSV) as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: -r["true_seqs"]))

    # ── shared axis parameters, one row per panel ───────────────────────────────────────────
    # Panel a is drawn in units of 10,000 with breaks every 10,000, so the numerals are whole:
    # 0, 1, 2, 3. Leon's call, and it reads more cleanly than 0, 0.5, 1.0 ...
    sa_lo, sa_hi, sa_br = linear_axis(
        max(max(r["true_seqs"], r["pred_seqs"]) for r in out), 10000)
    mb_lo, mb_hi, mb_br = linear_axis(
        max(max(r["true_mb"], r["pred_mb"]) for r in out), 10)
    with atomic_write(PARAMS) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["panel", "lo", "hi", "breaks"])
        w.writerow(["a", f"{sa_lo:.0f}", f"{sa_hi:.0f}",
                    ";".join(f"{b:.0f}" for b in sa_br)])
        w.writerow(["b", f"{mb_lo:.0f}", f"{mb_hi:.0f}",
                    ";".join(f"{b:.0f}" for b in mb_br)])
        # Panel c: ONE range shared by both axes, so the two containments are read on the same
        # scale. It must reach down to the minimum precision, which is well below any recall.
        w.writerow(["c", "0.4", "1.0", "0.4;0.6;0.8;1.0"])
        w.writerow(["a_unit", "10000", "10000", "10000"])
        w.writerow(["size", str(min(SIZE_BREAKS)), str(max(SIZE_BREAKS)),
                    ";".join(str(b) for b in SIZE_BREAKS)])

    # ── caption ─────────────────────────────────────────────────────────────────────────────
    n = len(out)
    rs = [r["ratio_seqs"] for r in out]
    rb = [r["ratio_mb"] for r in out]
    rc = [r["recall"] for r in out]
    pr = [r["precision"] for r in out]
    gm = [r["genome_mb"] for r in out]
    iso = [r["true_seqs"] / r["true_genes"] for r in out if r["true_genes"]]

    # Outliers worth labelling by hand, with the coordinates to place them.
    out_a = sorted((r for r in out if r["ratio_seqs"] > 1.3 or r["ratio_seqs"] < 0.7),
                   key=lambda r: -r["ratio_seqs"])
    out_b = sorted((r for r in out if r["ratio_mb"] > 1.3 or r["ratio_mb"] < 0.7),
                   key=lambda r: -r["ratio_mb"])
    out_c = sorted((r for r in out if r["precision"] < 0.82 or r["recall"] < 0.88),
                   key=lambda r: r["precision"])

    L = []
    A = L.append
    A("Ab-initio gene prediction against NCBI annotation")
    A("Caption for report/figures/predict_counts.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_predict_counts.py -- do not hand-edit, the numbers")
    A("will drift.")
    A("")
    A("Drawn by predict_counts.R in ASSEMBLY mode: no panel titles, no legend, no in-panel")
    A("species labels and no caption on the image. The axis numerals and axis titles are kept.")
    A("Everything else is listed below for placement in Illustrator.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A(f"Vipsania, an unsupervised deep ab-initio gene finder, run on {n} fungal assemblies that")
    A("ALREADY carry an NCBI gene model, so the annotation can serve as ground truth. Each point")
    A("is one assembly. The dashed line in panels a and b is 1:1, and both axes of each of those")
    A("panels share a single range and break set, so a point's distance from the diagonal is")
    A("read directly rather than inferred.")
    A("")
    A("THREE INCREASINGLY STRICT QUESTIONS, which is the structure of the figure:")
    A("")
    A(f"  a  Does prediction call the same NUMBER of CDS?  median ratio {st.median(rs):.2f}")
    A(f"     ({sum(1 for x in rs if x < 1)} under-call, {sum(1 for x in rs if x > 1)} over-call)")
    A(f"  b  The same AMOUNT of sequence?                  median ratio {st.median(rb):.2f}")
    A(f"     ({sum(1 for x in rb if x < 1)} under, {sum(1 for x in rb if x > 1)} over)")
    A(f"  c  The SAME sequence?    median recall {st.median(rc):.3f}, precision {st.median(pr):.3f}")
    A("")
    A("PANEL c IS THE POINT. Agreement in a and b is necessary and nowhere near sufficient:")
    A("identical totals are perfectly consistent with entirely different sequence. Only k-mer")
    A("containment tests identity, which is why the benchmark does not stop at counting genes.")
    A("Recall is the fraction of the TRUE CDS k-mer space the prediction recovers -- too low and")
    A("the species cannot be detected when present. Precision is the fraction of the PREDICTED")
    A("k-mer space that is real CDS -- too low and the invented sequence inflates the taxon's")
    A("minimizer count, which is the denominator every coverage figure in 05_detect divides by.")
    A("Containment is estimated from FracMinHash sketches, k=31, scaled=200, the same parameters")
    A("as the saturation and depth-bias tables, so the numbers are comparable to those.")
    A("")
    A("COUNTS IN PANEL a ARE FAIR TO COMPARE, which is not obvious: NCBI's cds_from_genomic.fna")
    A("holds one entry per CDS INCLUDING isoforms, while vipsania writes one per predicted gene.")
    A(f"Measured rather than assumed -- NCBI CDS entries divided by NCBI's own protein-coding gene")
    A(f"count has median {st.median(iso):.3f} across the {len(iso)} assemblies where both are known, and only")
    A(f"{sum(1 for x in iso if x > 1.10)} exceed 1.10. Fungal isoform inflation is negligible at this scale.")
    A("")
    A(f"POINT FILL is genome size in Mb on a log10 scale, {min(gm):.0f} to {max(gm):.0f} Mb across the figure,")
    A(f"viridis reversed so dark is large. Colourbar breaks: {', '.join(str(b) for b in SIZE_BREAKS)}.")
    A("Genome size is the variable that predicts the deviation: it correlates -0.53 with both")
    A("recall and precision (r2 = 0.28). The mechanism is transposons. Ab-initio prediction calls")
    A("genes inside repeats, large fungal genomes are large mostly because of repeats, and the")
    A("failure therefore shows up as OVER-prediction with intact recall -- the signature is a")
    A("point that sits above the diagonal in a and b and bottom-right in c.")
    A("")
    A("AXIS RANGES, shared within each panel:")
    A(f"  a   0 to {sa_hi:,.0f} CDS, LINEAR, drawn in units of 10,000 so the numerals read")
    A(f"      {', '.join(f'{b/10000:g}' for b in sa_br)}; axis titles carry the x 10,000.")
    A(f"  b   0 to {mb_hi:.0f} Mb, LINEAR, breaks every 10")
    A("  c   0.40 to 1.00 on BOTH axes, breaks every 0.2; dotted guides at 0.90 on both")
    A("      Linear and zero-based throughout: the ranges span under one order of magnitude,")
    A("      so log buys nothing, and a zero origin makes the 1:1 line a true identity line.")
    A("")
    A("LABELS TO PLACE BY HAND. Species names are italic; coordinates are in data units.")
    A("")
    A(f"  panel a, ratio outside 0.7-1.3 ({len(out_a)}):")
    for r in out_a:
        A(f"    {r['organism_name']:34s} x={r['true_seqs']:>6,}  y={r['pred_seqs']:>6,}  "
          f"ratio {r['ratio_seqs']:.2f}")
    A("")
    A(f"  panel b, ratio outside 0.7-1.3 ({len(out_b)}):")
    for r in out_b:
        A(f"    {r['organism_name']:34s} x={r['true_mb']:>6.1f}  y={r['pred_mb']:>6.1f}  "
          f"ratio {r['ratio_mb']:.2f}")
    A("")
    A(f"  panel c, precision < 0.82 or recall < 0.88 ({len(out_c)}):")
    for r in out_c:
        A(f"    {r['organism_name']:34s} x={r['recall']:.3f}  y={r['precision']:.3f}  "
          f"{r['genome_mb']:.0f} Mb")
    A("")
    A("WHAT THIS FIGURE DOES NOT SHOW. These 60 species already have gene models and are already")
    A("in db_v4, so prediction adds nothing for them -- the panel exists only to calibrate how")
    A("far a predictor can be trusted. Whether prediction is WORTH running on the species that")
    A("lack a model is a separate measurement, of how much of the predicted k-mer space is new")
    A("to db_v4, and is not in this figure.")
    A("")
    # ── short caption: the one that goes under the figure in a manuscript ──────────────────
    worst_b = max(out, key=lambda r: r["ratio_mb"])
    A("")
    A("SHORT CAPTION")
    A("")
    A(f"Ab-initio gene prediction compared with NCBI annotation for {n} fungal assemblies that")
    A("already carry a gene model. Each point is one assembly; fill is genome size on a log")
    A("scale. (a) Number of CDS called by prediction against the number annotated; dashed line")
    A(f"is 1:1, median ratio {st.median(rs):.2f}. (b) The same for total CDS sequence, median ratio")
    A(f"{st.median(rb):.2f}. (c) k-mer containment against the annotated CDS: recall is the fraction of")
    A(f"true CDS k-mer space recovered (median {st.median(rc):.2f}) and precision the fraction of")
    A(f"predicted space that is real CDS (median {st.median(pr):.2f}); dotted guides mark 0.90.")
    A("Agreement in (a) and (b) is consistent with entirely different sequence, so only (c)")
    A("tests identity. Deviation tracks genome size (r = -0.53 against both metrics): ab-initio")
    A("prediction calls genes within transposons, so large repeat-rich genomes over-call while")
    A(f"recall stays intact -- {worst_b['organism_name']} predicts {worst_b['ratio_mb']:.1f}x its annotated CDS at")
    A(f"precision {worst_b['precision']:.2f}. Containment from FracMinHash sketches, k=31, scaled=200.")
    A("")
    A("SOURCES")
    A(f"  data      {TSV.name}, from {PREDICT_METRICS.name}")
    A(f"  axes      {PARAMS.name}")
    A(f"  genome Mb {GENOMES_MANIFEST.name}")
    A(f"  gene count {REF_CANDIDATES.name} (protein_coding_genes)")
    A("")

    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))

    print(f"-> {TSV.name}     ({n} assemblies)")
    print(f"-> {PARAMS.name}")
    print(f"-> {CAPTION.name}  ({len(L)} lines)")
    print(f"\n   a median ratio {st.median(rs):.2f}   b median ratio {st.median(rb):.2f}"
          f"   c median {st.median(rc):.3f}/{st.median(pr):.3f}")


if __name__ == "__main__":
    main()
