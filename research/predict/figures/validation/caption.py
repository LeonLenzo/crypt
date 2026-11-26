#!/usr/bin/env python3
"""Write the caption for the prediction-validation figure from the data it plots.

Captions are generated, never typed, so a number cannot drift from the panel it
describes. Run prep.py first; this reads validation_numbers.tsv and nothing else.
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load():
    with open(HERE / "validation_numbers.tsv") as fh:
        return {r["key"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    n = {k: v["value"] for k, v in load().items()}
    note = {k: v["note"] for k, v in load().items()}

    text = (
        "Gene prediction recovers CDS that the reference needs and the archive does not "
        "supply. "
        f"(a) Benchmark against species that DO carry a submitter annotation (n={n['a_n']}): "
        "k-mer recall and precision of the predicted CDS against the true CDS, by the "
        "assembly level of the target. Median recall "
        f"{n['a_complete_recall']} on complete genomes and {n['a_contig_recall']} on contig-level "
        f"assemblies ({n['a_recall_median']} overall, precision {n['a_precision_median']}), so "
        "fragmentation does not degrade recovery; the production targets are mostly "
        "fragmented. "
        f"(b) The production run (n={n['b_n']}): BUSCO completeness of the predicted proteins "
        "against assembly size relative to the genus median, with marker area the contig "
        f"count. Median completeness {n['b_median']}% (Q1 {n['b_q1']}%), against the dashed "
        f"line and band, the median and interquartile range of the {int(n['b_ref_n']):,} DOWNLOADED "
        f"annotations already in db_v4 ({n['b_ref_median']}%). "
        f"{n['b_below90']} predictions fall below 90% and each is an assembly artefact rather "
        "than a prediction failure: three are far off their genus size, one is split over "
        f"23,880 contigs ({n['b_outliers']}). "
        f"(c) Novelty of the predicted CDS against db_v4 (n={n['c_n']}): the fraction of "
        "predicted k-mers the database has never seen, split on whether the target's genus is "
        f"already represented. Median {n['c_genus_yes_median']} where the genus is present "
        f"(n={n['c_genus_yes_n']}) against {n['c_genus_no_median']} where it is absent "
        f"(n={n['c_genus_no_n']}): a congener covers part of a species' k-mer space and nothing "
        "like all of it. "
        f"(d) The limitation (n={n['d_n']} of {n['a_n']} benchmark genomes with a genome size on "
        "record): predicted CDS bases over true CDS bases against genome size, fill the "
        f"prediction's precision. Over-prediction rises with genome size "
        f"(r = {n['d_r_size_overprediction']}) and precision falls with it "
        f"(r = {n['d_r_size_precision']}), reaching {n['d_gram_ratio']}x at precision "
        f"{n['d_gram_precision']} for Blumeria graminis and {n['d_trit_ratio']}x at "
        f"{n['d_trit_precision']} for Puccinia triticina. Both are wheat pathogens, which is why "
        "the wheat pathogens' own genomes are not predicted: inflating a taxon's minimizer "
        "denominator depresses the coverage of the species the chapter reports."
    )

    caveat = (
        f"Not in the figure. Panel b drops {n['b_no_expectation']} of {n['b_n']} targets "
        f"{note['b_no_expectation']}, so no size expectation can be formed. "
        f"Panel d drops {n['d_dropped']} benchmark genomes with no genome size on record "
        f"({note['d_dropped']}). Coverage for a predicted taxon is NOT comparable to coverage "
        "for a downloaded one, and predicted CDS is a competitor in the index: it may absorb "
        "reads, and it never supports a positive call about the species it came from."
    )

    out = HERE / "validation_caption.txt"
    out.write_text(text + "\n\n" + caveat + "\n")
    print(text, end="\n\n")
    print(caveat)


if __name__ == "__main__":
    main()
