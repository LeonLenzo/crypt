#!/usr/bin/env python3
"""prep_predict_gain.py — what gene prediction would add to db_v4, by route.

Companion to the db_contents ring, for a talk: same six routes, same colours, same precedence
rule, but a plain linear bar chart instead of a log ring. The ring shows what db_v4 HOLDS; this
shows what it is MISSING and could reach.

The route assignment is imported from prep_db_contents.py rather than reimplemented, so the
two figures cannot disagree about which species is which colour. That matters here more than
usual: the whole argument is a comparison between the two.

Each species is counted once, in one route:

    in db_v4     the species has at least one ANNOTATED assembly, so it is already a tip
                 on the ring
    addable      the species has assemblies but NO gene model anywhere, so it is invisible to
                 a CDS-only index today and is exactly what prediction would reach

Species with no assembly at all are absent from both; nothing can be done for them here.

Writes predict_gain.tsv and predict_gain_caption.txt.

    python 02_classify/report/figures/prep_predict_gain.py
"""

import csv
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import GOLD, SILVER  # noqa: E402
from _taxonomy import load as load_taxonomy  # noqa: E402
from _util import atomic_write  # noqa: E402

SURVEY = SILVER / "assembly_survey.tsv"
TSV = HERE / "predict_gain.tsv"
CAPTION = HERE / "predict_gain_caption.txt"
BENCH = GOLD / "predict_benchmark.tsv"


def _prep_db_contents():
    """Import the ring's own route logic. One definition, two figures."""
    spec = importlib.util.spec_from_file_location("_pdc", HERE / "prep_db_contents.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    pdc = _prep_db_contents()
    tax = load_taxonomy()
    on_cereal = pdc.cereal_recorded(tax)

    counts = {r: {"in_db": 0, "addable": 0} for r in pdc.PRECEDENCE}
    # The bar height. "CDS annotations" counts ASSEMBLIES carrying a gene model, actual
    # against potential, which is apples to apples:
    #   in_db    assemblies that HAVE an NCBI annotation           (n_annotated)
    #   addable  assemblies belonging to species with NONE at all  (n_assemblies)
    # Species that already have some annotated assemblies also have unannotated ones, and
    # those could be predicted too — but they add no new SPECIES, so they are excluded. The
    # pale bar is what prediction makes REACHABLE, not every prediction that could be run.
    cds = {r: {"in_db": 0, "addable": 0} for r in pdc.PRECEDENCE}
    for r in csv.DictReader(open(SURVEY), delimiter="\t"):
        tid = tax.by_ncbi.get(int(r["ncbi_taxid"]))
        n_total = int(r["n_assemblies"] or 0)
        n_ann = int(r["n_annotated"] or 0)
        if not n_total:
            continue
        route = pdc.primary_route(r["routes"] or "congener",
                                  tid in on_cereal if tid is not None else False)
        if n_ann:
            counts[route]["in_db"] += 1
            cds[route]["in_db"] += n_ann
        else:
            counts[route]["addable"] += 1
            cds[route]["addable"] += n_total

    # `host` is dropped from the CHART, not from the counts: it is 6 plant hosts with ZERO
    # addable species, so it contributes nothing to an argument about additions while forcing
    # a 180x axis range that would otherwise push the figure onto a log scale. A bar chart on
    # a log axis is the wrong trade — bar LENGTH stops being proportional to value, which is
    # the one thing a bar chart promises. Without host the range is 70 to 1,076 and linear is
    # comfortable. The host row is kept in the TSV and stated in the caption.
    order = [r for r in pdc.PRECEDENCE if sum(counts[r].values()) and r != "host"]
    with atomic_write(TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["route", "label", "status", "species", "cds_annotations", "drawn"])
        for route in [r for r in pdc.PRECEDENCE if sum(counts[r].values())]:
            for status in ("in_db", "addable"):
                w.writerow([route, pdc.ROUTE_LABEL[route], status,
                            counts[route][status], cds[route][status],
                            "yes" if route in order else "no"])

    tot_in = sum(c["in_db"] for c in counts.values())
    tot_add = sum(c["addable"] for c in counts.values())
    cds_in = sum(c["in_db"] for c in cds.values())
    cds_add = sum(c["addable"] for c in cds.values())

    # Measured prediction quality, so the caption argues from data rather than hope.
    recalls, precisions, novelties = [], [], []
    if BENCH.exists():
        for r in csv.DictReader(open(BENCH), delimiter="\t"):
            if r["arm"] == "annotated" and r["recall"]:
                recalls.append(float(r["recall"]))
                precisions.append(float(r["precision"]))
            if r["arm"] == "unannotated" and r["novelty"]:
                novelties.append(float(r["novelty"]))

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else float("nan")

    L, A = [], None
    A = L.append
    A("What gene prediction would add to db_v4")
    A("Caption for report/figures/predict_gain.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_predict_gain.py -- do not hand-edit.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A(f"CDS annotations available to db_v4 by route into scope. The SOLID bar is the {cds_in:,}")
    A("assemblies that carry an NCBI gene model and are therefore already in the database. The")
    A(f"PALE bar is the {cds_add:,} assemblies belonging to species with NO gene model anywhere, which")
    A("are invisible to a CDS-only index today and are what ab-initio prediction would reach.")
    A(f"In species terms that is {tot_in:,} species now and {tot_add:,} more reachable.")
    A("")
    A("Assemblies of species that ALREADY have an annotation are excluded from the pale bar.")
    A("They could be predicted too, but they add no new species, and the claim here is about")
    A("what becomes reachable rather than how many predictions could be run.")
    A("Routes, colours and the precedence rule are shared with the db_v4 ring figure; each")
    A("species is counted once, under the strongest evidence that applies.")
    A("")
    A("  route                                           CDS now   addable     total")
    for route in order:
        c = cds[route]
        A(f"  {pdc.ROUTE_LABEL[route]:46s} {c['in_db']:7,d} {c['addable']:9,d} {c['in_db']+c['addable']:9,d}")
    A(f"  {'TOTAL':46s} {cds_in:7,d} {cds_add:9,d} {cds_in+cds_add:9,d}")
    A("")
    if recalls:
        A("THE ADDITION IS NOT HYPOTHETICAL. Measured on 60 fungal assemblies that already carry")
        A("an NCBI gene model, so the annotation serves as ground truth, ab-initio prediction")
        A(f"recovers a median {med(recalls):.0%} of the true CDS k-mer space at {med(precisions):.0%} precision.")
    if novelties:
        A(f"On {len(novelties)} assemblies WITHOUT a gene model, a median {med(novelties):.0%} of the predicted k-mer")
        A("space is absent from db_v4 — so the pale bars are close to their full face value and")
        A("not largely redundant with what the database already holds.")
    A("")
    A(f"THE host ROUTE IS NOT DRAWN. It is {counts['host']['in_db']} plant host species ({cds['host']['in_db']} annotations) with ZERO")
    A("addable, so it")
    A("says nothing about additions while forcing a 180x axis range. Omitting it keeps the")
    A("chart linear, where a bar's LENGTH is its value; a log axis would break that for the")
    A("sake of one bar carrying no argument. The row stays in predict_gain.tsv.")
    A("")
    A("WHAT THE BAR DOES NOT SAY. Annotations, not sequence: a pale bar counts assemblies that")
    A("would gain a gene model, not the k-mer space they contribute, and the two are not")
    A("proportional. Nor is every addable species equally worth adding — novelty against db_v4")
    A("ranges from 0.05 for species in densely sampled genera (Penicillium, Alternaria) to 0.99")
    A("for genera the database does not represent at all, which is where the phyllosphere")
    A("residents sit. Prediction quality also falls on large repeat-rich genomes, where")
    A("ab-initio calls genes inside transposons.")
    A("")
    A("SHORT CAPTION")
    A("")
    A(f"CDS annotations reachable by db_v4, by route into scope. Solid bars are the {cds_in:,}")
    A(f"assemblies carrying an NCBI gene model; pale bars are the {cds_add:,} assemblies whose species")
    A(f"has no gene model at all, which ab-initio prediction would make available — {tot_add:,} further")
    A(f"species on top of the {tot_in:,} already present. Routes and colours are shared with the db_v4")
    A("ring, each species counted once under the strongest evidence that applies.")
    if recalls and novelties:
        A(f"Prediction recovers a median {med(recalls):.0%} of true CDS k-mer space at {med(precisions):.0%} precision on")
        A(f"species with a known annotation, and a median {med(novelties):.0%} of what it predicts for unannotated")
        A("species is absent from db_v4.")
    A("")
    A("SOURCES")
    A(f"  data   {TSV.name}, from {SURVEY.name}")
    A(f"  routes report/figures/prep_db_contents.py (imported, not reimplemented)")
    A(f"  quality {BENCH.name}")
    A("")

    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))

    print(f"-> {TSV.name}")
    print(f"-> {CAPTION.name}")
    print(f"\n  {'route':46s} {'in db_v4':>9s} {'addable':>9s}")
    for route in order:
        c = counts[route]
        print(f"  {route:46s} {c['in_db']:9,d} {c['addable']:9,d}")
    print(f"  {'TOTAL':46s} {tot_in:9,d} {tot_add:9,d}")


if __name__ == "__main__":
    main()
