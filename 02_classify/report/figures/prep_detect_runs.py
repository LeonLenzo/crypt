#!/usr/bin/env python3
"""prep_detect_runs.py — the decision, one sample at a time.

The aggregate surface shows how the criterion behaves over 300,061 candidates. This shows what
it DID to a handful of real samples, which is the only view in which "this sample carries three
pathogens" is checkable by eye.

Same axes as the saturation cloud -- redundancy against coverage, both normalised by the
taxon's own db_v4 minimizer total, so the null y = 1-exp(-x) is one curve for every taxon. What
is added is the CALL: which points 05_detect accepted, which it rejected, and on what grounds.

Runs are chosen to span the outcome rather than to flatter it: samples carrying zero, one, two
and many plant pathogens, drawn from different projects and crops where possible, and the
selection is seeded so the figure is reproducible.

    python 02_classify/report/figures/prep_detect_runs.py
    python 02_classify/report/figures/prep_detect_runs.py --runs ERR1607792,ERR3383632
"""

import argparse
import csv
import sqlite3
import random
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import (CLASSIFIED, COHORT_SAMPLES, COINFECTION,  # noqa: E402
                     HOLOBASE, LIT_RUNS)
from _util import atomic_write  # noqa: E402

TSV = HERE / "detect_runs.tsv"
SEED = 20261020
# One sample per rung of the outcome, so the panel row reads as a gradient from clean to
# heavily co-infected rather than as an arbitrary sample.
WANT = [(0, 1), (1, 2), (2, 2), (3, 1), (5, 1), (8, 1)]   # (n_pathogens, how many runs)


def plant_associated() -> set:
    """NCBI taxids with a USDA Fungus-Host record against a PLANT host.

    THE SAME DEFINITION ROUTE 4 USES, read from holobase rather than from the routes column,
    because that column is stale and cannot be refreshed cheaply. `select_refs.py` takes routes
    from `assemblies.tsv`, which is `survey.py`'s record of what NCBI held, and refreshing it
    means re-querying NCBI and pulling assemblies for 578 newly-seeded taxa -- a database
    composition change, not a relabel.

    The seed bug this works around is fixed in build_seeds.py and real: route 4 gated on
    assembly LEVEL, so a species was denied "recorded on a plant" because its only assembly was
    a contig, then picked up by route 1 on a genus-level FUNGuild label. Cladosporium allicinum
    holds 43 USDA plant-host records and was routed as a pathogen-genus taxon on that basis.
    Reading holobase asks the evidence directly and cannot go stale the same way.
    """
    if not HOLOBASE.exists():
        print(f"  no holobase at {HOLOBASE}; plant-association unknown", file=sys.stderr)
        return set()
    con = sqlite3.connect(HOLOBASE)
    return {int(t) for (t,) in con.execute(
        "SELECT DISTINCT tf.ncbi_taxid FROM host_association h "
        "JOIN taxon tf ON tf.taxon_id = h.fungus_id "
        "JOIN taxon th ON th.taxon_id = h.host_id "
        "WHERE h.source = 'usda_fungus_host' AND th.kingdom = 'Viridiplantae' "
        "AND tf.ncbi_taxid IS NOT NULL") if t is not None}


def guild_of(g: str) -> str:
    if not g:
        return "unlabelled"
    if "Plant Pathogen" in g:
        return "plant pathogen"
    if "Endophyte" in g:
        return "endophyte"
    if "Saprotroph" in g:
        return "saprotroph"
    if "Animal Pathogen" in g:
        return "animal pathogen"
    return "other"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", help="comma-separated run accessions, overriding the selection")
    args = ap.parse_args()

    summ = {s["run"]: s for s in csv.DictReader(open(COINFECTION), delimiter="\t")}
    # The run's own provenance, for the figure's title block. location, collection_date and
    # tissue are curated per run in the release; declared_pathogens lives in the literature
    # run table and is what the STUDY said it was sequencing, which is not the same as what
    # was detected and is the whole point of the comparison.
    PLANT_ASSOC = plant_associated()
    print(f"  {len(PLANT_ASSOC):,} taxids carry a USDA plant-host record")
    coh = {r["Run"]: r for r in csv.DictReader(open(COHORT_SAMPLES), delimiter="\t")}
    decl = {r["Run"]: (r.get("declared_pathogens") or "").strip()
            for r in csv.DictReader(open(LIT_RUNS), delimiter="\t")}

    if args.runs:
        chosen = [r.strip() for r in args.runs.split(",") if r.strip()]
    else:
        rng = random.Random(SEED)
        by_n = defaultdict(list)
        for s in summ.values():
            by_n[min(int(s["n_pathogens"]), 8)].append(s)
        chosen, used_proj = [], set()
        for n, k in WANT:
            pool = sorted(by_n.get(n, []), key=lambda s: s["run"])
            if not pool:
                continue
            rng.shuffle(pool)
            # Prefer a project not already shown, so six panels are not six views of one study.
            pool.sort(key=lambda s: s["BioProject"] in used_proj)
            for s in pool[:k]:
                chosen.append(s["run"])
                used_proj.add(s["BioProject"])

    rows = []
    for r in csv.DictReader(open(CLASSIFIED), delimiter="\t"):
        if r["run"] not in chosen:
            continue
        try:
            red, cov = float(r["redundancy"]), float(r["coverage"])
        except ValueError:
            continue
        if red <= 0 or cov <= 0:
            continue
        s = summ.get(r["run"], {})
        # DECLARED vs UNDECLARED, which is the paper's actual question and was not on the
        # figure. coinfection.tsv carries both lists per run; `undeclared` is a subset of
        # `pathogens`, so a called pathogen absent from `undeclared` is one the study named.
        # Anything that is not a called pathogen gets "" rather than "no", because "this study
        # did not declare the host plant" is not a meaningful statement.
        undeclared = {x for x in (s.get("undeclared") or "").split(";") if x}
        pathogens = {x for x in (s.get("pathogens") or "").split(";") if x}
        if r["call"] == "yes" and r["species"] in pathogens:
            declared = "undeclared" if r["species"] in undeclared else "declared"
        else:
            declared = ""
        # HOW THE TAXON ENTERED THE REFERENCE, in the Methods' own vocabulary, with PHI-BASE
        # AND FUNGUILD KEPT APART. Leon, 2026-10-27: those are different kinds of evidence and
        # the figure must not merge them. routes.py supplies the five membership routes; the
        # Methods then narrows "pathogen" to route 2 (PHI-base) or route 3 resolved at or below
        # species rank (FUNGuild), and splits the evidence three ways: documented on the host
        # sampled, documented on another plant, or a species-rank guild assignment with no host
        # record. `pathogen_source` is the seam: "phibase"/"both" against "funguild".
        route = (r.get("route") or "").strip()
        rts = {x for x in (r.get("routes") or "").split(";") if x}
        src = (r.get("pathogen_source") or "").strip()
        ph = (r.get("phibase_host") or "").strip()
        if (r["kingdom"] or "").strip() == "Viridiplantae":
            tier = "host"
        elif r.get("pathogen") == "yes" and ph == "this_host":
            tier = "phibase_crop"
        elif r.get("pathogen") == "yes" and ph == "other_plant":
            tier = "phibase_plant"
        elif r.get("pathogen") == "yes":
            tier = "funguild"          # species-rank guild, no PHI-base host record
        # NOT a reportable pathogen, so FALL BACK PAST `target`. The strongest route is the
        # right answer for membership but the wrong one for a detection label: "target" only
        # says a genus-level pathogen label exists somewhere, and precedence puts it above
        # every route that carries real species-level evidence. Asking what else is known
        # about the taxon is the honest question once the pathogen test has already failed.
        elif "host_associate" in rts:
            tier = "host_associate"
        elif "plant_associate" in rts or int(r["taxid"]) in PLANT_ASSOC:
            tier = "plant_associate"
        elif rts == {"congener"} or route == "congener":
            tier = "congener"
        else:
            # On the target route but NOT reportable: no PHI-base record and no species-rank
            # guild assignment. 2,148 of these 2,180 called rows carry an EMPTY guild string.
            # They are in scope to give the classifier competitors, not to be found.
            tier = "target_unresolved"
        rows.append({
            "tier": tier, "pathogen_source": src, "phibase_host": ph,
            "run": r["run"], "BioProject": r["BioProject"], "crop": s.get("crop", ""),
            "host": s.get("host", ""), "selection": s.get("sampling_selection", ""),
            "location": coh.get(r["run"], {}).get("location", ""),
            "year": (coh.get(r["run"], {}).get("collection_date", "") or "")[:4],
            "tissue": coh.get(r["run"], {}).get("tissue", ""),
            # "none" is a real curated value meaning the study declared nothing, not a gap.
            "declared_pathogens": "" if decl.get(r["run"], "") == "none"
                                  else decl.get(r["run"], ""),
            "n_pathogens": s.get("n_pathogens", ""),
            "n_undeclared": s.get("n_undeclared", ""), "declared": declared,
            "panel": f"{r['run']}  ·  {s.get('crop','?')}  ·  "
                     f"{s.get('n_pathogens','?')} pathogen"
                     f"{'' if s.get('n_pathogens')=='1' else 's'}",
            "species": r["species"], "guild": guild_of(r["guild"]),
            "kingdom": r["kingdom"], "reads": r["reads"], "distinct": r["distinct"],
            "redundancy": r["redundancy"], "coverage": r["coverage"],
            "saturation": r["saturation"], "call": r["call"], "reason": r["reason"],
        })

    # Panel order follows the pathogen count so the row reads clean -> heavily co-infected.
    order = sorted({(int(r["n_pathogens"]), r["run"], r["panel"]) for r in rows})
    rank = {run: i for i, (_, run, _) in enumerate(order)}
    for r in rows:
        r["panel_order"] = rank[r["run"]]

    with atomic_write(TSV) as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["panel_order"], -float(r["coverage"]))))

    print(f"-> {TSV.name}: {len(rows):,} candidates over {len(chosen)} runs")
    for n, run, panel in order:
        rs = [r for r in rows if r["run"] == run]
        called = [r for r in rs if r["call"] == "yes"]
        print(f"   {panel:52s} {len(rs):5d} candidates, {len(called):3d} called")


if __name__ == "__main__":
    main()
