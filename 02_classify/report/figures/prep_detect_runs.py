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
import random
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import CLASSIFIED, COINFECTION  # noqa: E402
from _util import atomic_write  # noqa: E402

TSV = HERE / "detect_runs.tsv"
SEED = 20261020
# One sample per rung of the outcome, so the panel row reads as a gradient from clean to
# heavily co-infected rather than as an arbitrary sample.
WANT = [(0, 1), (1, 2), (2, 2), (3, 1), (5, 1), (8, 1)]   # (n_pathogens, how many runs)


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
        elif r["call"] == "yes" and (r["kingdom"] or "").strip() == "Viridiplantae":
            declared = "host"
        elif r["call"] == "yes":
            declared = "other fungus"
        else:
            declared = ""
        rows.append({
            "run": r["run"], "BioProject": r["BioProject"], "crop": s.get("crop", ""),
            "host": s.get("host", ""), "selection": s.get("sampling_selection", ""),
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
