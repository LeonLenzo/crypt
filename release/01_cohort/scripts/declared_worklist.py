#!/usr/bin/env python3
"""declared_worklist.py — what is left before every cohort project has a declaration.

    python release/01_cohort/scripts/declared_worklist.py
    python release/01_cohort/scripts/declared_worklist.py --summary   # counts, write nothing

Reads  data/curation/declared_pathogens.tsv   the field being completed
       data/curation/relations.tsv            WHO GENERATED each BioProject
       data/curation/screen_fulltext.tsv      the full-text screen's own extraction
       data/raw/fulltext.tsv                  which papers we hold text for, and where
       data/outputs/supp_data_s2_cohort.tsv   runs per project, and the analysed set
       ../03_classify/data/outputs/supp_data_s5_detections.tsv   what rests on each project
Writes data/curation/declared_worklist.tsv    one row per project still to curate

## Why the field exists

`03_classify` can say *Zymoseptoria* is present in a run. The chapter needs the second half,
"and the study that generated that run never reported it", which is a property of the STUDY
and has to be read out of the paper like every other curated fact. 18 of the 31 cohort
projects have no declaration, so their detections read `unknown`, and `detect.py` keeps them
out of the undeclared count rather than letting 344 runs manufacture the headline.

## THIS SCRIPT EXTRACTS NOTHING, DELIBERATELY

It assembles what a curator needs in front of them and stops. `screen_target_pathogens` is
carried as an EXPECTATION, never as an answer: it is the full-text screen's regex-and-rules
extraction, not a curated reading, and a draft answer in front of a reader anchors them into
confirming rather than extracting. The value is set by reading the paper and quoting it, with
`evidence_doi`, `evidence_locus` and `evidence_quote` filled in together.

## The generating paper, not any paper that cites the accession

A BioProject can appear in several papers and usually does; PRJNA749614 is named by three.
Only one GENERATED the runs the cohort took, and that is the study whose declaration counts.
Curating against a reanalysis would record what somebody else looked for.

**`relations.tsv` IS AN EXCEPTION TABLE, AND ABSENCE FROM IT MEANS UNAMBIGUOUS.** It holds
the cases the deposit-section heuristic could not settle, so it carries 41 rows for 31
projects and most projects are simply not in it. An earlier version of this script treated a
missing `generated` row as "unresolved" and reported five projects as blocked that are not:
three are named by exactly one paper, and two carry a `reanalysed` row, which is POSITIVE
information because it eliminates the only competing candidate. The resolution is therefore:

    an explicit `generated` row                       -> that paper
    otherwise, S1 papers naming the project, less
    any marked `reanalysed` or `undecided` for it      -> if one remains, that paper
                                                       -> if several, genuinely ambiguous

## Read the chain before fetching anything

The last pass on this field lost real time to looking for each project's own paper in
`studies/` by DOI and treating a missing directory as a missing paper. The determination was
already written in the evidence chain the project carried for `setting` and `location`, from a
paper already read. `has_text` below says whether the body is on disk; `n_curated_runs` says
whether this project already has curated facts, which means it has a read paper behind it.
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]              # release/01_cohort
DATA = HERE / "data"
S5 = HERE.parent / "03_classify/data/outputs/supp_data_s5_detections.tsv"
OUT = DATA / "curation" / "declared_worklist.tsv"
csv.field_size_limit(10 ** 7)

COLS = ["BioProject", "status", "n_runs", "n_analysed_runs", "n_runs_with_wp_call",
        "n_wp_taxa", "generating_doi", "resolved_by", "other_dois", "has_text",
        "text_source",
        "screen_target_pathogens", "screen_mentioned_pathogens", "n_curated_runs",
        "affects_figures", "wp_taxa_at_stake"]


def rd(path, comment="#"):
    with open(path) as fh:
        rows = [ln for ln in fh if not ln.startswith(comment)]
    return list(csv.DictReader(rows, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--summary", action="store_true", help="print counts, write nothing")
    args = ap.parse_args()

    declared = rd(DATA / "curation/declared_pathogens.tsv")
    cohort = rd(DATA / "outputs/supp_data_s2_cohort.tsv")
    relations = rd(DATA / "curation/relations.tsv")
    screen = {r["doi"]: r for r in rd(DATA / "curation/screen_fulltext.tsv")}
    text = {r["doi"]: r for r in rd(DATA / "raw/fulltext.tsv")}
    curated = rd(DATA / "curation/curated_attributes.tsv")

    runs, analysed = defaultdict(set), defaultdict(set)
    for r in cohort:
        runs[r["BioProject"]].add(r["run"])
        if r["setting"] == "field" and r["host"] == "Triticum aestivum":
            analysed[r["BioProject"]].add(r["run"])

    curated_runs = {r["run"] for r in curated}
    n_curated = defaultdict(int)
    for r in cohort:
        if r["run"] in curated_runs:
            n_curated[r["BioProject"]] += 1

    wp_runs, wp_taxa = defaultdict(set), defaultdict(set)
    for r in rd(S5):
        if r["is_wheat_pathogen"] == "yes":
            wp_runs[r["BioProject"]].add(r["run"])
            wp_taxa[r["BioProject"]].add(r["taxon"])

    # the generating paper, and separately every paper that merely names the accession.
    # DOIs ARE CASE-INSENSITIVE and the tables disagree on case: relations.tsv carries
    # 10.1094/MPMI-09-21-0236-A where the screen carries the lowercase form, so a set
    # difference on the raw strings listed a project's own generating paper as a competing
    # one. Compare folded, report the spelling the source used.
    gen = {r["BioProject"]: r["doi"] for r in relations if r["relation"] == "generated"}
    eliminated = defaultdict(set)
    for r in relations:
        if r["relation"] in ("reanalysed", "undecided"):
            eliminated[r["BioProject"]].add(r["doi"].lower())
    mentions = defaultdict(dict)
    for r in screen.values():
        for b in filter(None, (x.strip() for x in r["bioprojects"].split(";"))):
            mentions[b][r["doi"].lower()] = r["doi"]

    def generating(p):
        """-> (doi, how). See the module docstring: absence from relations.tsv is not a gap."""
        if p in gen:
            return gen[p], "explicit"
        left = {k: v for k, v in mentions[p].items() if k not in eliminated[p]}
        if len(left) == 1:
            return next(iter(left.values())), "sole candidate"
        if not left:
            return "", "NO CANDIDATE"
        return "", f"AMBIGUOUS ({len(left)})"

    todo = [r for r in declared if r["status"] == "not_curated"]
    rows = []
    for r in sorted(todo, key=lambda x: -len(runs[x["BioProject"]])):
        p = r["BioProject"]
        g, how = generating(p)
        t = text.get(g, {}) or text.get(g.lower(), {})
        sc = screen.get(g, {}) or screen.get(g.lower(), {})
        rows.append({
            "BioProject": p, "status": r["status"],
            "n_runs": len(runs[p]), "n_analysed_runs": len(analysed[p]),
            "n_runs_with_wp_call": len(wp_runs[p]), "n_wp_taxa": len(wp_taxa[p]),
            "generating_doi": g or "UNRESOLVED", "resolved_by": how,
            "other_dois": ";".join(v for k, v in sorted(mentions[p].items())
                                   if k != g.lower()),
            "has_text": "yes" if t.get("path") else "no",
            "text_source": t.get("source", ""),
            "screen_target_pathogens": sc.get("target_pathogens", ""),
            "screen_mentioned_pathogens": sc.get("mentioned_pathogens", ""),
            "n_curated_runs": n_curated[p],
            # curating only changes a figure where the project has runs in the analysed set
            # AND a wheat-pathogen call in them; elsewhere it buys completeness, not a result
            "affects_figures": "yes" if (analysed[p] & wp_runs[p]) else "no",
            "wp_taxa_at_stake": ";".join(sorted(wp_taxa[p])),
        })

    done = [r for r in declared if r["status"] != "not_curated"]
    print(f"{len(done)} of {len(declared)} projects have a declaration "
          f"({sum(1 for r in done if r['status'] == 'none_declared')} of them 'none').")
    print(f"{len(rows)} to curate: {sum(r['n_runs'] for r in rows)} runs, "
          f"{sum(r['n_analysed_runs'] for r in rows)} of them in the analysed set, "
          f"{sum(r['n_runs_with_wp_call'] for r in rows)} with a wheat-pathogen call.")
    hot = [r for r in rows if r["affects_figures"] == "yes"]
    print(f"  {len(hot)} of them have analysed runs carrying a wheat-pathogen call, so "
          f"curating those moves a figure; the other {len(rows) - len(hot)} buy completeness: "
          + ", ".join(r["BioProject"] for r in hot))
    byhow = defaultdict(list)
    for r in rows:
        byhow[r["resolved_by"]].append(r["BioProject"])
    print("  generating paper: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(byhow.items())))
    no_gen = [r["BioProject"] for r in rows if r["generating_doi"] == "UNRESOLVED"]
    no_txt = [r["BioProject"] for r in rows if r["has_text"] == "no"]
    if no_gen:
        print(f"  GENERATING PAPER UNRESOLVED ({len(no_gen)}): {', '.join(no_gen)}")
        print(f"    -> resolve the relation first; curating against a reanalysis would "
              f"record what somebody else looked for.")
    if no_txt:
        print(f"  NO BODY TEXT on disk ({len(no_txt)}): {', '.join(no_txt)}")
        print(f"    -> check the project's existing evidence chain BEFORE fetching: a "
              f"cohort project is defined by manuscript-level evidence, so the paper has "
              f"already been read for some other field.")

    if args.summary:
        return
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"  -> {OUT}")


if __name__ == "__main__":
    main()
