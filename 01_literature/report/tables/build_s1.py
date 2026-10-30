#!/usr/bin/env python3
"""build_s1.py -- Supplementary Data S1: every publication the literature screen assessed.

    python 01_literature/report/tables/build_s1.py

WHY THIS TABLE EXISTS. Leon, 2026-10-29: the opening sentence of Results carries five counts
in one breath -- publications assessed, candidate manuscripts, manuscripts with sufficient
methodological resolution, BioProjects, runs and biological samples -- and a reader has no way
to check any of them. Each is now DERIVED from this table rather than asserted beside it, so a
number in the text and a number in the supplement cannot disagree: if they do, one of them was
hand-edited.

RESTRICTED TO MANUSCRIPTS SHOWN TO CONTAIN THE FOUR HOSTS, as the union of two independent
determinations, because neither alone is sound: the archive's own BioProject `organism` field
names one of the four cereals, or our curation established a cereal host for its samples.

THE ARCHIVE FIELD ALONE IS NOT ENOUGH, measured: it catches 27 of the 41 cohort projects and
misses 14, every one a rust survey registered under Puccinia striiformis, P. graminis or
P. triticina rather than under the wheat it was sampled from. Filtering on it would silently
drop the studies this paper most depends on. The curation alone is not enough either: it
exists only for projects that already reached the cohort, so using it as the candidate filter
collapses the funnel onto its own endpoint.

ONE ROW PER CANDIDATE MANUSCRIPT, and nothing else. Leon, 2026-10-29: S1 is the output of the
TITLE AND ABSTRACT screen, so it carries the manuscripts that screen admitted and the
BioProjects they reference. What happened to each project afterwards is S2's job, and the raw
yield of the search before screening is not reported at all.

COLUMNS A REVIEWER CAN USE, and no others. The build previously carried paper_key,
discovery_source, cohort_bioprojects, assessed, cereal_candidate, in_cohort and
undermind_decision. Every one of those is internal bookkeeping -- which tool found a paper,
which stage flag it tripped -- and none of it helps anyone check a number. A supplementary
table is not a dump of the working state.

THE SCREEN WAS NOT CEREAL-ONLY. The 421 triaged BioProjects include Vitis vinifera (23),
Arabidopsis thaliana (12), Solanum lycopersicum (10) and Malus domestica (9), so a sentence
saying the search was for cereals misdescribes what was done. Cereals are the curated subset.
"""

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _util import atomic_write  # noqa: E402

SILVER = MODULE / "data" / "silver"
RELEASE = MODULE / "release"
OUT = HERE / "Supplementary Data S1.tsv"

COLS = ["doi", "title", "n_bioprojects", "bioprojects"]


def clean_title(t: str) -> str:
    t = re.sub(r"\*+", "", t)
    t = t.replace("\t", " ").replace("\n", " ")
    return re.sub(r"\s+", " ", t).strip()


def main() -> None:
    papers = list(csv.DictReader(open(SILVER / "papers.tsv"), delimiter="\t"))
    pb = list(csv.DictReader(open(SILVER / "paper_bioproject.tsv"), delimiter="\t"))
    cohort_bp = {r["BioProject"] for r in
                 csv.DictReader(open(RELEASE / "cohort_bioprojects.tsv"), delimiter="\t")}
    # "Shown to contain the four hosts" is the UNION of two determinations; see the module
    # docstring for why neither is sufficient alone.
    HOST_GENERA = ("Triticum", "Zea", "Oryza", "Sorghum")
    org_cereal = {r["BioProject"] for r in
                  csv.DictReader(open(SILVER / "bioprojects.tsv"), delimiter="\t")
                  if any(g in (r["organism"] or "") for g in HOST_GENERA)}
    host_bp = cohort_bp | org_cereal
    samples = list(csv.DictReader(open(RELEASE / "cohort_samples.tsv"), delimiter="\t"))
    # cohort_samples.tsv carries no BioSample column; the Run -> BioSample map lives in
    # silver/runs.tsv. Counting Runs as BioSamples makes 7,688 of a figure that is 6,222.
    run_bs = {r["Run"]: r["BioSample"] for r in
              csv.DictReader(open(SILVER / "runs.tsv"), delimiter="\t")}

    # Per cohort BioProject: crops, runs and distinct BioSamples actually retained.
    bp_crop, bp_runs, bp_bs = {}, {}, {}
    for s in samples:
        b = s["BioProject"]
        bp_crop.setdefault(b, set()).add(s.get("crop", ""))
        bp_runs[b] = bp_runs.get(b, 0) + 1
        bp_bs.setdefault(b, set()).add(run_bs.get(s["Run"], s["Run"]))

    # paper -> BioProjects, from the explicit link table AND the paper's own column, because
    # the two disagree for umbrella projects and the union is what "this paper names" means.
    # paper_bioproject.paper_ref holds BOTH forms: 323 of 457 are the DOI that papers.tsv
    # calls paper_key, 134 are the short key it calls paper_ref. Joining on either alone
    # silently drops the other, which cost 14 papers from the cohort count.
    by_paper = {}
    for r in pb:
        by_paper.setdefault(r["paper_ref"], set()).add(r["BioProject"])

    rows, inco_papers, inco_bp = [], [], set()
    for p in papers:
        key = p["paper_key"]
        named = {x.strip() for x in (p["bioprojects"] or "").split(";") if x.strip()}
        named |= by_paper.get(key, set()) | by_paper.get(p.get("paper_ref", ""), set())
        inco = sorted(named & cohort_bp)
        crops = sorted({c for b in inco for c in bp_crop.get(b, set()) if c})
        # Only the host-containing projects are listed, not every accession a manuscript
        # happens to cite, so the column and the count describe the same set.
        hits = sorted(named & host_bp)
        if not hits:
            continue
        rows.append({
            "doi": p.get("doi", "") or key,
            # Titles arrive carrying the reference manager's markdown emphasis
            # (*et al.*, **bold**). Harmless internally, wrong in a reviewer-facing file.
            "title": clean_title(p.get("title") or ""),
            "n_bioprojects": len(hits),
            "bioprojects": ";".join(hits),
        })
        if inco:
            inco_papers.append(key)
            inco_bp.update(inco)
    rows.sort(key=lambda r: (-r["n_bioprojects"], r["doi"]))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with atomic_write(OUT) as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    allbp = {x for r in rows for x in r["bioprojects"].split(";") if x}
    print(f"-> {OUT.name}  ({len(COLS)} columns)")
    print(f"   candidate manuscripts                  {len(rows):,}")
    print(f"   BioProjects they reference             {len(allbp):,}")
    print(f"   ... that reached the curated cohort    {len(inco_bp):,}")
    print(f"   cohort runs                            {len(samples):,}")
    print(f"   cohort biosamples                      "
          f"{len({run_bs.get(s['Run'], s['Run']) for s in samples}):,}")


if __name__ == "__main__":
    main()
