#!/usr/bin/env python3
"""declared_worklist.py — what each BioProject needs before `declared_pathogens` can be set.

## Why this field exists

03_kraken can say "*Zymoseptoria* is present in this wheat sample". The chapter's claim needs
one more thing: "and the study that generated it never reported *Zymoseptoria*". That second
half is a property of the STUDY, not of the reads, so it has to be curated from the paper like
every other fact in this module.

It is also what makes the two evidence standards work. A detection is a SECOND pathogen where
one was already declared, and a FIRST pathogen where none was. Both are findings; they are
different findings, and the declared set is the only thing that tells them apart.

## Read relative to OUR scope, not in general

A study that targeted a VIRUS declares no fungus. Since viruses are out of scope, its declared
set is empty for our purposes and every fungus found in it is unreported. Two cohort projects
are like this (PRJNA796936 `Riboviria`, PRJNA1232637 `Wheat streak mosaic virus`), and reading
their declared set as "non-empty" would wrongly demote their detections.

## This script does NOT extract anything

It builds the worklist and stops. The value is filled in by reading the paper and quoting it,
the same way `setting` and `location` were. In particular it deliberately does **not** import
`llm_named_pathogens` from the old 02_literature frame. Leon's call 2026-10-15, and the reason
is stronger than that frame's 84.3% precision: a draft answer in front of a reader anchors
them into confirming rather than extracting, and an imported value would be the only thing in
this module with no quote behind it.

## The archive is the EXPECTATION, never the evidence

SRA's own `organism` field turns out to do most of the triage for free, and independently of
the LLM: 12 of the 14 disease-selected projects name their target pathogen there, and all 20
unselected projects name the plant. So `expectation` below says what to expect before reading.
It is a triage signal in the cheap-to-expensive cascade this module uses everywhere, and it
carries no evidential weight: the archive has misled this project in nine documented ways, and
`geo_loc_name` describing the isolate rather than the plant is the same class of error.

Run from the repository root:

    python 01_literature/curate/declared_worklist.py
    python 01_literature/curate/declared_worklist.py --summary   # counts only, write nothing
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def _root() -> Path:
    here = Path(__file__).resolve()
    for d in (here, *here.parents):
        if d.is_dir() and (d / "CLAUDE.md").exists():
            return d
    raise SystemExit("repo root not found")


ROOT = _root()
sys.path.insert(0, str(ROOT / "01_literature"))
from _layout import CURATION, GOLD, MODULE, SILVER, TEXT_CACHE  # noqa: E402

RELEASE = MODULE / "release"
BIOPROJECTS = RELEASE / "cohort_bioprojects.tsv"
SAMPLES = RELEASE / "cohort_samples.tsv"
OUT = SILVER / "declared_worklist.tsv"

FIELD = "declared_pathogens"

# Viruses are out of scope, so a virus-targeted study declares nothing WE can detect.
VIRUS_HINTS = ("riboviria", "virus", "viroid", "viridae", "virales")


def cached_fulltext_dois() -> set[str]:
    """DOIs whose full text is already on disk, so the paper needs no fetching."""
    out = set()
    if not TEXT_CACHE.exists():
        return out
    with open(TEXT_CACHE) as fh:
        for line in fh:
            _, _, rest = line.partition("\t")
            try:
                d = json.loads(rest)
            except json.JSONDecodeError:
                continue
            if (d.get("full_text") or "").strip():
                out.add((d.get("doi") or "").strip().lower())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", action="store_true", help="counts only, write nothing")
    a = ap.parse_args()

    for p in (BIOPROJECTS, SAMPLES):
        if not p.exists():
            raise SystemExit(f"missing {p}; run report/build_release.py first")

    projects = list(csv.DictReader(open(BIOPROJECTS), delimiter="\t"))
    samples = list(csv.DictReader(open(SAMPLES), delimiter="\t"))

    selection, crops, hosts = (defaultdict(Counter) for _ in range(3))
    for r in samples:
        selection[r["BioProject"]][r["sampling_selection"]] += 1
        crops[r["BioProject"]][r["crop"]] += 1
        hosts[r["BioProject"]][r["host"]] += 1

    # "Does the archive name a plant?" is answered against the CURATED host column rather than
    # a hand-written genus list, so it cannot drift as the cohort changes.
    plant_genera = {h.split()[0] for c in hosts.values() for h in c if h}

    cached = cached_fulltext_dois()

    # Rules already written for this field, so a rerun does not re-ask for finished work.
    existing = defaultdict(list)
    if CURATION.exists():
        for r in csv.DictReader(open(CURATION), delimiter="\t"):
            if (r.get("field") or "").strip() == FIELD:
                existing[r["BioProject"]].append(r.get("rule_id", ""))

    rows = []
    for p in projects:
        bp = p["BioProject"]
        org = (p.get("organism") or "").strip()
        sel = selection[bp].most_common(1)[0][0] if selection[bp] else ""
        low = org.lower()
        if any(h in low for h in VIRUS_HINTS):
            # Non-plant AND non-fungal: targeted a virus, so declares nothing in scope.
            expectation = "expect_empty_virus_target"
        elif org and org.split()[0] not in plant_genera:
            expectation = "expect_pathogen_named_in_archive"
        elif sel == "unselected":
            expectation = "expect_empty"
        else:
            expectation = "unknown_paper_required"

        dois = [d.strip() for d in (p.get("dois") or "").replace(";", ",").split(",")
                if d.strip()]
        n_cached = sum(1 for d in dois if d.lower() in cached)
        rows.append({
            "BioProject": bp,
            "runs": int(p.get("runs") or 0),
            "crop": crops[bp].most_common(1)[0][0] if crops[bp] else "",
            "sampling_selection": sel,
            "archive_organism": org,
            "expectation": expectation,
            "n_papers": p.get("papers") or "",
            "dois": ";".join(dois),
            "fulltext": f"{n_cached}/{len(dois)}" if dois else "0/0",
            "rule_state": ("done:" + ";".join(existing[bp])) if existing[bp] else "TODO",
            # Filled in by hand, from the paper. Left empty on purpose.
            "declared_pathogens": "",
            "evidence_id": "",
            "note": "",
        })

    # Ordered by runs, so the expensive attention goes where the data is. Same rule as
    # scaffold.py's worklist.
    rows.sort(key=lambda r: -r["runs"])

    print(f"{len(rows)} BioProjects, {sum(r['runs'] for r in rows):,} runs\n")
    print(f"{'expectation':34s} {'projects':>8s} {'runs':>7s} {'fulltext ready':>15s}")
    for exp, _ in Counter(r["expectation"] for r in rows).most_common():
        sel = [r for r in rows if r["expectation"] == exp]
        ready = sum(1 for r in sel if r["fulltext"].split("/")[0] != "0")
        print(f"{exp:34s} {len(sel):8d} {sum(r['runs'] for r in sel):7,} "
              f"{ready:>8d}/{len(sel)}")
    todo = [r for r in rows if r["rule_state"] == "TODO"]
    print(f"\nTODO: {len(todo)} projects, {sum(r['runs'] for r in todo):,} runs")

    if a.summary:
        print("\n--summary: nothing written")
        return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    print("Fill `declared_pathogens` from the PAPER, add a quote to seeds/provenance.tsv, then "
          "write the rule into seeds/curation.tsv. See curate/CURATING.md.")


if __name__ == "__main__":
    main()
