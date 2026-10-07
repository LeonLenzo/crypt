#!/usr/bin/env python3
"""Every path in 01a_Literature, and which layer may be regenerated.

Added 2026-10-06 after three incidents in one day with a single cause: nothing in the layout
distinguished a file that may be rebuilt from one that must never be overwritten. `papers.tsv`
was hand-annotated and the next `triage.py` run discarded it; an imported accession list
survived only as a shell argument; and the file written to fix that was gitignored, because
tracking depended on `git add -f` habits rather than on what kind of file it was.

    SEEDS   HAND-WRITTEN and irreplaceable. The module's whole intellectual asset: every
            curated value traces to a rule, a join and a quote in here. TRACKED by default.
            Never written by a script except by appending a reviewed row.
    BRONZE  immutable captures of somebody else's data: runinfo, BioSample attributes, GEO
            brief records, ENA reports, and the papers under STUDIES. Deleting one costs a
            refetch and nothing else. Ignored.
    DATA    generated tables (silver): runs, bioprojects, papers, paper_bioproject, registry,
            worklist. Rebuildable from BRONZE + SEEDS, so ignored. Hand-editing one is always
            a mistake.
    GOLD    small human-readable outputs. Generated, but tracked anyway: these are the summary
            numbers a reader wants to see change over the life of the project.

The rule in one line: **if losing a file costs a refetch it is BRONZE or DATA; if it costs
re-reading papers it is SEEDS.**

Named `_layout` and NOT `_paths` because the repo root already owns `_paths.py` for
cross-module paths. A second module of that name shadows it for whichever script imports
first, which immediately broke review.py: it resolved `_paths` to the root module and lost
SEEDS. Nine scripts here import ROOT from the root module, so the two must stay distinct.

Scripts must take paths from here rather than deriving them from `Path(__file__).parent`.
That was the actual debt: ~50 hardcoded `HERE / "data"` expressions meant the layer split
existed in the directory listing but not in the code, and it made moving any script unsafe.
"""

from pathlib import Path

MODULE = Path(__file__).resolve().parent      # 01a_Literature/
ROOT = MODULE.parent                          # phd/01-review/

SEEDS = MODULE / "seeds"
BRONZE = MODULE / "bronze"
DATA = MODULE / "data"
GOLD = MODULE / "gold"
STUDIES = MODULE / "studies"

# --- bronze: external captures, one file per accession ----------------------------
RUNINFO = BRONZE / "runinfo"
BIOSAMPLE = BRONZE / "biosample_attrs"
GEO = BRONZE / "geo"
ENA = BRONZE / "ena"

# --- data: generated, rebuildable ------------------------------------------------
RUNS = DATA / "runs.tsv"
BIOPROJECTS = DATA / "bioprojects.tsv"
PAPERS = DATA / "papers.tsv"
PAPER_BIOPROJECT = DATA / "paper_bioproject.tsv"
REGISTRY = DATA / "registry.tsv"
WORKLIST = DATA / "worklist.tsv"
ACCESSION_PAPERS = DATA / "accession_papers.tsv"
RESOLVED_ACCESSIONS = DATA / "resolved_accessions.tsv"
COHORT_STUDIES = DATA / "kraken_cohort_studies.tsv"
RUN_SPECIES = DATA / "run_species.tsv"
TRUNCATION = DATA / "truncation.tsv"

# --- seeds: hand-written, irreplaceable ------------------------------------------
CURATION = SEEDS / "curation.tsv"
JOINS = SEEDS / "joins.tsv"
PROVENANCE = SEEDS / "provenance.tsv"
DECISIONS = SEEDS / "decisions.tsv"
FOUND = SEEDS / "found.tsv"
RUN_SCOPE = SEEDS / "run_scope.tsv"
BACKFILL = SEEDS / "backfill_accessions.tsv"

# --- gold: generated, small, tracked ---------------------------------------------
COHORT_TSV = GOLD / "cohort.tsv"
HOST_SUMMARY = GOLD / "host_summary.tsv"
OFFTARGET = GOLD / "offtarget.tsv"
