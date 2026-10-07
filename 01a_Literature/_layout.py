#!/usr/bin/env python3
"""Every path in 01a_Literature, and which layer may be regenerated.

Added 2026-10-06 after three incidents in one day with a single cause: nothing in the layout
distinguished a file that may be rebuilt from one that must never be overwritten. `papers.tsv`
was hand-annotated and the next `triage.py` run discarded it; an imported accession list
survived only as a shell argument; and the file written to fix that was gitignored, because
tracking depended on `git add -f` habits rather than on what kind of file it was.

    seeds/            HAND-WRITTEN and irreplaceable. The module's whole intellectual asset: every
            curated value traces to a rule, a join and a quote in here. TRACKED by default.
            Never written by a script except by appending a reviewed row.
    data/bronze/      immutable captures of somebody else's data: runinfo, BioSample attributes, GEO
            brief records, ENA reports, and the papers under STUDIES. Deleting one costs a
            refetch and nothing else. Ignored.
    data/silver/      generated tables: runs, bioprojects, papers, paper_bioproject, registry,
            worklist. Rebuildable from BRONZE + SEEDS, so ignored. Hand-editing one is always
            a mistake.
    data/gold/        small human-readable outputs. Generated, but tracked anyway: these are the summary
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

# Two top-level ideas: seeds/ is AUTHORED, data/ is DERIVED. Everything derived is a
# medallion layer under data/, each named for its layer, so the middle one is no longer the
# only directory named after nothing. studies/ stays outside both: it is a per-study human
# drop zone for papers and supplements, not a pipeline layer.
SEEDS = MODULE / "seeds"
DATA = MODULE / "data"
BRONZE = DATA / "bronze"
SILVER = DATA / "silver"
GOLD = DATA / "gold"
STUDIES = MODULE / "studies"

# --- bronze: external captures, one file per accession ----------------------------
RUNINFO = BRONZE / "runinfo"
BIOSAMPLE = BRONZE / "biosample_attrs"
GEO = BRONZE / "geo"
BIOPROJECT_XML = BRONZE / "bioproject"
ENA = BRONZE / "ena"

# --- silver: generated, rebuildable ------------------------------------------------
RUNS = SILVER / "runs.tsv"
BIOPROJECTS = SILVER / "bioprojects.tsv"
PAPERS = SILVER / "papers.tsv"
PAPER_BIOPROJECT = SILVER / "paper_bioproject.tsv"
REGISTRY = SILVER / "registry.tsv"
WORKLIST = SILVER / "worklist.tsv"
ACCESSION_PAPERS = SILVER / "accession_papers.tsv"
RESOLVED_ACCESSIONS = SILVER / "resolved_accessions.tsv"
COHORT_STUDIES = SILVER / "kraken_cohort_studies.tsv"
RUN_SPECIES = SILVER / "run_species.tsv"
TRUNCATION = SILVER / "truncation.tsv"

# --- seeds: hand-written, irreplaceable ------------------------------------------
# The Undermind search returns. Authored elsewhere and NOT reproducible - re-running a
# search months later returns a different set - so they are seeds, not bronze. They are also
# the provenance for 171 candidate accessions.
UNDERMIND = SEEDS / "undermind"

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
