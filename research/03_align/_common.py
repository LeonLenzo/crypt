#!/usr/bin/env python3
"""_common.py — the tables research/03_align reads, named once.

Everything here is produced by `release/`. **This module consumes the deliverable and never
reaches around it into another module's scratch**, which is the same contract the figures run
on: the shipped tables are the interface, so a number quoted here can be traced to a published
file rather than to an intermediate that may have moved.

Paths that cross module boundaries in the research tree live in `research/_paths.py`; these do
not cross anything, they are this module's own inputs.
"""
from pathlib import Path

_MARKERS = ("CLAUDE.md", ".gitignore")


def _root(start=None):
    here = Path(start or __file__).resolve()
    for d in (here, *here.parents):
        if d.is_dir() and any((d / m).exists() for m in _MARKERS):
            return d
    raise RuntimeError(f"repo root not found above {here}")


ROOT = _root()
RELEASE = ROOT / "release"

# ── 03_classify: what was called, and in which runs ─────────────────────────────────────────
DETECTIONS = RELEASE / "03_classify/data/outputs/supp_data_s5_detections.tsv"
RUNS = RELEASE / "03_classify/data/outputs/supp_data_s6_runs.tsv"

# ── 02_reference: the taxa, the assemblies, and which assembly was already chosen ───────────
REFERENCE = RELEASE / "02_reference/data/outputs/supp_data_s3_reference.tsv"
ASSEMBLIES = RELEASE / "02_reference/data/curation/screen_assemblies.tsv"
BUILD_PLAN = RELEASE / "02_reference/data/curation/build_plan.tsv"
EXTRA_CDS = RELEASE / "02_reference/data/curation/extra_cds_db2.tsv"
HOST_CDS_MANIFEST = RELEASE / "02_reference/data/raw/host_cds_manifest.tsv"

# ── 01_cohort: the runs themselves ──────────────────────────────────────────────────────────
COHORT = RELEASE / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"

# CDS as 02_reference stores it. One directory per accession on scratch; the same tree is in
# the db1p library backup on Acacia, which is where it comes from if scratch has been purged.
CDS_DIR = RELEASE / "02_reference/data/raw/cds"
HOST_CDS_DIR = RELEASE / "02_reference/data/raw/host_cds"

# The host, decided 2026-11-25: Chinese Spring alone, NOT db2's 19-cultivar sink.
HOST_ENSEMBL = "triticum_aestivum_refseqv2"
HOST_ACCESSION = "GCA_018294505.1"
HOST_TAXID = "4565"

ARM = "db2"
