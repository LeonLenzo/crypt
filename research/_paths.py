#!/usr/bin/env python3
"""_paths.py — the artefacts that cross module boundaries, named once.

Scope, deliberately narrow: this file holds ONLY paths that one module produces and another
consumes. Those are the pipeline's interfaces, and they are what breaks when the tree is
renumbered. Paths used inside a single module belong in that module's `_common.py`; putting
them here would mean a root file had to know every module's internal layout.

Two problems this fixes, both observed in the tree rather than imagined:

1. **Bare relative paths.** Several scripts open `'03_kraken/05_filter/data/three.tsv'` with no
   root anchor, so they only work when invoked from the repo root and silently write output to
   the wrong place from anywhere else. Everything here is absolute.

2. **`parents[N]` with a hand-counted N.** Scripts reach the repo root with
   `Path(__file__).resolve().parents[3]`, where N depends on how deeply the script is nested.
   The 2026-10-01 renumbering changed nesting for some scripts and the counts had to be chased.
   `ROOT` below finds the root by marker instead, so a script can move between levels and keep
   working.

Usage, from any depth:

    import sys; sys.path.insert(0, str(_find_root()))   # or see <module>/_common.py
    from _paths import CLASSIFIED, REPORTS, HOST_CALLS

Related: `_util.py` holds domain-free helpers, `manifest.py` is a CLI over one of them, and
`<module>/_common.py` holds a module's own paths and domain logic.
"""
from pathlib import Path

# Markers that identify the repo root. CLAUDE.md and .gitignore are both tracked and both sit
# only at the root, so either one is sufficient and two means a missing file does not break it.
_ROOT_MARKERS = ("CLAUDE.md", ".gitignore")


def _find_root(start=None):
    """Walk up from `start` (default: this file) until a directory holding a marker is found."""
    here = Path(start or __file__).resolve()
    for d in (here, *here.parents):
        if d.is_dir() and any((d / m).exists() for m in _ROOT_MARKERS):
            return d
    # a checkout without the markers is a broken checkout, not something to paper over
    raise RuntimeError(f"repo root not found above {here}; expected one of {_ROOT_MARKERS}")


ROOT = _find_root()

# ══ RETIRED MODULES, moved to superseded/ on 2026-10-20 ══════════════════════════════════════
#
# 01_stat, 02_literature, 03_kraken and 04_kallisto are retired and untracked; the directories
# stay on disk under superseded/. See SUPERSEDED.md for what replaced each one.
#
# The constants below are re-rooted rather than deleted, because four of them are still READ by
# live code and the rest are how a retired output is located if it is ever needed again:
#
#   STAT_FILTERED, SAMPLES   01_literature/retrieve/triage.py, report/truncation_check.py,
#                            curate/import_provenance.py -- the curation audit chain
#   REF_CANDIDATES           02_classify/01_refs/compare_refs.py, as the db_v3 baseline
#   CLASSIFIED               02_classify/report/figures/prep_rarefaction.py, as the db_v3 arm
#
# ONE prefix, so a second move needs one edit here and none anywhere else. That is the whole
# reason this file exists; the 2026-10-01 renumbering is what proved it.
SUPERSEDED     = ROOT / "superseded"

# ── 01_stat: NCBI's k-mer profiles. Retired; the STAT frame is not the cohort frame ─────────
STAT_FILTERED  = SUPERSEDED / "01_stat/03_filter/data"

# ── 02_literature: the papers and the sample universe. Retired; 01_literature replaced it ──
SAMPLES        = SUPERSEDED / "02_literature/03_classify/data/samples.tsv"

# ── 03_kraken/01_search: reference selection and the CDS tree. 02_classify/01_refs replaced ──
CDS_PATHOGEN   = SUPERSEDED / "03_kraken/01_search/data/cds/pathogen"
CDS_HOST       = SUPERSEDED / "03_kraken/01_search/data/cds/host"
REF_CANDIDATES = SUPERSEDED / "03_kraken/01_search/data/ref_candidates_v3.tsv"
HOST_CANDIDATES = SUPERSEDED / "03_kraken/01_search/data/host_candidates.tsv"

# ── 03_kraken/02_build: db_v3. 02_classify/02_build built db_v4 ─────────────────────────────
DB_INSPECT     = SUPERSEDED / "03_kraken/02_build/data/db_v3_inspect.tsv"

# ── 03_kraken/03_select: the reads. 02_classify/03_cohort owns them now ─────────────────────
# READS_DIR is the FLAT pre-migration tree and still holds the ~9.1 TB of non-cohort
# STAT-gated runs on scratch, which are meant to expire. The live cohort tree is
# `02_classify/_layout.READS`, one directory per BioProject.
READS_DIR      = SUPERSEDED / "03_kraken/03_select/data/reads"
READS_MANIFEST = SUPERSEDED / "03_kraken/03_select/data/manifest.tsv"
RUN_LIST       = SUPERSEDED / "03_kraken/03_select/data/run_list.tsv"

# ── 03_kraken/04_assign: per-run classification. 02_classify/04_classify replaced it ───────
REPORTS        = SUPERSEDED / "03_kraken/04_assign/data/reports"
HOST_CALLS     = SUPERSEDED / "03_kraken/04_assign/host/data/host_calls.tsv"
DETECTIONS     = SUPERSEDED / "03_kraken/04_assign/data/detections.tsv"

# ── 03_kraken/05_filter: the db_v3 presence calls. 02_classify/05_detect replaced it ───────
# CLASSIFIED here is the db_v3 table. Its `score`/`cls` columns are the SUPERSEDED statistic
# (a fitted accumulation curve, no coverage floor at cls == "real"); 05_detect's own
# _classified.tsv carries `call`/`reason` against saturation instead. Coverage is NOT
# comparable across the two, because the denominator comes from each build's inspect table.
CLASSIFIED     = SUPERSEDED / "03_kraken/05_filter/data/_classified.tsv"
PATHO_LINEAGE  = SUPERSEDED / "03_kraken/05_filter/data/patho_lineage.tsv"
HOST_LINEAGE   = SUPERSEDED / "03_kraken/05_filter/data/host_lineage.tsv"
FILTER_DATA    = SUPERSEDED / "03_kraken/05_filter/data"

# ── 03_kraken/utilities: QC records, not pipeline steps ─────────────────────────────────────
BUSCO_SCORES   = SUPERSEDED / "03_kraken/utilities/busco/data/busco_scores.tsv"
SKETCH_INDEX   = SUPERSEDED / "03_kraken/utilities/saturation/data/sketch_index.tsv"
SKETCHES       = SUPERSEDED / "03_kraken/utilities/saturation/data/sketches"
SATURATION     = SUPERSEDED / "03_kraken/utilities/saturation/data/saturation_summary.tsv"

# ── 04_kallisto: competitive read assignment. Superseded by the planned 02_classify/04_align ─
STRATA         = SUPERSEDED / "04_kallisto/01_select/data/strata.tsv"
NEIGHBOURS     = SUPERSEDED / "04_kallisto/02_build/data/neighbours.tsv"
KALLISTO_IDX   = SUPERSEDED / "04_kallisto/02_build/data/idx"
QUANT_DATA     = SUPERSEDED / "04_kallisto/03_quant/data"
CALLS          = SUPERSEDED / "04_kallisto/04_confirm/data/calls.tsv"

# External reference data, outside the repository. See _refdata/README.md and the per-database
# PROVENANCE.md. Only core_nt's inspect table is held locally; the 316 GB index is not.
REFDATA        = ROOT.parents[1] / "_refdata"
CORENT_INSPECT = REFDATA / "databases/kraken2_core_nt_20251015/inspect.txt"

# 01_literature's studies/ drop zone holds published supplements. Colvin & Chopra 2026's
# per-sample fungal abundance tables (doi:10.64898/2026.05.21.727028) are the empirical
# resident-community seed for reference scope, so they cross from 01_literature into 03_kraken.
COLVIN_DATA    = ROOT / "01_literature/studies/doi_10.64898_2026.05.21.727028"
