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

# ── 01_stat: NCBI's k-mer profiles ───────────────────────────────────────────────────────────
STAT_FILTERED  = ROOT / "01_stat/03_filter/data"

# ── 02_literature: the papers and the sample universe ───────────────────────────────────────
SAMPLES        = ROOT / "02_literature/03_classify/data/samples.tsv"

# ── 03_kraken/01_search: reference selection and the CDS tree ───────────────────────────────
CDS_PATHOGEN   = ROOT / "03_kraken/01_search/data/cds/pathogen"
CDS_HOST       = ROOT / "03_kraken/01_search/data/cds/host"
REF_CANDIDATES = ROOT / "03_kraken/01_search/data/ref_candidates_v3.tsv"
HOST_CANDIDATES = ROOT / "03_kraken/01_search/data/host_candidates.tsv"

# ── 03_kraken/02_build: the Kraken2 database ────────────────────────────────────────────────
DB_INSPECT     = ROOT / "03_kraken/02_build/data/db_v3_inspect.tsv"

# ── 03_kraken/03_select: the reads ──────────────────────────────────────────────────────────
READS_DIR      = ROOT / "03_kraken/03_select/data/reads"
READS_MANIFEST = ROOT / "03_kraken/03_select/data/manifest.tsv"
RUN_LIST       = ROOT / "03_kraken/03_select/data/run_list.tsv"

# ── 03_kraken/04_assign: per-run classification ─────────────────────────────────────────────
REPORTS        = ROOT / "03_kraken/04_assign/data/reports"
HOST_CALLS     = ROOT / "03_kraken/04_assign/host/data/host_calls.tsv"
DETECTIONS     = ROOT / "03_kraken/04_assign/data/detections.tsv"

# ── 03_kraken/05_filter: the presence calls this pipeline turns on ──────────────────────────
CLASSIFIED     = ROOT / "03_kraken/05_filter/data/_classified.tsv"
PATHO_LINEAGE  = ROOT / "03_kraken/05_filter/data/patho_lineage.tsv"
HOST_LINEAGE   = ROOT / "03_kraken/05_filter/data/host_lineage.tsv"
FILTER_DATA    = ROOT / "03_kraken/05_filter/data"

# ── 03_kraken/utilities: QC records, not pipeline steps ─────────────────────────────────────
BUSCO_SCORES   = ROOT / "03_kraken/utilities/busco/data/busco_scores.tsv"
SKETCH_INDEX   = ROOT / "03_kraken/utilities/saturation/data/sketch_index.tsv"
SKETCHES       = ROOT / "03_kraken/utilities/saturation/data/sketches"
SATURATION     = ROOT / "03_kraken/utilities/saturation/data/saturation_summary.tsv"

# ── 04_kallisto: competitive read assignment ────────────────────────────────────────────────
STRATA         = ROOT / "04_kallisto/01_select/data/strata.tsv"
NEIGHBOURS     = ROOT / "04_kallisto/02_build/data/neighbours.tsv"
KALLISTO_IDX   = ROOT / "04_kallisto/02_build/data/idx"
QUANT_DATA     = ROOT / "04_kallisto/03_quant/data"
CALLS          = ROOT / "04_kallisto/04_confirm/data/calls.tsv"
