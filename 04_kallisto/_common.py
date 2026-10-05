#!/usr/bin/env python3
"""Shared paths and helpers for the 04_kallisto steps.

Unnumbered, so support rather than a pipeline step. Every step inserts the module root on
sys.path via Path(__file__).resolve().parents[1] and imports from here, which is why the
numbered directories are safe: nothing imports by module name.
"""
import csv
import collections
import sys
from pathlib import Path

MODULE = Path(__file__).resolve().parent
ROOT   = MODULE.parent

# inputs from 03_kraken
CLASSIFIED = ROOT / "03_kraken/05_filter/data/_classified.tsv"
CDS_DIR    = ROOT / "03_kraken/01_search/data/cds/pathogen"
HOST_CDS   = ROOT / "03_kraken/01_search/data/cds/host"
REF_V3     = ROOT / "03_kraken/01_search/data/ref_candidates_v3.tsv"
HOST_V3    = ROOT / "03_kraken/01_search/data/host_candidates.tsv"
SKETCH_IX  = ROOT / "03_kraken/utilities/saturation/data/sketch_index.tsv"
SKETCHES   = ROOT / "03_kraken/utilities/saturation/data/sketches"
SATURATION = ROOT / "03_kraken/utilities/saturation/data/saturation_summary.tsv"
PATHO_LIN  = ROOT / "03_kraken/05_filter/data/patho_lineage.tsv"
READS_MAN  = ROOT / "03_kraken/03_select/data/manifest.tsv"
READS_ROOT = ROOT / "03_kraken/03_select/data"      # strata.tsv read_files are relative to this

# this module's own step outputs
NEIGHBOURS = MODULE / "02_build/data/neighbours.tsv"
IDX_DIR    = MODULE / "02_build/data/idx"
STRATA     = MODULE / "01_select/data/strata.tsv"
QUANT_DIR  = MODULE / "03_quant/data"
CALLS      = MODULE / "04_confirm/data/calls.tsv"

# the filter's presence call, which this module exists to calibrate
SCORE_MIN, COV_MIN = 0.7, 0.01


def slug(name):
    return name.lower().replace(" x ", "_x_").replace(" ", "_").replace("/", "_")


def sp2(name):
    """First two words of an organism name, i.e. the species."""
    p = name.split()
    return " ".join(p[:2]) if len(p) >= 2 else name


def num(row, key, default=0.0):
    try:
        return float(row[key])
    except (ValueError, KeyError, TypeError):
        return default


def classified(host=None):
    """Every row of the filter's output, optionally for one host."""
    rows = list(csv.DictReader(open(CLASSIFIED), delimiter="\t"))
    return [r for r in rows if r["host"] == host] if host else rows


def flagged(host=None):
    """Detections passing the filter's presence call."""
    return [r for r in classified(host)
            if num(r, "score") >= SCORE_MIN and num(r, "coverage") >= COV_MIN]


def flagged_hosts():
    """Hosts carrying at least one flagged detection. 01_build covers all of them, rather
    than only the ones 02_select happens to sample, so the steps stay independent."""
    return sorted({r["host"] for r in flagged()})


def busco_ranked():
    """species -> [accession] best-first on BUSCO completeness, then scaffold N50.

    Completeness rather than contiguity: a fragmented but complete CDS set still carries
    every gene a read can map to, which is all that matters for pseudoalignment.
    """
    out = collections.defaultdict(list)
    for r in csv.DictReader(open(SKETCH_IX), delimiter="\t"):
        out[sp2(r["organism_name"])].append(
            (-num(r, "complete_pct", -1.0), -num(r, "scaffold_n50_kb", -1.0), r["accession"]))
    return {s: [a for _, _, a in sorted(v)] for s, v in out.items()}


def step_log(step_file, name):
    """Start the project's log convention for a step, returning the _Tee to close.

    Wraps the three calls every other script in the repo makes by hand (make_log_dir, _Tee,
    link_latest) so a step here needs one line. This is what `_common.py` is for: it composes
    the root `_util.py` rather than re-implementing it, and keeps the repo-root sys.path insert
    in one place instead of once per script.

        log = step_log(__file__, "kallisto_select")
        try:
            ...
        finally:
            log.close()
    """
    sys.path.insert(0, str(ROOT))
    from _util import _Tee, link_latest, make_log_dir

    base = Path(step_file).resolve().parent / "logs"
    log_dir = make_log_dir(base)
    log = _Tee(log_dir / f"{name}.log")
    link_latest(base, log_dir / f"{name}.log")
    sys.stdout = log
    return log
