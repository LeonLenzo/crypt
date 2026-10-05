#!/usr/bin/env python3
"""Shared paths and domain logic for the 03_kraken steps.

Unnumbered, so support rather than a pipeline step. Composes the root `_util.py` (domain-free
helpers) and `_paths.py` (cross-module artefacts); what lives here is 03_kraken's own vocabulary.

The reason this file exists is `detections()`. The predicate that defines a pathogen detection
was written out longhand in six scripts:

    rank == "S", distinct is not None, taxid in pathogens, taxid not in hosts,
    reads_clade >= 100, distinct >= 1

That is a scientific definition, not boilerplate. Changing the 100-read floor meant six
coordinated edits, and missing one would have produced a figure that disagreed with its own
table without failing. Now the floor is `MIN_READS_CLADE` here and nowhere else.

Two scripts deliberately do NOT use it and should not be made to:
  03_kraken/05_filter/_scratch/prep_detection_panels.py keeps non-pathogens as their own class
  03_kraken/04_assign/host/kraken_host_call.py selects hosts, which is the opposite question
"""
import csv
import functools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _paths import DB_INSPECT, HOST_CALLS, REPORTS, ROOT  # noqa: E402

# Kraken2 report lines can carry very long taxon names; every consumer of parse_report() hit
# this and raised the limit by hand. Done once, on import.
csv.field_size_limit(10 ** 7)

# The detection predicate, in one place.
MIN_READS_CLADE = 100   # below this the accumulation shape cannot be judged
MIN_DISTINCT    = 1     # a species with no distinct minimizers is not evidence
SPECIES_RANK    = "S"


def _report_api():
    """The 04_assign report API, imported without a relative path hack.

    Every caller did `sys.path.insert(0, '03_kraken/04_assign')` with a bare relative string,
    which only resolved when the process started at the repo root. Anchored on _paths.ROOT here
    so it works from any depth and any working directory.
    """
    sys.path.insert(0, str(ROOT / "03_kraken/04_assign"))
    from kraken_report_analysis import load_inspect, load_reference, parse_report
    return load_reference, parse_report, load_inspect


def load_inspect(path=None):
    """db_v3's per-taxon minimizer ceilings. Re-exported so callers need no path hack."""
    _, _, _load_inspect = _report_api()
    return _load_inspect(path or DB_INSPECT)


@functools.lru_cache(maxsize=1)
def reference():
    """(pathogen taxids, host taxids) as sets of str, loaded once per process."""
    _load_reference, _, _ = _report_api()
    return _load_reference()


@functools.lru_cache(maxsize=1)
def host_of():
    """run accession -> host species, from 04_assign's host calls."""
    return {r: h for r, (h, _) in run_meta().items()}


@functools.lru_cache(maxsize=1)
def run_meta():
    """run accession -> (host, biosample). Both come from the same file, so read it once."""
    out = {}
    with open(HOST_CALLS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run = r["run"].strip()
            out[run] = ((r.get("host") or "").strip(), (r.get("biosample") or "").strip())
    return out


def passes(rec, pathogens, hosts, min_reads=MIN_READS_CLADE, min_distinct=MIN_DISTINCT):
    """The pathogen-detection predicate, as a function so callers that need the surrounding
    records (detect.py sums reads over sub-taxa) can apply it without reimplementing it."""
    if rec["rank"] != SPECIES_RANK or rec["distinct"] is None:
        return False
    tid = str(rec["taxid"])
    if tid in hosts or tid not in pathogens:
        return False
    return rec["reads_clade"] >= min_reads and rec["distinct"] >= min_distinct


def reports(reports_dir=None, host_fallback=""):
    """Yield (run, host, recs) per Kraken2 report, recs being every parsed line.

    Sorted, so anything built on it is reproducible. Callers needing only the passing records
    should use detections(); this exists for callers that need a record's siblings too.
    """
    _, parse_report, _ = _report_api()
    meta = run_meta()
    for p in sorted(Path(reports_dir or REPORTS).glob("*.txt")):
        host = meta.get(p.stem, ("", ""))[0] or host_fallback
        yield p.stem, host, list(parse_report(p))


def detections(reports_dir=None, min_reads=MIN_READS_CLADE, min_distinct=MIN_DISTINCT,
               require_host=False, host_fallback=""):
    """Yield (run, host, rec) for every pathogen detection passing the standard predicate."""
    pathogens, hosts = reference()
    for run, host, recs in reports(reports_dir, host_fallback):
        if require_host and not host:
            continue
        for rec in recs:
            if passes(rec, pathogens, hosts, min_reads, min_distinct):
                yield run, host, rec
