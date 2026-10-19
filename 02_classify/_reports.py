#!/usr/bin/env python3
"""_reports.py — read a Kraken2 report, and turn a taxon line into a coverage figure.

Module-level support, used by `04_classify`, `05_detect` and the rarefaction utility. Ported
from `03_kraken/04_assign/kraken_report_analysis.py` rather than imported from it, because
importing would make the replacement module depend on the module it replaces, and the whole
point of the rewrite was to cut those live dependencies.

## The two report layouts, told apart by COLUMN COUNT and never assumed

    6 cols   pct  reads  taxReads  rank  taxid  name
    8 cols   pct  reads  taxReads  minimizers  distinct  rank  taxid  name

The 8-column form is what `--report-minimizer-data` produces. Indexing blind reads the
minimizer count as the rank code, which does not error: it silently yields rank="1234567" for
every line, so nothing matches rank "S" and EVERY SPECIES DISAPPEARS. That is why the count is
checked per line.

## Coverage, and what its denominator is

    coverage = distinct minimizers observed for a taxon
               ---------------------------------------
               that taxon's CLADE minimizers in the database, from kraken2-inspect

Measured against the DATABASE, not the genome, so reference incompleteness largely cancels: a
species represented by one assembly has a small ceiling and a read that hits it still scores
proportionally. What does NOT cancel is that the ceiling itself varies enormously. Measured on
db_v4's 1,547 species:

    min           2,922 minimizers  ->  1% is       29 k-mers
    median    4,049,049             ->  1% is   40,490
    max      70,426,087             ->  1% is  704,261

So a 1% relative floor is 1,547 different floors spanning 24,000x, and 20 species have ceilings
under 100,000 where 1% is fewer than 1,000 k-mers -- reachable by cross-mapping noise alone.
`absolute_floor_ok()` exists for that reason and must be applied ALONGSIDE the relative floor,
not instead of it.

CLADE and not DIRECT minimizers, deliberately. db_v4 tags sequence at mixed ranks (a species
node for most assemblies, a strain node where NCBI resolved one), so a species' own k-mers are
split between its node and its children. The direct column would count only the former and
under-report every species that has a strain in the database. Both sides of the published db_v3
ratios were clade figures and remain comparable.
"""

import csv
import math
from pathlib import Path

# Kraken2 report lines carry very long taxon names; every consumer of parse_report hit this and
# raised the limit by hand. Done once, on import.
csv.field_size_limit(10 ** 7)

SPECIES_RANK = "S"


def parse_report(path: Path):
    """Yield one dict per taxon line, with its lineage from the name column's indentation."""
    stack = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            mini = len(parts) >= 8
            raw_name = parts[7] if mini else parts[5]
            depth = (len(raw_name) - len(raw_name.lstrip(" "))) // 2
            try:
                rec = {
                    "reads_clade": int(parts[1]),
                    "reads_direct": int(parts[2]),
                    "minimizers": int(parts[3]) if mini and parts[3].strip().isdigit() else None,
                    "distinct": int(parts[4]) if mini and parts[4].strip().isdigit() else None,
                    "rank": (parts[5] if mini else parts[3]).strip(),
                    "taxid": int((parts[6] if mini else parts[4]).strip()),
                    "name": raw_name.strip(),
                    "depth": depth,
                }
            except ValueError:
                continue
            while stack and stack[-1]["depth"] >= depth:
                stack.pop()
            rec["lineage"] = [(s["rank"], s["taxid"], s["name"]) for s in stack]
            stack.append(rec)
            yield rec


def load_inspect(path: Path) -> dict:
    """taxid -> {'clade': n, 'direct': n, 'rank': r, 'name': s} from kraken2-inspect output.

    This is the coverage denominator for every figure downstream, so it is read from the db that
    actually produced the reports. db_v3's inspect table still exists at the repo root as
    `_paths.DB_INSPECT` and has the SAME filename; a script needing both must alias one.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                out[int(f[4])] = {"clade": int(f[1]), "direct": int(f[2]),
                                  "rank": f[3].strip(), "name": f[5].strip()}
            except ValueError:
                continue
    return out


def coverage(distinct: int | None, taxid: int, inspect: dict) -> float | None:
    """Distinct minimizers observed over the taxon's clade total. None if either is unusable."""
    if distinct is None:
        return None
    e = inspect.get(taxid)
    if not e or not e["clade"]:
        return None
    return distinct / e["clade"]


def saturation(distinct: int | None, reads: int, taxid: int, inspect: dict,
               kmers_per_read: float) -> float | None:
    """How much of a taxon's k-mer space the reads assigned to it actually reached.

        saturation = distinct observed / distinct EXPECTED if those reads were spread
                     uniformly over the taxon's k-mer space
        expected   = M * (1 - exp(-N * kmers_per_read / M))      M = the taxon's clade total

    THIS IS THE ANSWER TO "HOW DO I KNOW THE CALL IS REAL", and it is the cheap form of what
    `detect.py` fits as its accumulation-curve `score`. Leon asked it on 2026-10-19 about
    99,509 reads assigned to *Oryza sativa* in a wheat sample: if that were real rice, those
    reads should have produced ~6.6 M distinct minimizers. They produced 24,388 -- 270x fewer.
    Real rice cannot look like that; only a few hammered loci can.

    The cause is wheat's reference being INCOMPLETE, not the classifier being unreliable. A
    read from a wheat transcript absent from db_v4's wheat CDS, whose rice ortholog is present,
    has only rice k-mers to match and gets the best home that exists. Were both present the LCA
    would resolve it up to Poaceae.

    Measured over 4 runs of PRJEB33109 (63 taxa):

        clears 1% coverage          n= 9   median 0.2795   range 0.0700 - 0.5595
        leakage (>=1k reads, <1%)   n=54   median 0.0017   range 0.0001 - 0.0457

    No overlap, and the separation is stable per run (241x, 359x, 230x, 162x) BECAUSE the
    expectation scales with read count. Do not substitute reads/distinct, which is the same
    intuition without the correction and is useless: across those same four runs its
    separation decayed 20.8x -> 7.5x -> 3.1x -> 1.0x as depth rose, because a deeply sequenced
    REAL organism also accumulates redundant reads. A 3.0 cut-off on it flagged 56% of leakage
    and 22% of genuine calls.

    Real calls land at 0.07-0.56 rather than near 1.0, and that is expected: RNA-seq expression
    is wildly non-uniform, so reads pile onto highly expressed genes. The discriminating fact is
    that artefacts sit two further orders of magnitude below.

    `kmers_per_read` is (L-K+1) summed over the mates of one spot, so 2*(150-34) for a 2x150
    library. It must come from the LIBRARY, not a constant: at k=35 a 51 nt read yields 17
    k-mers and a 150 nt read 116, which moves the expectation by a factor of 7.

    CAVEAT, unresolved as of 2026-10-19: this was measured on four wheat runs of ONE project,
    all Puccinia-dominated. Re-test across the full cohort before any threshold is set on it.
    """
    if distinct is None or not reads:
        return None
    e = inspect.get(taxid)
    if not e or not e["clade"]:
        return None
    m = e["clade"]
    # Clamped because exp(-x) underflows for a deeply sequenced taxon and the expectation is
    # simply M there: every minimizer should have been seen.
    expected = m * (1.0 - math.exp(-min(700.0, reads * kmers_per_read / m)))
    return distinct / expected if expected else None


def absolute_floor_ok(distinct: int | None, floor: int) -> bool:
    """The companion to the relative floor, and not optional.

    A relative floor alone is not a floor: on db_v4's smallest species 1% is 29 k-mers, which
    cross-mapping supplies for free. Test POSITIVELY -- `distinct` can be None and None < floor
    raises while `not (None >= floor)` would admit it.
    """
    return distinct is not None and distinct >= floor
