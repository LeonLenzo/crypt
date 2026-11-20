#!/usr/bin/env python3
"""Read a Kraken2 report, and turn a taxon line into a coverage and a saturation figure.

Shared by classify.py, calibrate_floors.py and detect.py. Ported from
research/02_classify/_reports.py rather than imported, so a release module runs on its own.

## The two report layouts, told apart by COLUMN COUNT and never assumed

    6 cols   pct  reads  taxReads  rank  taxid  name
    8 cols   pct  reads  taxReads  minimizers  distinct  rank  taxid  name

The 8-column form is what `--report-minimizer-data` produces. Indexing blind reads the
minimizer count as the rank code, which does not error: it yields rank="1234567" on every
line, nothing matches rank "S", and EVERY SPECIES ROW DISAPPEARS while the script exits 0.
The count is therefore checked per line.

## The denominators

    coverage   = distinct minimizers observed / the taxon's CLADE minimizers in the build
    saturation = coverage / (1 - exp(-redundancy))

CLADE and not DIRECT. Both arms contain zero sub-species nodes, so the two columns are equal
for every species in them and the choice is currently a no-op. It stays `clade` because a
future build that admits a strain-level assembly would silently break a `direct` reader: the
species node would read 0 while its minimizers sat under the strain.

Measured against the DATABASE and not the genome, so reference incompleteness largely cancels.
What does not cancel is that the ceiling itself varies: in db1p the species ceilings run from
2,921 to 70,294,921 clade minimizers, so a relative floor is 2,308 different floors. That is
why `absolute_floor_ok` is applied ALONGSIDE the relative floor and not instead of it.
"""

import csv
import math
from pathlib import Path

# Kraken2 report lines carry very long taxon names. Raised once, on import.
csv.field_size_limit(10 ** 7)

SPECIES_RANK = "S"
KRAKEN_K = 35


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
    """taxid -> {'clade', 'direct', 'rank', 'name'} from kraken2-inspect output.

    Read from the arm that produced the reports. The two arms ship the same filename pattern
    under 02_reference/data/outputs/, and a script holding both must keep them apart by arm.
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


def kmers_per_spot(base_count, read_count, layout: str, k: int = KRAKEN_K) -> float | None:
    """K-mers one spot contributes, from the library's own bases per spot.

    It must come from the LIBRARY and never from a constant: at k=35 a 51 nt mate yields 17
    k-mers and a 150 nt mate 116, a factor of 7 on the saturation expectation, and the cohort
    holds both. ENA's base_count/read_count is bases per SPOT, so a paired library's per-mate
    length is half of it and both mates contribute.
    """
    try:
        bases, reads = float(base_count), float(read_count)
    except (TypeError, ValueError):
        return None
    if reads <= 0 or bases <= 0:
        return None
    per_spot = bases / reads
    mates = 2 if (layout or "").upper().startswith("PAIR") else 1
    per_mate = per_spot / mates
    return mates * max(0.0, per_mate - k + 1)


def coverage(distinct, taxid: int, inspect: dict) -> float | None:
    """Distinct minimizers observed over the taxon's clade total. None if either is unusable."""
    if distinct is None:
        return None
    e = inspect.get(taxid)
    if not e or not e["clade"]:
        return None
    return distinct / e["clade"]


def redundancy(reads: int, taxid: int, inspect: dict, kmers_per_spot: float) -> float | None:
    """How many times over the reads assigned to a taxon could cover its k-mer space."""
    e = inspect.get(taxid)
    if not e or not e["clade"] or not reads or not kmers_per_spot:
        return None
    return reads * kmers_per_spot / e["clade"]


def saturation(distinct, reads: int, taxid: int, inspect: dict,
               kmers_per_spot: float) -> float | None:
    """Distinct minimizers observed over what uniform sampling at that depth would give.

        expected = M * (1 - exp(-N * kmers_per_spot / M))    M = the taxon's clade total

    Normalising both axes by M makes the null ONE curve for every taxon, so a 2.9 Mb
    Colletotrichum and 33 Mb wheat are comparable on one panel.

    Real calls land at 0.07 to 0.56 rather than near 1.0, because RNA-seq expression is not
    uniform and reads pile onto highly expressed genes. The discriminating fact is that
    pile-ups sit two further orders of magnitude below.

    Do not substitute reads/distinct. It is the same intuition without the correction and it
    does not work: measured across four runs its separation decayed 20.8x, 7.5x, 3.1x, 1.0x as
    depth rose, because a deeply sequenced real organism accumulates redundant reads too.
    """
    if distinct is None or not reads:
        return None
    e = inspect.get(taxid)
    if not e or not e["clade"]:
        return None
    m = e["clade"]
    # Clamped because exp(-x) underflows for a deeply sequenced taxon, where the expectation
    # is simply M: every minimizer should have been seen.
    expected = m * (1.0 - math.exp(-min(700.0, reads * kmers_per_spot / m)))
    return distinct / expected if expected else None


def absolute_floor_ok(distinct, floor: int) -> bool:
    """The companion to the relative floor, and not optional.

    A relative floor alone is not a floor: db2's smallest species ceiling is 897 minimizers, so
    1% of it is 9 k-mers, which cross-mapping supplies for free. Written POSITIVELY, because
    `distinct` can be None and `None < floor` raises while `not (None >= floor)` would admit it.
    """
    return distinct is not None and distinct >= floor
