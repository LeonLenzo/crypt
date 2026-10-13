#!/usr/bin/env python3
"""build_seeds.py — regenerate the GENERATED half of 03_kraken/seeds/ from holobase.

seeds/ is the single source of truth for what is in and what is out of the CDS database.
Half of it is AUTHORED and must never be overwritten by a script (sources.tsv, hosts.tsv,
exclusions.tsv, residents.tsv); half is a reproducible extraction from holobase. This script
writes only the second half. See README.md for which is which and why.

THE SCOPE RULE, and it is the thing most likely to be broken by a one-line change:

    Plant-associated, ANY host. NEVER scope by the cohort's own hosts.

USDA Fungus-Host Distributions is host-ASSOCIATION data. It is used here to ask "is this
fungus ever recorded on any plant", which is a scope gate. It must never be used to ask "is
it recorded on wheat/maize/rice/sorghum", which is a plausibility filter, because the study
exists to find UNREPORTED host associations and such a filter makes them undetectable by
construction. CLAUDE.md already rejects host-plausibility filtering of CALLS for this reason
(Leon's objection). Moving the same filter earlier, into reference selection, is worse: a
filter on calls can be relaxed afterwards, but a missing reference is invisible and
unrecoverable without a rebuild, because the reads go silently to the nearest taxon that IS
present. That failure mode is documented in 01_search/README.md.

Two traps in holobase, both of which returned confidently wrong answers before being caught:

  * `taxon.taxon_id` is an INTERNAL SURROGATE, not an NCBI taxid — the schema says so, and
    the real key is `taxon.ncbi_taxid`. Joining db_v3's taxids against taxon_id returns a
    confident 0% overlap instead of an error. Always join on ncbi_taxid.
  * `assembly.is_annotated` is 1 on all 36,709 rows, so holobase CANNOT tell you whether an
    assembly has a gene model. The database is CDS-only, so annotation is the binding
    constraint, and the only way to size a build is a real --select-only run. Availability
    counts overestimated the resident set by 3x (508 assemblies against the selector's 161).

Also: `taxon.kingdom` spells plants `Viridiplantae`, not `Plants`. Filtering on 'Plants'
silently returns zero.

Run from the repository root:

    python 03_kraken/seeds/build_seeds.py            # -> fungi_usda.tsv
    python 03_kraken/seeds/build_seeds.py --check    # verify, write nothing
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from pathlib import Path


def _root() -> Path:
    here = Path(__file__).resolve()
    for d in (here, *here.parents):
        if d.is_dir() and (d / "CLAUDE.md").exists():
            return d
    raise SystemExit("repo root not found")


ROOT = _root()
sys.path.insert(0, str(ROOT))
from _paths import DB_INSPECT, REFDATA  # noqa: E402

HERE = Path(__file__).resolve().parent
HOLOBASE = REFDATA / "holobase.db"
USDA_TSV = HERE / "fungi_usda.tsv"

# Assembly levels a reference may be built from. Contig-only taxa are excluded here because
# 01_search/kraken_search.py's own gate is scaffold_plus(), so listing them as seeds would
# promise taxa the selector then silently drops. If that gate is ever relaxed (there is an
# open decision about it for six resident genera whose only annotated assembly is a contig,
# including Papiliotrema, rice's most abundant genus), relax it here too and in step.
SCAFFOLD_PLUS = ("Scaffold", "Chromosome", "Complete Genome")

# Ranks worth seeding. formae speciales are kept deliberately: the rusts' f.sp. split is a
# real signal in 05_filter, and `_expand_taxids()` must pass intermediate_nodes=True for the
# same reason. Dropping them would lose 38 taxa including wheat rust lineages.
RANKS = ("species", "forma specialis", "varietas", "forma")


def db_v3_fungal_species() -> set[str]:
    """NCBI taxids of db_v3's fungal species nodes, from its own inspect table."""
    out, stack = set(), []
    with open(DB_INSPECT) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            raw = p[5]
            depth = (len(raw) - len(raw.lstrip(" "))) // 2
            while stack and stack[-1][0] >= depth:
                stack.pop()
            anc = {t for _, t in stack}
            stack.append((depth, p[4].strip()))
            if "4751" in anc and p[3].strip() == "S":
                out.add(p[4].strip())
    return out


QUERY = f"""
SELECT tf.ncbi_taxid                          AS ncbi_taxid,
       tf.scientific_name                     AS scientific_name,
       tf.rank                                AS rank,
       COUNT(DISTINCT th.ncbi_taxid)          AS n_plant_hosts,
       (SELECT COUNT(*) FROM assembly a
          WHERE a.taxon_id = tf.taxon_id
            AND a.assembly_level IN {SCAFFOLD_PLUS}) AS n_assemblies
  FROM host_association h
  JOIN taxon tf ON tf.taxon_id = h.fungus_id
  JOIN taxon th ON th.taxon_id = h.host_id
 WHERE h.source = 'usda_fungus_host'
   AND tf.ncbi_taxid IS NOT NULL
   AND th.kingdom = 'Viridiplantae'          -- NOT 'Plants'; that returns zero
   AND tf.rank IN {RANKS}
 GROUP BY tf.ncbi_taxid
HAVING n_assemblies > 0
 ORDER BY n_assemblies DESC, n_plant_hosts DESC
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    a = ap.parse_args()

    if not HOLOBASE.exists():
        raise SystemExit(f"holobase not found at {HOLOBASE}")
    con = sqlite3.connect(f"file:{HOLOBASE}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    rows = [dict(r) for r in con.execute(QUERY)]
    con.close()

    v3 = db_v3_fungal_species()
    for r in rows:
        r["in_db_v3"] = "yes" if str(r["ncbi_taxid"]) in v3 else "no"

    new = [r for r in rows if r["in_db_v3"] == "no"]
    print(f"USDA plant-associated fungi with a scaffold+ assembly : {len(rows):,}")
    print(f"  already in db_v3                                   : {len(rows) - len(new):,}")
    print(f"  NEW to db_v3                                       : {len(new):,}")
    print(f"  assemblies they hold (pre-floor, pre-annotation)    : "
          f"{sum(r['n_assemblies'] for r in new):,}")
    print("\nNOTE: n_assemblies ignores annotation, which holobase cannot report "
          "(is_annotated is 1 on\nevery row). The database is CDS-only, so the real "
          "addition needs a --select-only run.")

    if a.check:
        print("\n--check: nothing written")
        return
    with open(USDA_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "ncbi_taxid", "scientific_name", "rank", "n_plant_hosts",
            "n_assemblies", "in_db_v3"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {USDA_TSV.relative_to(ROOT)}  ({len(rows):,} rows)")


if __name__ == "__main__":
    main()
