#!/usr/bin/env python3
"""survey.py — ask NCBI what assemblies and gene models ACTUALLY exist for the routes list.

    python 02_classify/01_refs/survey.py --limit 20     # smoke test first
    python 02_classify/01_refs/survey.py                # the real survey, ~1,400 queries

This is the measurement step, and it exists because **every estimate this project has made
from availability counts has been too high by roughly 3x.** holobase cannot answer the only
question that decides membership in a CDS-only database: `assembly.is_annotated` is 1 on all
36,709 rows, so it cannot filter for annotation at all. The two earlier sizings that were
checked against reality both collapsed:

    1,054 new USDA taxa "addable"   ->  179 once NCBI was asked
    508 resident assemblies         ->  161 once the selector actually ran

So nothing here is projected. One `datasets summary genome taxon <taxid>` per query unit,
aggregated per species.

## What it measures

Per species: how many assemblies exist, how many carry a gene model, how many are
scaffold-plus, the best assembly level, and the release-date span. Annotation is reported
separately at every level because it is the binding constraint rather than a quality detail:
measured over `fungi_annotated.tsv`, only **25% of scaffold-plus assemblies carry a gene
model**, and `kraken_build.py` skips any assembly without one, so an unannotated assembly
contributes exactly nothing to the index until miniprot has run.

A species with assemblies and zero annotation is therefore not a failure. It is the
**miniprot worklist**, and it is written out as one.

## Why it queries by taxid, and by the smallest covering clade

By TAXID because the `datasets` CLI refuses an ambiguous NAME and the clade then returns
nothing with no warning: *Cryptococcus* (fungus 5206 / scale insect 79213) and *Rhizophagus*
(fungus 1129544 / beetle 295919) cost 42 assemblies and 17 species that way.

By the smallest covering clade because one query returns a whole subtree. `routes.py` drops
any unit whose ancestor is also a unit, which takes ~4,100 members down to ~1,400 queries.
Querying both a genus and a species inside it is one wasted call and two chances to
double-count.

Every assembly's own taxid is normalised to its SPECIES node before counting. NCBI labels many
assemblies at strain, varietas or forma rank, and counting per raw taxid splits one species
across several nodes: that happened to 87 species in db_v3, worst *Zymoseptoria tritici* across
6 nodes.
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _layout import SILVER, LOGS, NCBI_SUMMARY                   # noqa: E402
from _util import _Tee, atomic_write, make_log_dir, link_latest  # noqa: E402
import _assembly as A                                            # noqa: E402
import _taxonomy                                                 # noqa: E402
import routes as R                                               # noqa: E402

SURVEY_TSV = SILVER / "assembly_survey.tsv"
MINIPROT_TSV = SILVER / "miniprot_worklist.tsv"
ASSEMBLIES_TSV = SILVER / "assemblies.tsv"
# One file per query unit, named by taxid. BRONZE, because it is an immutable capture of
# somebody else's data: deleting it costs a refetch and nothing else.
NCBI_CACHE = NCBI_SUMMARY / "datasets_summary"

SURVEY_COLS = [
    "ncbi_taxid", "name", "rank", "routes",
    "n_assemblies", "n_annotated", "n_scaffold_plus", "n_scaffold_plus_annotated",
    "best_level", "first_release", "last_release", "n_countries", "found_under",
]


def cached_query(ncbi: int, refresh: bool = False) -> tuple:
    """(assembly records, was_cached) for one taxon, caching NCBI's answer in BRONZE.

    The survey is the module's one expensive discovery step (~14 min over 1,321 clades) and it
    gets re-run whenever `seeds/` changes or a labelling bug is found. Without a cache, fixing
    a one-line aggregation bug costs another full pass against NCBI, which is both slow and
    rude. With it, the second run is seconds.

    It also makes the run RESUMABLE: a crash or a kill leaves every completed unit on disk.

    Written through `atomic_write`, so a crash mid-write cannot leave a truncated cache file
    that the next run reads as a complete answer. That failure mode is not hypothetical here:
    it took runs.tsv from 33,958 rows to 451 on 2026-10-15.
    """
    path = NCBI_CACHE / f"{ncbi}.jsonl"
    if path.exists() and not refresh:
        with open(path) as fh:
            return [json.loads(line) for line in fh if line.strip()], True
    rows = A.datasets_query(ncbi)
    NCBI_CACHE.mkdir(parents=True, exist_ok=True)
    with atomic_write(path) as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    return rows, False


def routes_for(rec: dict, internal, route_of: dict) -> str:
    """The membership routes for one species.

    Its OWN routes if it is a member in its own right; otherwise the routes of the clade it was
    found under. Routes 2 and 5 are expressed as GENERA, so almost no species is a member by
    those routes directly -- without this the whole congener and host_associate population
    lands in the table as `via_clade`, which says a taxon is present without saying why. That
    is the db_v3 failure this column exists to end, and `05_detect` depends on the answer: a
    taxon present only as a `congener` is there to absorb reads, not to be reported as a find.
    """
    own = route_of.get(internal) if internal else None
    if own:
        return own
    inherited = rec.get("unit_routes") or set()
    return ";".join(sorted(inherited)) if inherited else "via_clade"


def survey(queries: dict, route_of: dict, tax, workers: int, refresh: bool = False) -> tuple:
    """Query every unit, fold to species, and keep every assembly RECORD.

    Returns (per-species rows, per-assembly rows, stats). The per-assembly table is the
    reason this is the only step that talks to NCBI: `select_refs.py` then filters a table on
    disk instead of re-querying, so the selection is reproducible from a file and a catalogue
    change between the two steps cannot make them disagree.
    """
    per_species = defaultdict(lambda: {
        "assemblies": [], "found_under": set(), "unit_routes": set()})
    stats = defaultdict(int)

    units = [(tax.ncbi_of(tid), tax.label(tid), tid) for tid in queries]
    units = [(nc, nm, tid) for nc, nm, tid in units if nc]
    stats["units_without_ncbi_taxid"] = len(queries) - len(units)

    def work(unit):
        ncbi, name, tid = unit
        return unit, cached_query(ncbi, refresh=refresh)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(work, u) for u in units]
        for done, fut in enumerate(as_completed(futs), 1):
            (ncbi, name, tid), (rows, was_cached) = fut.result()
            stats["units_from_cache" if was_cached else "units_queried"] += 1
            if not rows:
                stats["units_returning_nothing"] += 1
            # The routes that made THIS query unit a member. A species found under it inherits
            # them when it is not itself a member, which is the whole point of the column.
            unit_routes = route_of.get(tid, "")
            for a in rows:
                organism = a.get("organism") or {}
                raw = organism.get("tax_id")
                if not raw:
                    stats["assemblies_without_taxid"] += 1
                    continue
                internal = tax.by_ncbi.get(raw)
                if internal is None:
                    # Newer than holobase's taxonomy snapshot. Counted under its own taxid so
                    # it is never silently lost, and reported so the gap is visible.
                    stats["assemblies_taxid_not_in_holobase"] += 1
                    sid_ncbi, sid_name = raw, organism.get("organism_name", str(raw))
                else:
                    sp = tax.species_of(internal)
                    if sp:
                        sid_ncbi, sid_name = tax.ncbi_of(sp) or raw, tax.label(sp)
                    else:
                        stats["assemblies_above_species_rank"] += 1
                        sid_ncbi, sid_name = raw, organism.get("organism_name", str(raw))
                rec = per_species[(sid_ncbi, sid_name)]
                rec["assemblies"].append(a)
                rec["found_under"].add(name)
                rec["unit_routes"].update(r for r in unit_routes.split(";") if r)
            if done % 25 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} queries … "
                      f"({len(per_species):,} species, "
                      f"{sum(len(v['assemblies']) for v in per_species.values()):,} assemblies)",
                      end="\r", flush=True)
    print()

    rows = []
    for (ncbi, name), rec in per_species.items():
        # One accession can be returned by two overlapping queries; dedupe before counting.
        seen, asms = set(), []
        for a in rec["assemblies"]:
            acc = a.get("accession", "")
            if acc and acc not in seen:
                seen.add(acc)
                asms.append(a)
        annotated = [a for a in asms if A.has_annotation(a)]
        scaffold = [a for a in asms if A.scaffold_plus(a)]
        sp_ann = [a for a in scaffold if A.has_annotation(a)]
        dates = sorted(d for a in asms if (d := A.release_date(a)))
        countries = {c for a in asms if (c := A.get_country(a))}
        internal = tax.by_ncbi.get(ncbi)
        rows.append({
            "ncbi_taxid": ncbi,
            "name": name,
            "rank": tax.rank.get(internal, "") if internal else "not_in_holobase",
            "routes": routes_for(rec, internal, route_of),
            "n_assemblies": len(asms),
            "n_annotated": len(annotated),
            "n_scaffold_plus": len(scaffold),
            "n_scaffold_plus_annotated": len(sp_ann),
            "best_level": max((A.assembly_level(a) for a in asms),
                              key=lambda l: A.LEVEL_RANK.get(l, 0), default=""),
            "first_release": dates[0] if dates else "",
            "last_release": dates[-1] if dates else "",
            "n_countries": len(countries),
            "found_under": "; ".join(sorted(rec["found_under"])[:3]),
        })
    rows.sort(key=lambda r: (-r["n_annotated"], -r["n_assemblies"], r["name"]))

    # Every assembly record, so the selection never needs NCBI again.
    asm_rows = []
    for (ncbi, name), rec in per_species.items():
        internal = tax.by_ncbi.get(ncbi)
        route = routes_for(rec, internal, route_of)
        seen = set()
        for a in rec["assemblies"]:
            acc = a.get("accession", "")
            if not acc or acc in seen:
                continue
            seen.add(acc)
            row = A.to_row(a, ncbi, name, "fungal", route, 0, "surveyed")
            row["routes"] = route
            asm_rows.append(row)
    asm_rows.sort(key=lambda r: (r["organism_name"].lower(), r["accession"]))
    return rows, asm_rows, stats


def report(rows: list, stats: dict) -> None:
    n_asm = sum(r["n_assemblies"] for r in rows)
    n_ann = sum(r["n_annotated"] for r in rows)
    n_sp_ann = sum(r["n_scaffold_plus_annotated"] for r in rows)
    with_ann = [r for r in rows if r["n_annotated"]]
    no_ann = [r for r in rows if r["n_assemblies"] and not r["n_annotated"]]

    print(f"\n── What exists, measured ────────────────────────────────────────")
    print(f"  species with >=1 assembly        {len(rows):>8,}")
    print(f"  species with >=1 GENE MODEL      {len(with_ann):>8,}   <- these can enter db_v4")
    print(f"  species with assemblies, NO model{len(no_ann):>8,}   <- the miniprot worklist")
    print()
    print(f"  assemblies total                 {n_asm:>8,}")
    print(f"  ... annotated                    {n_ann:>8,}   ({n_ann/n_asm:.0%})")
    print(f"  ... annotated AND scaffold-plus  {n_sp_ann:>8,}   ({n_sp_ann/n_asm:.0%})")
    print(f"\n  db_v3 for scale: 2,063 assemblies / 877 fungal species / 4.80 B minimizers")

    print(f"\n  depth distribution over species that can enter:")
    buckets = defaultdict(int)
    for r in with_ann:
        n = r["n_annotated"]
        key = "1" if n == 1 else "2-3" if n < 4 else "4-9" if n < 10 else "10-49" if n < 50 else "50+"
        buckets[key] += 1
    for key in ("1", "2-3", "4-9", "10-49", "50+"):
        if buckets[key]:
            print(f"    {key:>6} annotated assemblies: {buckets[key]:>6,} species")
    print(f"  the median is the reason a floor was never the binding constraint on breadth")

    print(f"\n  deepest species (annotated assemblies):")
    for r in with_ann[:10]:
        print(f"    {r['name'][:44]:<44} {r['n_annotated']:>5,} of {r['n_assemblies']:>5,}")

    for k in ("units_queried", "units_from_cache",
              "units_returning_nothing", "units_without_ncbi_taxid",
              "assemblies_taxid_not_in_holobase", "assemblies_above_species_rank",
              "assemblies_without_taxid"):
        if stats.get(k):
            print(f"\n  NOTE {k}: {stats[k]:,}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--congener-rank", default="genus", choices=("genus", "family"))
    ap.add_argument("--workers", type=int, default=8,
                    help="NCBI query concurrency. 8, not 16: 16 measured worse")
    ap.add_argument("--limit", type=int, default=0,
                    help="smoke test over the first N query units; writes nothing")
    ap.add_argument("--refresh", action="store_true",
                    help="re-query NCBI even where the bronze cache has an answer. The "
                         "catalogue moves, so use this to resurvey deliberately")
    args = ap.parse_args()

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "survey.log")
    link_latest(LOGS, log_dir / "survey.log")
    sys.stdout = log
    try:
        print("Building the routes list …\n")
        _, queries, _, tax, route_of = R.build(args.congener_rank)
        if args.limit:
            queries = dict(list(queries.items())[:args.limit])
            print(f"\n*** --limit {args.limit}: SMOKE TEST. Writes nothing. ***")
        # _assembly's to_row needs the phylum map (BUSCO lineage) and the species names.
        # Filled from the Taxonomy already in memory rather than by two more full passes over
        # the 3 M-row table, which is why _taxonomy owns the single read.
        print("\n  filling the BUSCO phylum map from the taxonomy already loaded …")
        A.PHYLUM_OF = tax.phylum_map()
        A.SPECIES_NAME = {tax.ncbi_of(t): tax.label(t)
                          for t, rk in tax.rank.items() if rk == "species" and t in tax.ncbi}

        print(f"\nQuerying NCBI for {len(queries):,} clades ({args.workers} workers) …")
        rows, asm_rows, stats = survey(queries, route_of, tax, args.workers, args.refresh)
        report(rows, stats)

        if args.limit:
            print("\n(--limit given: nothing written)")
            return
        SURVEY_TSV.parent.mkdir(parents=True, exist_ok=True)
        with open(SURVEY_TSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=SURVEY_COLS, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        worklist = [r for r in rows if r["n_assemblies"] and not r["n_annotated"]]
        with open(MINIPROT_TSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=SURVEY_COLS, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(worklist)
        with open(ASSEMBLIES_TSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=A.CANDIDATE_COLS + ["routes"], delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(asm_rows)
        print(f"\n  -> {SURVEY_TSV}        ({len(rows):,} species)")
        print(f"  -> {ASSEMBLIES_TSV}        ({len(asm_rows):,} assembly records)")
        print(f"  -> {MINIPROT_TSV}  ({len(worklist):,} species with assemblies but no model)")
        print(f"\nNext: python 02_classify/01_refs/select_refs.py --select-only   "
              f"(filters the table above; no further NCBI queries)")
    finally:
        log.close()


if __name__ == "__main__":
    main()
