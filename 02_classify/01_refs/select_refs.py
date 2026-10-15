#!/usr/bin/env python3
"""select_refs.py — choose db_v4's reference assemblies from seeds/, and nothing else.

Step 1 of 02_classify. Replaces `03_kraken/01_search/kraken_search.py`, whose SELECTION
MECHANICS are carried across unchanged in `_assembly.py` and whose POLICY is replaced here.

    python 02_classify/01_refs/select_refs.py --select-only --limit 20   # smoke test first
    python 02_classify/01_refs/select_refs.py --select-only              # the real selection
    python 02_classify/01_refs/select_refs.py --download                 # select, then fetch CDS

## What selects what

`seeds/` is the membership authority, so this script asks it three questions and nothing else:

    fungi_annotated.tsv  route=present (432) + route=add (179) = 611 taxa selected DIRECTLY,
                         one NCBI query each. route=miniprot (849) have a scaffold-plus
                         assembly but NO gene model, so they cannot enter a CDS-only database
                         until miniprot has run; route=none (26) have nothing. Both are
                         reported and skipped, never silently dropped.
    residents.tsv        30 genera swept, species not already covered by the above. Column 1
                         may be a TAXID rather than a name, and must be for a homonym.
    exclusions.tsv       a blocklist by taxid, applied to both. Nothing leaves the build
                         without a row there saying why.

## The four policy changes from db_v3, each with its reason

**1. USDA replaces PHI-base, and breadth stops coming from a sweep.** db_v3 seeded on
PHI-base's ~340 fungi and bought breadth by sweeping every ORDER containing a seed, because
the seed list was too narrow to be the breadth. USDA's `recorded_on` relation reaches 1,486
plant-associated fungi with a scaffold-plus assembly directly, with a record behind each one,
so the taxa are selected on evidence instead of on taxonomic proximity to a pathogen. The
order sweep is therefore gone: it would now add non-plant-associated congeners, which the
scope rule excludes. The residents sweep stays, because `residents.tsv` names genera and the
species under them have to be found.

**2. No STAT, including the gate nobody had counted.** `load_detected_genera()` read
`01_stat/03_filter/data/runs.tsv` and derived 38 STAT-detected genera; `select_seed` then gave
those genera UNLIMITED depth and everything else `min(floor, available)`. That is a third STAT
gate, in the selector, on top of the two the module README names, and it decides database
COMPOSITION. Worse, it fails silently in the direction of disaster: the loader catches
FileNotFoundError and treats ALL genera as detected, so once 01_stat moves to legacy the cap
lifts everywhere and *Aureobasidium* arrives at 135 assemblies instead of 11.

It is not re-pointed at `declared_pathogens`, which would be the obvious literature-first
analogue, because the asymmetry runs the wrong way for this study. Deep sampling of one
species takes k-mer space FROM its shallow congeners: a k-mer in both pangenomes but absent
from the shallow species' assemblies is stored as the deep one, so the shallow species' reads
score for it. Giving declared pathogens unlimited depth would therefore bias the database
toward confirming declarations and against their undeclared congeners, which are the cryptic
findings the chapter exists to detect. So the floor applies UNIFORMLY: every species gets
`min(floor, available)`, with no genus privileged by what a previous classifier happened to
see.

**3. No `--retain`.** db_v3 inherited its predecessor's selection through a pointer at a
generated table, which defaulted to the db_v2-era `ref_candidates_legacy.tsv` and had already
contaminated one side run with 785 assemblies including 156 *P. nodorum*. A generated table was
acting as an authority, which is the exact confusion the seeds/data split exists to end. This
script selects from the rules alone; `compare_refs.py` then diffs the result against any prior
table, and anything worth keeping is promoted into `seeds/retain.tsv` with a reason per row.

Dropping it has a real cost and it is the reason `compare_refs.py` is not optional: the floor
is a cap for every species here, so a fresh run trims db_v3's deep species. The measured
trade is 0.42 pp of congener specificity gained against a median 16% of a species' observed
pangenome lost, up to 62%, and pangenomes are open (Heaps gamma ~0.16). Expect to retain.

**4. The contig-level annotation gate is relaxed** (Leon, 2026-10-15). db_v3 required
`scaffold_plus AND has_annotation`, and for six genera those never co-occur because the only
annotated assembly is contig-level: *Boeremia*, *Dioszegia*, *Papiliotrema*, *Clohesyomyces*,
*Leucosporidium*, *Phaeosphaeria*. That excluded *Papiliotrema*, rice's single most abundant
genus at 14.5%, and *Dioszegia*, a wheat core yeast. Relaxed as a FALLBACK rather than a
removal: scaffold-plus annotated assemblies are still preferred, and a contig-level one is
taken only where its species has no better. So nothing already in db_v3 changes, and the 9
recovered assemblies across 6 species arrive where there was previously nothing. `--strict-gate`
restores db_v3's behaviour for comparison.

A resident exists to be a COMPETITOR, somewhere for a read to park, and that needs real k-mers
rather than a good assembly. Reference incompleteness largely cancels in the coverage ratio
because a fragmented assembly shrinks numerator and denominator together.

## What this script does NOT do

  * **Hosts.** `seeds/hosts.tsv` names 6 plant species and they need genomic, not CDS,
    assemblies; `03_kraken/01_search/kraken_hosts.py` still owns that and has not been moved.
  * **miniprot.** The 849 route=miniprot taxa are where the remaining breadth lives and they
    need gene prediction first. `utilities/annotate` owns it.
  * **BUSCO.** Scored in `utilities/`, reported, and NOT used to drop anything: at the 50%
    threshold 2 of 1,018 assemblies fail, and a per-taxid fallback means 0 of 309 taxa were
    ever lost to it. The real attrition is `no_cds`, 20x more common.
"""

import argparse
import csv
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _layout import (SEEDS, SILVER, CDS_PATHOGEN, REF_CANDIDATES, LOGS,   # noqa: E402
                     FUNGI_ANNOTATED, RESIDENTS, EXCLUSIONS)
from _seeds import read_seed_tsv, read_seed_rows                           # noqa: E402
from _util import _Tee, make_log_dir, link_latest                          # noqa: E402
import _assembly as A                                                      # noqa: E402

DEFAULT_FLOOR = 4          # measured 2026-09-24: a floor, never a cap. See the docstring
DEFAULT_WORKERS = 8        # --workers 8, not 16; 16 is worse for NCBI (measured)
ADDABLE_ROUTES = ("present", "add")
KINGDOM = "fungal"         # fungi only, Leon 2026-10-13. Drives the BUSCO lineage fallback


# ── seeds ─────────────────────────────────────────────────────────────────────

def load_excluded() -> dict:
    """{taxid: kind} from exclusions.tsv. Nothing may leave the build without a row here.

    Read through `_seeds` rather than a bare DictReader: exclusions.tsv opens with a 9-line
    `#` preamble, which a DictReader takes as the header, and the loader then blocks nothing
    while reporting success. See `_seeds.py` for what that would have cost.
    """
    out = {}
    for row in read_seed_tsv(EXCLUSIONS, require=("ncbi_taxid", "kind")):
        tid = (row.get("ncbi_taxid") or "").strip()
        if tid.isdigit():
            out[int(tid)] = (row.get("kind") or "").strip() or "unspecified"
    if not out:
        raise ValueError(f"{EXCLUSIONS} parsed but yielded no taxids to block")
    return out


def load_usda_taxa(excluded: dict) -> tuple:
    """Addable USDA taxa as [(taxid, name)], plus the route census.

    route=present and route=add are selected. route=miniprot and route=none are counted and
    reported: 80% of what USDA offers has no gene model, and a CDS-only database cannot hold
    it, so the number is the headline constraint on breadth rather than an aside.
    """
    taxa, census, blocked = [], defaultdict(int), []
    for row in read_seed_tsv(FUNGI_ANNOTATED,
                             require=("ncbi_taxid", "scientific_name", "route")):
        route = (row.get("route") or "").strip()
        census[route] += 1
        if route not in ADDABLE_ROUTES:
            continue
        tid = (row.get("ncbi_taxid") or "").strip()
        if not tid.isdigit():
            continue
        tid = int(tid)
        if tid in excluded:
            blocked.append((tid, row.get("scientific_name", ""), excluded[tid]))
            continue
        taxa.append((tid, (row.get("scientific_name") or str(tid)).strip()))
    return taxa, dict(census), blocked


def load_residents() -> list:
    """[(clade, note)] from residents.tsv, skipping anything not marked `add`.

    Column 1 may be a TAXID rather than a name and must be for a taxonomic homonym: the
    `datasets` CLI refuses an ambiguous name, which used to select zero assemblies silently.
    *Cryptococcus* is a fungus (5206) and a scale insect (79213); *Rhizophagus* a fungus
    (1129544) and a beetle (295919). Between them that cost 42 assemblies and 17 species.
    """
    out = []
    for parts in read_seed_rows(RESIDENTS, min_cols=3):
        if len(parts) < 3 or parts[2].strip() != "add":
            continue
        out.append((parts[0].strip(), parts[3].strip() if len(parts) > 3 else ""))
    if not out:
        raise ValueError(f"{RESIDENTS} parsed but yielded no clades marked `add`")
    return out


# ── selection ─────────────────────────────────────────────────────────────────

def collapse_to_species(taxa: list, to_species: dict) -> tuple:
    """Normalise every seed taxid to its SPECIES node, merging those that collapse together.

    USDA records host association at the rank the record was made at, so `fungi_usda.tsv`
    carries 37 forma specialis, 2 varietas and 2 forma rows alongside its 1,445 species; 16 of
    those fall in the addable set. *Fusarium oxysporum* is the case that matters: holobase
    lists it ONLY as formae speciales (f. sp. melonis, lycopersici, cubense, 12 of them in the
    addable set) and never at the species node 5507, which db_v3 holds.

    Selecting per seed taxid would then repeat the trap `select_rank_sweep` documents, in the
    one place nothing was guarding it: the floor would count per forma specialis, so one
    species arrives as 12 depth-4 taxa, each passing the floor while the species itself is
    sampled 12x over and reported as twelve things. *F. oxysporum* is also one of the
    worst-measured misassignment partners, so getting its depth wrong is not cosmetic.

    Querying the SPECIES node instead returns the whole clade, including every forma
    specialis' assemblies, and the floor then applies once to the species. Detections are
    reported at species rank and `declared_species_set()` already promotes a declared forma
    specialis to its species node, so both sides agree.
    """
    merged, members = {}, defaultdict(list)
    for tid, name in taxa:
        sid = to_species.get(tid, tid)
        members[sid].append(tid)
        # Prefer a name recorded at species rank; a forma specialis name would otherwise
        # label the whole species after whichever member happened to be first.
        merged[sid] = A.SPECIES_NAME.get(sid) or merged.get(sid) or name
    out = [(sid, merged[sid]) for sid in merged]
    collapsed = [(sid, merged[sid], m) for sid, m in members.items()
                 if len(m) > 1 or m[0] != sid]
    return out, sorted(collapsed, key=lambda c: -len(c[2]))



def annotated_pool(assemblies: list, strict: bool) -> list:
    """Annotated assemblies for one species, best tier available.

    Annotation is non-negotiable: without a gene model there is no CDS and the assembly
    cannot enter the database at all. Assembly LEVEL is a preference, applied as a fallback
    so a contig-level annotated assembly is taken only where its species has nothing better.
    `strict=True` restores db_v3's `scaffold_plus AND has_annotation` gate.
    """
    annotated = [a for a in assemblies if A.has_annotation(a)]
    scaffolded = [a for a in annotated if A.scaffold_plus(a)]
    if strict:
        return scaffolded
    return scaffolded or annotated


def select_taxon(taxid: int, name: str, assemblies: list, floor: int,
                 strict: bool) -> list:
    """Select for one USDA taxon: min(floor, available), diversity-ordered.

    The floor is a target and not a filter, so a species with one assembly in existence is
    included at one. It still earns its place by giving Kraken2's LCA a competitor, which is
    what stops a read being handed to whichever species happens to be in the build. Its own
    species-level counts are not reliable and should be aggregated to genus when reported.
    """
    pool = annotated_pool(assemblies, strict)
    if not pool:
        return []
    pool.sort(key=A.quality_key, reverse=True)
    pairs = A.greedy_ordered(pool)
    if floor > 0:
        pairs = pairs[:floor]
    level = "" if any(A.scaffold_plus(a) for a, _ in pairs) else "; contig_fallback"
    return [A.to_row(a, taxid, name, KINGDOM, "usda", rank + 1, reason + level)
            for rank, (a, reason) in enumerate(pairs)]


def select_sweep(clade: str, assemblies: list, covered: set, excluded: dict,
                 to_species: dict, floor: int, strict: bool) -> list:
    """Sweep one resident clade: every species in it not already covered, min(floor, avail).

    Taxids are normalised to their SPECIES node first. NCBI labels many assemblies at strain,
    varietas or forma rank, and grouping on the assembly's own taxid breaks this two ways that
    both defeat the design: a covered species' own strains stop being recognised as covered
    (159 re-selected at order scope), and the floor counts per taxid so a species split across
    several strain nodes passes the floor while staying shallow (87 species, worst
    *Zymoseptoria tritici* across 6 nodes).
    """
    by_species = defaultdict(list)
    for a in assemblies:
        organism = a.get("organism") or {}
        tid = organism.get("tax_id")
        if not tid:
            continue
        sid = to_species.get(tid, tid)
        if sid in covered or tid in covered or sid in excluded or tid in excluded:
            continue
        if A.has_annotation(a):
            by_species[sid].append(a)

    rows = []
    for sid, pool in by_species.items():
        pool = annotated_pool(pool, strict)
        if not pool:
            continue
        pool.sort(key=A.quality_key, reverse=True)
        # Prefer a name recorded at species rank; strain names carry isolate suffixes that
        # would make the same species look like several.
        name = A.SPECIES_NAME.get(sid) or \
            (pool[0].get("organism") or {}).get("organism_name", str(sid))
        pairs = A.greedy_ordered(pool)
        if floor > 0:
            pairs = pairs[:floor]
        level = "" if any(A.scaffold_plus(a) for a, _ in pairs) else "; contig_fallback"
        for rank, (a, reason) in enumerate(pairs):
            rows.append(A.to_row(a, sid, name, KINGDOM, "resident_sweep",
                                 rank + 1, reason + level))
    return rows


def run_select(floor: int, workers: int, strict: bool, limit: int = 0) -> list:
    excluded = load_excluded()
    taxa, census, blocked = load_usda_taxa(excluded)
    residents = load_residents()

    print(f"Floor: {floor or 'none'} per species, applied UNIFORMLY (no STAT no-cap set)")
    print(f"Annotation gate: {'strict, scaffold-plus only (db_v3)' if strict else 'relaxed, contig-level as a fallback'}")
    print(f"Exclusions: {len(excluded)} taxids blocked")
    print("\nseeds/fungi_annotated.tsv routes:")
    for route in ("present", "add", "miniprot", "none"):
        n = census.get(route, 0)
        mark = "selected" if route in ADDABLE_ROUTES else "SKIPPED, no gene model"
        print(f"  {route:9} {n:5,}   {mark}")
    if blocked:
        print(f"  {len(blocked)} addable taxa blocked by exclusions.tsv:")
        for tid, name, reason in blocked[:10]:
            print(f"      {tid:>8}  {name:<40} {reason}")
    if limit:
        taxa, residents = taxa[:limit], residents[:2]
        print(f"\n*** --limit {limit}: SMOKE TEST over {len(taxa)} taxa and "
              f"{len(residents)} clades. Do not treat the output as a selection. ***")

    print("\nReading holobase (species map + phylum map) …")
    to_species = A.init_taxonomy()
    taxa, collapsed = collapse_to_species(taxa, to_species)
    covered = {tid for tid, _ in taxa}
    if to_species:
        covered |= {to_species[t] for t in covered if t in to_species}
        print(f"  normalised {len(to_species):,} taxids to species rank")
        if collapsed:
            print(f"  collapsed {len(collapsed)} sub-species seed taxa into their species "
                  f"node, so the floor applies per SPECIES:")
            for sid, name, members in collapsed[:8]:
                print(f"      {sid:>8}  {name:<36} <- {len(members)} seed taxa")
        print(f"  querying {len(taxa)} species nodes; sweep coverage {len(covered)} taxa")

    print(f"\n[1/2] Querying {len(taxa):,} USDA taxa ({workers} workers) …")
    rows, no_pool = [], []

    def work_taxon(taxid, name):
        got = select_taxon(taxid, name, A.datasets_query(taxid), floor, strict)
        return taxid, name, got

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(work_taxon, t, n) for t, n in taxa]
        for done, fut in enumerate(as_completed(futs), 1):
            taxid, name, got = fut.result()
            if got:
                rows.extend(got)
            else:
                no_pool.append((taxid, name))
            if done % 25 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} taxa … ({len(rows):,} assemblies)",
                      end="\r", flush=True)
    print()
    usda_count = len(rows)
    print(f"  USDA assemblies selected: {usda_count:,} across "
          f"{len({r['taxid'] for r in rows}):,} species")
    if no_pool:
        # Not an error and not silent: seeds said addable, NCBI returned nothing usable. The
        # catalogue moves between a --annotate run and a selection, so this is the drift.
        print(f"  {len(no_pool)} taxa marked addable returned no annotated assembly:")
        for tid, name in no_pool[:10]:
            print(f"      {tid:>8}  {name}")

    print(f"\n[2/2] Sweeping {len(residents)} resident clades …")

    def work_clade(clade):
        return select_sweep(clade, A.datasets_query(clade), covered, excluded,
                            to_species, floor, strict)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work_clade, c): c for c, _ in residents}
        for done, fut in enumerate(as_completed(futs), 1):
            rows.extend(fut.result())
            if done % 5 == 0 or done == len(futs):
                print(f"  {done}/{len(residents)} clades … "
                      f"({len(rows) - usda_count:,} swept)", end="\r", flush=True)
    print()
    print(f"  Resident assemblies selected: {len(rows) - usda_count:,}")

    # One accession can qualify under a USDA taxon and again under a resident sweep. Keep the
    # USDA row: it carries the taxid the seeds file authorised.
    seen, deduped, dropped = set(), [], 0
    for r in sorted(rows, key=lambda r: 0 if r["source"] == "usda" else 1):
        if r["accession"] in seen:
            dropped += 1
            continue
        seen.add(r["accession"])
        deduped.append(r)
    if dropped:
        print(f"  Dropped {dropped} duplicate accessions (same assembly under both sources)")

    deduped.sort(key=lambda r: (r["source"], r["organism_name"].lower(),
                                int(r["selection_rank"] or 0)))
    return deduped


def write_candidates(rows: list, out: Path = REF_CANDIDATES) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=A.CANDIDATE_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    by_source = defaultdict(int)
    contig = 0
    for r in rows:
        by_source[r["source"]] += 1
        if "contig_fallback" in (r["selection_reason"] or ""):
            contig += 1
    print(f"\n── Selection ────────────────────────────────────────────────────")
    print(f"  assemblies      {len(rows):,}")
    print(f"  species         {len({r['taxid'] for r in rows}):,}")
    for src, n in sorted(by_source.items()):
        print(f"    {src:16} {n:,}")
    print(f"  with CDS        {sum(1 for r in rows if r['fasta_type'] == 'cds'):,}")
    print(f"  contig fallback {contig:,}   (0 under --strict-gate)")
    print(f"  -> {out}")
    print("\nNext: compare against the previous selection before downloading anything:")
    print("  python 02_classify/01_refs/compare_refs.py "
          "--against 03_kraken/01_search/data/ref_candidates_v3.tsv")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--floor", type=int, default=DEFAULT_FLOOR,
                    help="assemblies per species; 0 takes every annotated assembly. "
                         "Applied uniformly, unlike db_v3")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                    help="NCBI query concurrency. 8, not 16: 16 measured worse")
    ap.add_argument("--strict-gate", action="store_true",
                    help="require scaffold-plus AND annotation, as db_v3 did. For comparison; "
                         "costs 9 assemblies across 6 species incl. Papiliotrema")
    ap.add_argument("--limit", type=int, default=0,
                    help="smoke test: select over the first N taxa and 2 clades only")
    ap.add_argument("--out", default=None, help=f"output TSV (default {REF_CANDIDATES})")
    ap.add_argument("--select-only", action="store_true",
                    help="select and write the table, download nothing")
    ap.add_argument("--download", action="store_true",
                    help="after selecting, download CDS for every candidate")
    args = ap.parse_args()

    if not (args.select_only or args.download):
        ap.error("pass --select-only or --download")

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "select_refs.log")
    link_latest(LOGS, log_dir / "select_refs.log")
    sys.stdout = log
    try:
        rows = run_select(args.floor, args.workers, args.strict_gate, args.limit)
        out = Path(args.out) if args.out else REF_CANDIDATES
        if args.limit:
            # A smoke run must never land on the path the build reads. The db_v3-era
            # equivalent was a symlink that a side run could repoint, which silently made a
            # partial table "the database"; here the output path simply moves instead.
            out = SILVER / "ref_candidates_smoke.tsv"
            print(f"\n--limit given: writing to {out}, NOT {REF_CANDIDATES}")
        write_candidates(rows, out)
        if args.download:
            A.run_download(rows, CDS_PATHOGEN, args.workers)
    finally:
        log.close()


if __name__ == "__main__":
    main()
