#!/usr/bin/env python3
"""compare_refs.py — diff db_v4's fresh selection against a previous one, and say WHY.

    python 02_classify/01_refs/compare_refs.py                      # v4 vs db_v3
    python 02_classify/01_refs/compare_refs.py --probe-ncbi         # resolve the unknowns
    python 02_classify/01_refs/compare_refs.py --against <table>    # any prior selection

`select_refs.py` deliberately has no `--retain`: it selects from the rules alone, so nothing is
inherited through a pointer at a generated table. **This script is what makes that safe, and it
is not optional.** The floor is a cap for every species in a fresh run, so a fresh selection
TRIMS depth that db_v3 held, and the depth policy measured on 2026-09-24 says that costs more
than it buys: capping a deep species was worth 0.42 pp of congener specificity against a median
16% of that species' observed pangenome, up to 62%, and pangenomes are open (Heaps gamma ~0.16,
no asymptote). Without this diff, dropping `--retain` is a silent regression. With it, retention
becomes a reviewed decision with a reason per row.

## Why a plain set difference is not enough

An accession in the old table and absent from the new one has at least four unrelated causes,
and they call for opposite responses. NCBI's catalogue moves between a selection and the next
one, so some of the difference is not our decision at all:

    floor_cap             its species IS selected in v4, but this accession fell outside
                          min(floor, available). Our rules trimmed depth. RETAIN CANDIDATE.
    excluded              the taxid is in seeds/exclusions.tsv. A decision already taken,
                          with a reason on the row. Not a candidate.
    route_no_cds          the taxon is in USDA but routed miniprot (no gene model) or none.
                          Cannot enter a CDS-only database yet. Waits on utilities/annotate.
    taxon_out_of_scope    the species is not in seeds/ at all, by NAME or by species node.
                          db_v3 seeded on PHI-base and swept by ORDER, so it holds congeners
                          with no plant-association record. This is the bucket that carries
                          the open question: those congeners are out of SCOPE but they were
                          in the database as LCA COMPETITORS, which is a different job.
    oomycete_out_of_scope removed by the fungi-only decision of 2026-10-13, not by absence
                          from USDA. A decision, with its reason in exclusions.tsv.
    ncbi_absent           NCBI no longer returns the accession for that taxon: suppressed,
                          replaced, or re-assigned. Retaining it would retain a dead record.
    annotation_withdrawn  still returned, but no longer carries a gene model.

The last two cannot be told apart from the two tables alone, so they are reported as
`unresolved_needs_probe` unless `--probe-ncbi` is passed, which re-queries the affected taxa.
Guessing would be worse than labelling: a retain list containing suppressed accessions fails at
download time, and one that drops live assemblies loses depth silently.

## Output

    data/silver/refs_diff.tsv          every accession that moved, with its bucket
    data/silver/retain_candidates.tsv  the floor_cap rows only, the review worklist

`retain_candidates.tsv` is SILVER, so it is generated and disposable. What survives review gets
promoted by hand into `seeds/retain.tsv`, which is AUTHORED and tracked, carrying accession,
species and a reason. That is the whole point: db_v3 inherited through
`--retain ref_candidates_legacy.tsv`, a generated table acting as an authority, which had
already pulled 785 assemblies including 156 *P. nodorum* into a side run. A seeds file with a
reason per row cannot go stale without somebody seeing it.
"""

import argparse
import csv
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _layout import SILVER, REF_CANDIDATES, LOGS, FUNGI_ANNOTATED, SEEDS   # noqa: E402
from _seeds import read_seed_tsv                                           # noqa: E402
from _util import _Tee, make_log_dir, link_latest                          # noqa: E402
import _assembly as A                                                      # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _paths import REF_CANDIDATES as DBV3_CANDIDATES                       # noqa: E402

DIFF_TSV = SILVER / "refs_diff.tsv"
RETAIN_TSV = SILVER / "retain_candidates.tsv"
RETAIN_SEED = SEEDS / "retain.tsv"

DIFF_COLS = ["direction", "bucket", "accession", "taxid", "organism_name",
             "assembly_level", "has_annotation", "old_depth", "new_depth",
             "old_reason", "note"]


def load_table(path: Path) -> list:
    """A candidate table as a list of dicts. These are GENERATED, so no preamble handling."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path} has no rows")
    if "accession" not in rows[0]:
        raise ValueError(f"{path}: no `accession` column. Found {list(rows[0])}")
    return rows


def seed_routes() -> dict:
    """{taxid: route} over every USDA taxon, addable or not, so a drop can name its route."""
    out = {}
    for row in read_seed_tsv(FUNGI_ANNOTATED,
                             require=("ncbi_taxid", "route")):
        tid = (row.get("ncbi_taxid") or "").strip()
        if tid.isdigit():
            out[int(tid)] = (row.get("route") or "").strip()
    return out


def _taxid(row) -> int:
    tid = str(row.get("taxid", "")).strip()
    return int(tid) if tid.isdigit() else 0


def bucket_drop(row, new_taxa: set, routes: dict, route_species: dict,
                excluded: dict, to_species: dict) -> str:
    """Why is this old accession not in the new selection? Cheap causes only.

    Ordered most-specific first. `floor_cap` is the only bucket that means "our depth rule
    trimmed it", so it must not absorb drops that have a scope explanation: a species excluded
    by a seeds decision is not a retain candidate even though its accession also fell outside
    the floor.

    The taxid is normalised to its species node on BOTH sides before any comparison, and that
    is not a detail. The first version of this function did not, and it labelled 972 assemblies
    `taxon_out_of_scope` — which reads as "a scope decision already taken" — with the list
    headed by *Pyricularia oryzae* (45), *Fusarium oxysporum* (57) and *Verticillium dahliae*
    (11). Those are major cereal pathogens and no such decision was ever taken about them. The
    cause was RANK, not scope: holobase records *F. oxysporum* only as formae speciales and
    never at the species node 5507 that db_v3 holds, so the species looked absent from seeds/.
    A bucket that confident about a number that wrong is worse than no bucket at all.

    Oomycetes are split out for the same reason: "removed by the fungi-only decision" and
    "absent from the USDA list" are different statements and only one of them is a decision.
    """
    tid = _taxid(row)
    sid = to_species.get(tid, tid)
    if tid in excluded or sid in excluded:
        return "excluded"
    route = routes.get(tid) or route_species.get(sid)
    if route in ("miniprot", "none"):
        return "route_no_cds"
    if tid in new_taxa or sid in new_taxa:
        # The species is in the new selection, so seeds authorised it and NCBI served it.
        # This accession simply fell outside min(floor, available).
        return "floor_cap"
    if route is None:
        # Fungi only, Leon 2026-10-13: db_v3's 29 oomycete species come out by decision, and
        # that is not the same statement as "absent from the USDA list".
        if (row.get("kingdom") or "").strip() == "oomycete":
            return "oomycete_out_of_scope"
        return "taxon_out_of_scope"
    # route is present/add but the species produced nothing: either NCBI dropped the taxon's
    # assemblies or select_refs reported it under `no annotated assembly`. Needs a probe.
    return "unresolved_needs_probe"


def probe(taxa: list, workers: int) -> dict:
    """{taxid: {accession: annotated?}} straight from NCBI, for the unresolved drops only.

    Bounded by the number of distinct species with unresolved drops, not by accessions, because
    one query returns the whole taxon. `datasets_query` warns on a non-zero exit rather than
    returning [] silently, which matters here: a homonym would otherwise look like deletion.
    """
    out = {}

    def work(tid):
        return tid, {a.get("accession", ""): A.has_annotation(a)
                     for a in A.datasets_query(tid)}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(work, t) for t in taxa]
        for done, fut in enumerate(as_completed(futs), 1):
            tid, got = fut.result()
            out[tid] = got
            if done % 10 == 0 or done == len(futs):
                print(f"  probed {done}/{len(futs)} taxa …", end="\r", flush=True)
    print()
    return out


def run_compare(new_path: Path, old_path: Path, do_probe: bool, workers: int) -> None:
    new_rows, old_rows = load_table(new_path), load_table(old_path)
    routes = seed_routes()

    excluded = {}
    for row in read_seed_tsv(SEEDS / "exclusions.tsv", require=("ncbi_taxid", "kind")):
        tid = (row.get("ncbi_taxid") or "").strip()
        if tid.isdigit():
            excluded[int(tid)] = (row.get("kind") or "").strip()

    print("Reading holobase to normalise taxids on both sides …")
    to_species = A.init_taxonomy()
    # A route recorded against a forma specialis has to be findable from the species node the
    # old table used, and vice versa. Without this the diff reports a rank mismatch as a scope
    # decision; see bucket_drop.
    route_species = {}
    for tid, route in routes.items():
        route_species.setdefault(to_species.get(tid, tid), route)

    new_by_acc = {r["accession"]: r for r in new_rows}
    old_by_acc = {r["accession"]: r for r in old_rows}
    new_depth, old_depth = defaultdict(int), defaultdict(int)
    for r in new_rows:
        new_depth[_taxid(r)] += 1
    for r in old_rows:
        old_depth[_taxid(r)] += 1
    # Depth is compared per SPECIES, so an old table holding a species at its species node and
    # a new one reaching it through a forma specialis are not read as two different taxa.
    new_taxa = set(new_depth) | {to_species.get(t, t) for t in new_depth}

    print(f"  new  {new_path}")
    print(f"       {len(new_rows):,} assemblies / {len(new_taxa):,} species")
    print(f"  old  {old_path}")
    print(f"       {len(old_rows):,} assemblies / {len(old_depth):,} species")
    print(f"  shared accessions: {len(set(new_by_acc) & set(old_by_acc)):,}")

    diff = []
    for acc, r in sorted(new_by_acc.items()):
        if acc in old_by_acc:
            continue
        tid = _taxid(r)
        diff.append({"direction": "added", "bucket": r.get("source", ""),
                     "accession": acc, "taxid": tid,
                     "organism_name": r.get("organism_name", ""),
                     "assembly_level": r.get("assembly_level", ""),
                     "has_annotation": r.get("has_annotation", ""),
                     "old_depth": old_depth.get(tid, 0), "new_depth": new_depth.get(tid, 0),
                     "old_reason": "", "note": r.get("selection_reason", "")})

    dropped = []
    for acc, r in sorted(old_by_acc.items()):
        if acc in new_by_acc:
            continue
        tid = _taxid(r)
        b = bucket_drop(r, new_taxa, routes, route_species, excluded, to_species)
        rec = {"direction": "dropped", "bucket": b, "accession": acc, "taxid": tid,
               "organism_name": r.get("organism_name", ""),
               "assembly_level": r.get("assembly_level", ""),
               "has_annotation": r.get("has_annotation", ""),
               "old_depth": old_depth.get(tid, 0), "new_depth": new_depth.get(tid, 0),
               "old_reason": r.get("selection_reason", ""),
               "note": excluded.get(tid, "") or routes.get(tid, "")}
        dropped.append(rec)

    unresolved = sorted({r["taxid"] for r in dropped
                         if r["bucket"] == "unresolved_needs_probe" and r["taxid"]})
    if unresolved and do_probe:
        print(f"\nProbing NCBI for {len(unresolved)} taxa with unresolved drops …")
        live = probe(unresolved, workers)
        for r in dropped:
            if r["bucket"] != "unresolved_needs_probe":
                continue
            held = live.get(r["taxid"], {})
            if r["accession"] not in held:
                r["bucket"], r["note"] = "ncbi_absent", "not returned for this taxon"
            elif not held[r["accession"]]:
                r["bucket"], r["note"] = "annotation_withdrawn", "returned, no gene model"
            else:
                r["bucket"] = "floor_cap"
                r["note"] = "live and annotated; outside min(floor, available)"
    elif unresolved:
        print(f"\n  {len(unresolved)} taxa have drops this diff CANNOT explain from the "
              f"tables alone.\n  Re-run with --probe-ncbi to separate a suppressed accession "
              f"from a trimmed one.\n  Until then they are NOT counted as retain candidates.")

    diff.extend(dropped)
    DIFF_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(DIFF_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=DIFF_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(diff)

    census = defaultdict(int)
    for r in diff:
        census[(r["direction"], r["bucket"])] += 1
    print(f"\n── Diff ─────────────────────────────────────────────────────────")
    for (direction, b), n in sorted(census.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {direction:8} {b:24} {n:6,}")
    print(f"  -> {DIFF_TSV}")

    retain = [r for r in diff if r["bucket"] == "floor_cap"]
    by_species = defaultdict(list)
    for r in retain:
        by_species[(r["taxid"], r["organism_name"])].append(r)
    with open(RETAIN_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=DIFF_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(retain)

    print(f"\n── Retain candidates: depth our rules trimmed ───────────────────")
    print(f"  {len(retain):,} assemblies across {len(by_species):,} species")
    if by_species:
        worst = sorted(by_species.items(), key=lambda kv: -len(kv[1]))[:15]
        print(f"  {'species':<42} {'old':>5} {'new':>5} {'trimmed':>8}")
        for (tid, name), rows in worst:
            print(f"  {name[:42]:<42} {old_depth.get(tid, 0):>5} "
                  f"{new_depth.get(tid, 0):>5} {len(rows):>8}")
    print(f"  -> {RETAIN_TSV}")
    print(f"\nReview it, then promote what is worth keeping into {RETAIN_SEED}")
    print("(authored and tracked: accession, species, reason per row). The CDS for these is "
          "already on Acacia under db/search, so retaining costs no download.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--new", default=None,
                    help=f"the fresh selection (default {REF_CANDIDATES})")
    ap.add_argument("--against", default=None,
                    help=f"the previous selection (default {DBV3_CANDIDATES})")
    ap.add_argument("--probe-ncbi", action="store_true",
                    help="re-query NCBI for taxa whose drops the tables cannot explain, to "
                         "separate a suppressed accession from one our floor trimmed")
    ap.add_argument("--workers", type=int, default=8,
                    help="NCBI query concurrency. 8, not 16")
    args = ap.parse_args()

    new_path = Path(args.new) if args.new else REF_CANDIDATES
    old_path = Path(args.against) if args.against else DBV3_CANDIDATES
    if not new_path.exists():
        ap.error(f"{new_path} not found — run select_refs.py --select-only first")
    if not old_path.exists():
        ap.error(f"{old_path} not found")

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "compare_refs.log")
    link_latest(LOGS, log_dir / "compare_refs.log")
    sys.stdout = log
    try:
        run_compare(new_path, old_path, args.probe_ncbi, args.workers)
    finally:
        log.close()


if __name__ == "__main__":
    main()
