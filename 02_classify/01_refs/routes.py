#!/usr/bin/env python3
"""routes.py — the five ways a taxon gets into db_v4, resolved to taxids, as a query list.

    python 02_classify/01_refs/routes.py                 # counts + overlap, write nothing
    python 02_classify/01_refs/routes.py --write          # -> data/silver/routes.tsv

Leon's model, 2026-10-15, and it replaces the USDA-only membership of the first `01_refs`
draft. Five orthogonal routes, unioned, with no floors and no ceilings:

    1 target           a known plant pathogen            FUNGuild + PHI-base, via holobase
    2 congener         same GENUS as a known pathogen    the genus of every route-1 member
    3 host             a plant host of interest          seeds/hosts.tsv, 6 species
    4 plant_associate  recorded on any plant             seeds/fungi_usda.tsv, 1,486 taxa
    5 host_associate   seen in a phytobiome study        seeds/residents.tsv, 30 genera

The routes are better than the thing they replace in two ways beyond simplicity. **Route 2
makes LCA competition a first-class reason for membership.** db_v3 got competitors as a side
effect of an ORDER sweep from PHI-base seeds, which is why the reason a taxon was present could
not be recovered from the table. Breadth exists to give Kraken2's LCA somewhere to park a read:
db_v2 held 326 taxa, most reads had exactly one plausible match, and that is why soybean rust
was "detected" in every wheat sample. Naming that route says so out loud.

**And the route IS the label `05_detect` needs.** `pathogen | congener | host | plant_associate
| host_associate` is membership provenance and detection label at once, which is one mechanism
instead of two. It carries a hard requirement with it: **a taxon in the database only as a
congener must not be reportable as a detection the way a target is.** It is there to absorb
reads, and announcing "*Aspergillus tubingensis* found in wheat" about such a taxon would be an
artefact of our own database design. `05_detect` must read the route column.

## Why this writes a QUERY list, not a species list

Route 2 at genus covers 30,668 holobase species, of which only ~1,700 have any assembly. So
enumerating species here would be mostly fiction, and it would also be the wrong unit: one
`datasets summary genome taxon <genus>` returns the whole clade, so a genus is ONE query that
discovers its own species. Which species actually exist, and which carry a gene model, is a
question for NCBI and is answered by `survey.py`.

So each row here is a **query unit**: a taxid to ask NCBI about, plus the routes that put it
there. A unit whose ancestor is already a unit is dropped, because the ancestor's query already
returns it. That is what takes the list from ~4,400 units to something far smaller.

**Everything is resolved to a taxid and queried BY taxid.** Not a convenience: the `datasets`
CLI refuses an ambiguous NAME and the clade then returns nothing with no warning, which is how
*Cryptococcus* (fungus 5206 / scale insect 79213) and *Rhizophagus* (fungus 1129544 / beetle
295919) silently cost 42 assemblies and 17 species. A taxid cannot be ambiguous.

## What this does NOT decide

Nothing about assemblies. No annotation state, no counts, no selection. Membership is "this
taxon is of interest and here is why"; what exists for it is measured, not assumed, because
every estimate this project has made from availability counts has been too high by ~3x.
"""

import argparse
import csv
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _layout import (SILVER, LOGS, HOLOBASE, HOSTS, RESIDENTS,   # noqa: E402
                     FUNGI_USDA, EXCLUSIONS)
from _seeds import read_seed_tsv, read_seed_rows                 # noqa: E402
from _util import _Tee, make_log_dir, link_latest                # noqa: E402
import _taxonomy                                                 # noqa: E402

ROUTES_TSV = SILVER / "routes.tsv"
ROUTE_COLS = ["ncbi_taxid", "name", "rank", "kingdom", "routes", "evidence"]

# The five routes, in the order they are reported.
ROUTE_NAMES = ("target", "congener", "host", "plant_associate", "host_associate")


# ── route 1: known plant pathogens ────────────────────────────────────────────

def route_targets(tax) -> dict:
    """{internal tid: evidence} for known plant pathogens, from both sources in holobase.

    FUNGuild (`taxon_guild`) is the broad source: 2,849 taxa whose guild mentions Plant
    Pathogen. **1,646 of those are at GENUS rank and only 1,197 at species**, because FUNGuild
    assigns a guild per genus far more often than per species. That is not a defect here, it is
    most of route 2 arriving already done: expanding FUNGuild to genus adds only 212 species
    with assemblies over FUNGuild as given.

    PHI-base (`host_association` where `association_type='pathogenic'` AND the host is
    Viridiplantae) is narrow by comparison and is unioned in because it is free and curated per
    species. **The plant-host condition is load-bearing**: PHI-base scores pathogenicity
    against whatever host the experiment used, and 147 of its fungi are pathogens of Metazoa
    against 201 of plants. It is NOT the scope source any more: `pathogenic`
    reaches only the taxa somebody did a pathogenicity experiment on, which structurally cannot
    include the saprotrophs, endophytes and epiphytes routes 4 and 5 exist to reach.
    """
    out = {}
    con = sqlite3.connect(HOLOBASE)
    for tid, guild in con.execute(
            "SELECT taxon_id, guild FROM taxon_guild WHERE guild LIKE '%Plant Pathogen%'"):
        if tid in tax.rank:
            out[tid] = f"funguild:{(guild or '').strip()[:40]}"
    n_funguild = len(out)
    # The host must be a PLANT. PHI-base records pathogenicity against whatever host the
    # experiment used, and 147 of its fungi are pathogens of METAZOA -- mice, nematodes,
    # Galleria larvae, Drosophila -- against 201 with a plant host. Filtering only on the
    # FUNGUS being a fungus labelled every one of them "known plant pathogen", which put
    # *Aspergillus fumigatus*, *Candidozyma auris*, *Cryptococcus neoformans* and *Beauveria
    # bassiana* in route 1. That is not a cosmetic mislabel: the route IS the detection label,
    # so 05_detect would have reported a human clinical yeast as a plant pathogen finding.
    # Found 2026-10-16 from the contents figure, which is what the figure was for.
    # A null host kingdom is excluded too: fungi-only and plant-host-only both need a positive
    # answer, not the absence of a negative one.
    for (tid,) in con.execute(
            "SELECT DISTINCT a.fungus_id FROM host_association a "
            "JOIN taxon h ON h.taxon_id = a.host_id "
            "WHERE a.association_type='pathogenic' AND a.fungus_id IS NOT NULL "
            "AND h.kingdom = 'Viridiplantae'"):
        if tid in tax.rank and tax.kingdom.get(tid) == "Fungi":
            out[tid] = (out[tid] + "; phibase:plant_host") if tid in out else "phibase:plant_host"
    con.close()
    print(f"  route 1 target:          {n_funguild:,} FUNGuild + "
          f"{len(out) - n_funguild:,} PHI-base-only = {len(out):,} taxa")
    return out


# ── route 2: congeners of a known pathogen ────────────────────────────────────

def route_congeners(tax, targets: dict, rank: str = "genus") -> dict:
    """{internal tid: evidence} for the clade at `rank` holding each route-1 taxon.

    Measured 2026-10-15 over holobase, species with at least one assembly:

        FUNGuild as given           2,736 species / 12,161 assemblies
        expanded to GENUS           2,948 species / 13,419 assemblies
        expanded to FAMILY          4,731 species / 20,553 assemblies

    **Genus, not family.** Family adds 1,783 species for 7,134 more assemblies and reaches
    whole orders with no plant association at all, which is the breadth routes 4 and 5 already
    supply WITH evidence. `--congener-rank family` is available because the measurement should
    be reproducible, not because it is the default.

    This route deliberately includes taxa with no plant association, which is a change to the
    scope rule in `seeds/README.md` ("plant-associated, ANY host"). The change is justified
    because scope and competition are different jobs, and it has to be written into that file
    rather than left as an inconsistency.
    """
    out = {}
    unresolved = 0
    for tid in targets:
        anc = tax.ancestor_at(tid, rank)
        if anc:
            out.setdefault(anc, f"congener_of_pathogen:{rank}")
        else:
            unresolved += 1
    print(f"  route 2 congener:        {len(out):,} {rank} clades hold them"
          + (f" ({unresolved} route-1 taxa have no {rank})" if unresolved else ""))
    return out


# ── routes 3-5: the authored seeds files ──────────────────────────────────────

def route_hosts(tax) -> dict:
    out = {}
    for row in read_seed_tsv(HOSTS, require=("ncbi_taxid", "species", "role")):
        tid = tax.resolve(row["ncbi_taxid"])
        if tid:
            out[tid] = f"host:{row.get('role','').strip()}"
        else:
            print(f"    WARNING host unresolved: {row.get('species')} "
                  f"({row.get('ncbi_taxid')})")
    print(f"  route 3 host:            {len(out):,} plant species")
    return out


def route_plant_associates(tax) -> dict:
    """USDA Fungus-Host: recorded on ANY plant. The scope gate, never a plausibility filter.

    USDA is host-ASSOCIATION data and is used only to ask *is this fungus ever recorded on any
    plant*. It must never be used to ask *is it recorded on wheat, maize, rice or sorghum*: the
    study exists to find unreported host associations, and a missing reference is invisible and
    unrecoverable without a rebuild, because the reads go silently to the nearest taxon present.
    """
    out, unresolved = {}, 0
    for row in read_seed_tsv(FUNGI_USDA, require=("ncbi_taxid", "scientific_name", "rank")):
        tid = tax.resolve(row["ncbi_taxid"])
        if tid:
            out[tid] = f"usda:{row.get('n_plant_hosts','?')}_plant_hosts"
        else:
            unresolved += 1
    print(f"  route 4 plant_associate: {len(out):,} USDA taxa"
          + (f" ({unresolved} taxids not in holobase's tree)" if unresolved else ""))
    return out


def route_host_associates(tax) -> dict:
    """The phytobiome-abundance genera for wheat, rice, maize and sorghum.

    These exist because USDA is a herbarium records database, so a leaf yeast is barely
    recorded as "fungus on host X": *Papiliotrema*, *Dioszegia*, *Cystofilobasidium* and
    *Udeniomyces* appear at 1-2 records, and *Candida* and *Mrakia* not at all. The abundance
    literature says seven of the top 20 wheat genera are yeasts and *Papiliotrema* is rice's
    single most abundant genus at 14.5%. Two sources because neither alone is sufficient.

    **Both `add` and `miniprot` rows are members.** The 8 miniprot genera have assemblies but
    no gene models, so they cannot contribute CDS yet; that makes them a WORKLIST, not
    non-members. Gating membership on the action column dropped them from the list entirely,
    which is the opposite of "include everything and flag what has no annotation", and the
    8 are not marginal: *Nigrospora* is rice core at 13.2% and wheat core rank 9.06, and
    *Vishniacozyma* and *Cystofilobasidium* are both wheat core yeasts.
    """
    out, unresolved, by_action = {}, 0, defaultdict(int)
    for parts in read_seed_rows(RESIDENTS, min_cols=3):
        if len(parts) < 3:
            continue
        action = parts[2].strip()
        if action not in ("add", "miniprot"):
            continue
        by_action[action] += 1
        token = parts[0].strip()
        tid = tax.resolve(token)
        if tid:
            note = (parts[3].strip() if len(parts) > 3 else "")[:40]
            out[tid] = f"phytobiome:{action}:{note}"
        else:
            unresolved += 1
            print(f"    WARNING resident unresolved: {token!r} — put a TAXID in column 1 if "
                  f"the name is a homonym")
    print(f"  route 5 host_associate:  {len(out):,} phytobiome clades "
          f"({by_action['add']} annotated, {by_action['miniprot']} awaiting gene prediction)"
          + (f", {unresolved} unresolved" if unresolved else ""))
    return out


# ── union, exclusions, and the query-unit reduction ───────────────────────────

def enforce_kingdom(tax, members: dict, evidence: dict) -> tuple:
    """Every member must be FUNGI, except route-3 hosts which must be VIRIDIPLANTAE.

    A positive check, not a cleanup. FUNGuild's plant-pathogen set is not fungi-only:

        Fungi 2,638 | Oomycota 127 | Metazoa 66 | (null) 8 | Viridiplantae 6
        Stramenopiles 2 | Bacteria 2

    FUNGuild matches its guild assignments by NAME, so it inherits every homonym between a
    fungal genus and an animal or bacterial one — the same failure as *Cryptococcus* (fungus /
    scale insect) and *Rhizophagus* (fungus / beetle), but inside the source data rather than
    in our query. Route 2 then takes the GENUS of each, and a genus query returns everything in
    it, which is how *Ancylobacter* sp. (an alphaproteobacterium, 15 annotated assemblies)
    reached a smoke-test survey of a fungal database.

    A null kingdom is dropped too. holobase could not place it, so there is no evidence it is
    fungal, and the fungi-only decision needs a positive answer rather than the absence of a
    negative one.

    `kingdom` spells plants **Viridiplantae**, never 'Plants'; filtering on 'Plants' silently
    returns zero.
    """
    drops = defaultdict(list)
    for tid in list(members):
        want = "Viridiplantae" if "host" in members[tid] else "Fungi"
        got = tax.kingdom.get(tid)
        if got != want:
            drops[f"{got or '(null)'} (wanted {want})"].append(tax.label(tid))
            del members[tid]
    total = sum(len(v) for v in drops.values())
    print(f"\n  kingdom guard dropped {total} taxa (fungi only, hosts Viridiplantae)")
    for key, names in sorted(drops.items(), key=lambda kv: -len(kv[1])):
        shown = ", ".join(sorted(names)[:4])
        print(f"    {key:<34} {len(names):>5,}  e.g. {shown}")
    return members, drops


def load_excluded(tax) -> tuple:
    """({internal tid: kind} for the exclusion NODES, and the list of unresolvable rows).

    An exclusion must block its SUBTREE, not just the node. `4762 Oomycota` is a PHYLUM and the
    fungi-only decision of 2026-10-13 removes every oomycete, not one node that holds no
    assemblies itself. Node-only matching reported "0 blocked" and looked correct.

    Six rows carry `-` for the taxid and cannot be enforced here at all. Two are
    `synonym_already_held` (*Magnaporthe*, *Phaeosphaeria*), which name resolution handles
    properly by collapsing them onto the accepted node — a better fix than a blocklist. The
    other four are `no_assembly`, which is a statement about NCBI rather than a decision, and
    the survey will show it as zero. They are listed so the count is honest.
    """
    nodes, unenforceable = {}, []
    for row in read_seed_tsv(EXCLUSIONS, require=("ncbi_taxid", "kind")):
        kind = (row.get("kind") or "").strip() or "unspecified"
        tid = tax.resolve(row["ncbi_taxid"])
        if tid:
            nodes[tid] = kind
        else:
            unenforceable.append((row.get("name", ""), kind))
    if not nodes:
        raise ValueError(f"{EXCLUSIONS} parsed but resolved no taxids at all")
    return nodes, unenforceable


def apply_exclusions(tax, members: dict, nodes: dict) -> dict:
    """Drop every member inside an excluded subtree. Returns {dropped tid: kind}."""
    blocked = {}
    for tid in list(members):
        for ex, kind in nodes.items():
            if tid == ex or tax.is_ancestor(ex, tid):
                blocked[tid] = kind
                del members[tid]
                break
    return blocked


def reduce_to_queries(tax, members: dict) -> tuple:
    """Drop any unit whose ANCESTOR is also a unit: the ancestor's query already returns it.

    `datasets summary genome taxon <taxid>` returns the whole clade, so querying *Fusarium*
    and then *Fusarium oxysporum* separately is one wasted query and two chances to
    double-count. The reduction is done by walking each unit's lineage against the unit set
    rather than by comparing ranks, so it is correct for a forma specialis under a species
    under a genus as well.

    Returns (query units, {dropped unit: the ancestor covering it}).
    """
    unit_set = set(members)
    kept, covered = {}, {}
    for tid in members:
        anc, cur, seen = 0, tax.parent.get(tid), set()
        while cur is not None and cur not in seen:
            seen.add(cur)
            if cur in unit_set:
                anc = cur
                break
            nxt = tax.parent.get(cur)
            if nxt == cur:
                break
            cur = nxt
        if anc:
            covered[tid] = anc
        else:
            kept[tid] = members[tid]
    return kept, covered


def build(congener_rank: str = "genus") -> tuple:
    tax = _taxonomy.load()
    print(f"holobase: {len(tax.rank):,} taxa, {len(tax.by_ncbi):,} with an NCBI taxid, "
          f"{len(tax.by_name):,} names incl. synonyms\n")

    per_route = {}
    targets = route_targets(tax)
    per_route["target"] = targets
    per_route["congener"] = route_congeners(tax, targets, congener_rank)
    per_route["host"] = route_hosts(tax)
    per_route["plant_associate"] = route_plant_associates(tax)
    per_route["host_associate"] = route_host_associates(tax)

    # Union, keeping EVERY route that applies: a taxon is usually in by more than one, and
    # which ones decide how a detection on it may be reported.
    members, evidence = defaultdict(list), defaultdict(list)
    for name in ROUTE_NAMES:
        for tid, ev in per_route[name].items():
            members[tid].append(name)
            evidence[tid].append(ev)

    members, kingdom_drops = enforce_kingdom(tax, dict(members), evidence)

    ex_nodes, unenforceable = load_excluded(tax)
    blocked = apply_exclusions(tax, members, ex_nodes)
    by_kind = defaultdict(int)
    for kind in blocked.values():
        by_kind[kind] += 1
    print(f"\n  exclusions.tsv: {len(ex_nodes)} enforceable rules blocked {len(blocked)} taxa "
          f"(subtrees included)")
    for kind, n in sorted(by_kind.items(), key=lambda kv: -kv[1]):
        print(f"    {kind:<24} {n:>6,}")
    if unenforceable:
        print(f"    {len(unenforceable)} rows carry no taxid and cannot be enforced by id: "
              + ", ".join(f"{n} ({k})" for n, k in unenforceable))

    print(f"\n── Union ────────────────────────────────────────────────────────")
    print(f"  {len(members):,} distinct taxa across the five routes")
    overlap = defaultdict(int)
    for routes in members.values():
        overlap[len(routes)] += 1
    for n in sorted(overlap):
        print(f"    in {n} route{'s' if n > 1 else ' '}: {overlap[n]:>6,} taxa")
    print(f"  by rank:")
    rank_census = defaultdict(int)
    for tid in members:
        rank_census[tax.rank.get(tid, "?")] += 1
    for rk, n in sorted(rank_census.items(), key=lambda kv: -kv[1])[:8]:
        print(f"    {rk:<16} {n:>6,}")

    queries, covered = reduce_to_queries(tax, dict(members))
    print(f"\n── Query units ──────────────────────────────────────────────────")
    print(f"  {len(queries):,} taxa to ask NCBI about")
    print(f"  {len(covered):,} dropped: an ancestor in the list already returns them")
    print(f"  every one is queried BY TAXID, so a homonym cannot return zero silently")

    route_of = {tid: ";".join(r for r in ROUTE_NAMES if r in routes)
                for tid, routes in members.items()}
    rows = []
    for tid, routes in sorted(members.items(), key=lambda kv: tax.label(kv[0])):
        rows.append({
            "ncbi_taxid": tax.ncbi_of(tid),
            "name": tax.label(tid),
            "rank": tax.rank.get(tid, ""),
            "kingdom": tax.kingdom.get(tid) or "",
            "routes": route_of[tid],
            "evidence": " | ".join(evidence[tid])[:300],
        })
    return rows, queries, covered, tax, route_of


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--congener-rank", default="genus", choices=("genus", "family"),
                    help="rank at which route 2 expands. genus by default; family adds 1,783 "
                         "species and reaches clades with no plant association")
    ap.add_argument("--write", action="store_true",
                    help=f"write {ROUTES_TSV}. Without it, report and write nothing")
    args = ap.parse_args()

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "routes.log")
    link_latest(LOGS, log_dir / "routes.log")
    sys.stdout = log
    try:
        rows, queries, covered, tax, route_of = build(args.congener_rank)
        if not args.write:
            print(f"\n(--write not given: nothing written)")
            return
        ROUTES_TSV.parent.mkdir(parents=True, exist_ok=True)
        with open(ROUTES_TSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=ROUTE_COLS, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        qpath = SILVER / "route_queries.tsv"
        with open(qpath, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["ncbi_taxid", "name", "rank", "routes"])
            for tid in sorted(queries, key=lambda t: tax.label(t)):
                w.writerow([tax.ncbi_of(tid), tax.label(tid), tax.rank.get(tid, ""),
                            route_of.get(tid, "")])
        print(f"\n  -> {ROUTES_TSV}  ({len(rows):,} taxa)")
        print(f"  -> {qpath}  ({len(queries):,} query units)")
        print(f"\nNext: python 02_classify/01_refs/survey.py   "
              f"(asks NCBI what assemblies and annotations actually exist)")
    finally:
        log.close()


if __name__ == "__main__":
    main()
