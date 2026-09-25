#!/usr/bin/env python3
"""
kraken_db_search.py — select and download candidate reference assemblies for the
Kraken2 pathogen DB. Submodule 1, step 1 of 3 (search → busco → build).

Three modes, run in order:

  --scope     Optional reconnaissance report before selecting candidates: how many
              assemblies exist at NCBI per PHI-base seed (pan-genome) and per genus
              (fill-in headroom). Informational only — doesn't affect selection.
              Output: data/pangenome.tsv, data/genus_fill.tsv

  (default)   SEED (pan-genome): for each PHI-base euk seed taxid, collect ALL
              scaffold-plus+annotated assemblies. Detected-genera seeds get all of
              them ordered by geographic/temporal diversity; undetected-genera seeds
              get --floor assemblies (db_v2 took one). Falls back through
              scaffold-plus, then contig-level — always take something over nothing.

              RANK SWEEP (breadth): for every taxon at --sweep-rank containing a
              seed species, add up to --floor annotated scaffold-plus assemblies for
              each species in it not already covered by a seed. db_v2 swept at genus
              with no floor; db_v3 sweeps at order.

              Candidate taxids are normalised to their SPECIES node first. NCBI
              labels many assemblies at strain, varietas or forma rank, and grouping
              on the assembly's own taxid both re-selects a seed's own strains and
              splits one species across several nodes so the floor passes while the
              species stays shallow. See select_rank_sweep.

              --retain guarantees no assembly in the previous selection is dropped,
              which is what makes the no-cap policy below hold: min(floor, available)
              is a cap for any species the sweep owns, and db_v2 took non-seed species
              with no limit.

              Output: data/ref_candidates_<version>.tsv, plus a ref_candidates.tsv
              symlink to it so kraken_db_busco.py picks up the current version
              without the previous version's table being overwritten.

Depth policy (measured, see kraken_db_depth_policy.py):

  A FLOOR, not a cap. Misassignment between congeners is driven by the DIFFERENCE in
  their sampling depth: a k-mer in both species' pangenomes but absent from the
  shallow one's included assemblies is stored as the deep species, so the shallow
  species' reads score for the deep one. Over the 26 congener pairs above 1%
  misassignment, raising the shallow species to 4 assemblies was worth a median
  1.70 pp against 0.42 pp for capping the deep species at 5, and unlike the cap it
  costs nothing: it only adds assemblies. Pangenomes are open (Heaps gamma ~0.16,
  no asymptote), so a cap discards real novel sequence — up to 62% of a species
  complex's observed k-mer content. Hence --floor with no cap.

  The floor is a target, not a filter. A species with one assembly in existence is
  included at one; `min(floor, available)`. It still earns its place by giving the
  LCA a competitor, which is what stops a read being handed to whichever species
  happens to be in the build, but its own species-level counts are not reliable and
  should be aggregated to genus when reported. At order scope that applies to about
  70% of the swept species.

  --download  Download CDS FASTA for every candidate with fasta_type=cds into
              data/cds/pathogen/{accession}/. This is the ONLY place in the kraken_db_*
              pipeline that downloads CDS — kraken_db_busco.py and
              kraken_db_build.py both read from here, never fetch their own copies.
              Resumable: skips accessions with .fna files already on disk.

Typical run (select + download in one pass):
    python kraken/search.py --download

Reproduce db_v2's composition exactly (genus sweep, no floor, broad genera excluded):
    python kraken/search.py --version v2 --sweep-rank genus \
        --floor 0 --exclude-broad --select-only

Then: python kraken/utilities/busco.py   (on Setonix)

Output:
    kraken/search/data/ref_candidates.tsv   (tracked)
    kraken/search/data/pangenome.tsv        (tracked, --scope only)
    kraken/search/data/genus_fill.tsv       (tracked, --scope only)
    kraken/search/data/cds/pathogen/{accession}/  (gitignored, --download only)
"""

import argparse
import csv
import gzip
import json
import shutil
import sqlite3
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _util import _Tee, make_log_dir, link_latest

DB_PATH      = Path("stat/build/data/phibase_db.json")
RUNS_TSV     = Path("stat/filter/data/runs.tsv")
OUT_DIR      = Path("kraken/search")
DATA_DIR     = OUT_DIR / "data"
CANDIDATES_TSV = DATA_DIR / "ref_candidates.tsv"      # symlink to the current version
DEFAULT_GENOMES_DIR = DATA_DIR / "cds" / "pathogen"

# holobase supplies the taxonomy tree used to resolve each seed species to its
# containing order/family/class. NCBI's own assembly records carry no rank above
# species, so the sweep targets cannot be derived from the datasets output alone.
HOLOBASE = Path.home() / "data_analysis" / "_refdata" / "holobase.db"

DEFAULT_VERSION    = "v3"
DEFAULT_FLOOR      = 4        # 0 = no floor, take everything (db_v2's behaviour)
DEFAULT_SWEEP_RANK = "order"
SWEEP_RANKS        = ("genus", "family", "order", "class")

YEAR_BIN_SIZE = 5
LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

# Genera excluded from the sweep by --exclude-broad: broad/saprophytic, not primarily
# plant pathogens. db_v2 always excluded these. db_v3 does NOT, by default: the sweep
# exists to give Kraken2's LCA competitors, and a saprophyte in a target order is as
# good a competitor as a pathogen. The floor bounds what their inclusion costs.
BROAD_GENERA = {
    "Aspergillus", "Penicillium", "Trichoderma", "Beauveria", "Metarhizium",
    "Claviceps", "Epichloe", "Ceratocystis", "Leptographium", "Ciboria",
}

# Basidiomycete PHI-base genera → basidiomycota_odb10
BASIDIOMYCETE_GENERA = {
    "Puccinia", "Melampsora", "Phakopsora", "Hemileia", "Uromyces",
    "Tranzschelia", "Phragmidium", "Gymnosporangium",        # rusts
    "Ustilago", "Tilletia", "Sporisorium", "Testicularia",   # smuts
    "Mycosarcoma",                                            # Ustilaginomycotina
    "Rhizoctonia", "Heterobasidion", "Moniliophthora", "Crinipellis",
}

# Chytrid/other early-diverging fungi with no phylum-level BUSCO lineage → fungi_odb10
CHYTRID_GENERA = {"Synchytrium"}

CANDIDATE_COLS = [
    "taxid", "organism_name", "kingdom", "source",
    "accession", "assembly_level", "release_date", "has_annotation",
    "country", "protein_coding_genes", "scaffold_n50_kb", "total_length_mb",
    "fasta_type", "busco_lineage", "selection_rank", "selection_reason",
]


def _ts() -> str:
    return time.strftime("%H:%M:%S")


# ── NCBI datasets query ───────────────────────────────────────────────────────

def datasets_query(taxon) -> list:
    r = subprocess.run(
        ["datasets", "summary", "genome", "taxon", str(taxon),
         "--as-json-lines", "--limit", "all"],
        capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        return []
    rows = []
    for line in r.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if "reports" in d:
            rows.extend(d["reports"])
        elif "accession" in d:
            rows.append(d)
    return rows


# ── assembly field helpers ────────────────────────────────────────────────────

def _ai(a): return a.get("assembly_info") or {}
def _ann(a): return a.get("annotation_info") or {}
def _stats(a): return a.get("assembly_stats") or {}

def assembly_level(a): return _ai(a).get("assembly_level", "")
def release_date(a): return _ai(a).get("release_date", "") or ""
def has_annotation(a): return bool(_ann(a))

def get_country(a) -> str:
    for attr in (_ai(a).get("biosample") or {}).get("attributes", []) or []:
        if attr.get("name") == "geo_loc_name":
            v = (attr.get("value") or "").strip()
            if v and v.lower() not in ("not applicable", "missing", "not collected", ""):
                return v.split(":")[0].strip()
    return ""

def get_year_bin(a) -> str:
    y = release_date(a)[:4]
    try:
        return str((int(y) // YEAR_BIN_SIZE) * YEAR_BIN_SIZE)
    except ValueError:
        return ""

def level_rank(a) -> int:
    return LEVEL_RANK.get(assembly_level(a), 0)

def quality_key(a) -> tuple:
    ann   = _ann(a)
    genes = ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", 0) or 0 if ann else 0
    acc   = a.get("accession", "")
    return (
        1 if ann else 0,
        min(int(genes), 50000) // 5000,
        level_rank(a),
        1 if acc.startswith("GCF_") else 0,
        release_date(a),
    )

def scaffold_plus(a) -> bool:
    return level_rank(a) >= 2

# Phylum -> BUSCO lineage. Populated from holobase by holobase_phylum_map(); the genus
# sets above are the fallback when holobase is unavailable.
PHYLUM_LINEAGE = {
    "Basidiomycota":   "basidiomycota_odb10",
    "Ascomycota":      "ascomycota_odb10",
    "Chytridiomycota": "fungi_odb10",
    "Blastocladiomycota": "fungi_odb10",
    "Mucoromycota":    "fungi_odb10",
    "Zoopagomycota":   "fungi_odb10",
    "Oomycota":        "stramenopiles_odb10",
}
PHYLUM_OF = {}


def busco_lineage_for(kingdom: str, organism_name: str, taxid=None) -> str:
    """Pick the BUSCO lineage database for an assembly.

    Taxonomy first, genus lists second. BASIDIOMYCETE_GENERA holds 17 genera, all
    pathogens, which covered db_v2 because db_v2 only ever saw pathogens. An
    order-level sweep does not stay inside it: 160 of db_v3's 228 Basidiomycota
    candidates fall in 63 genera outside the list (Lentinula, Tulasnella,
    Ceratobasidium, Amanita, Agaricus …) and would be scored against
    ascomycota_odb10, producing spuriously low completeness for exactly the breadth
    the sweep was added to provide.
    """
    if kingdom == "oomycete":
        return "stramenopiles_odb10"
    if taxid is not None:
        lineage = PHYLUM_LINEAGE.get(PHYLUM_OF.get(int(taxid)))
        if lineage:
            return lineage
    genus = organism_name.split()[0]
    if genus in BASIDIOMYCETE_GENERA:
        return "basidiomycota_odb10"
    if genus in CHYTRID_GENERA:
        return "fungi_odb10"
    return "ascomycota_odb10"


def holobase_phylum_map(db_path: Path = HOLOBASE) -> dict:
    """{ncbi taxid: phylum name}, by walking holobase's tree up to phylum rank.

    Returns {} if holobase is absent, so busco_lineage_for falls back to the genus
    sets. Guards against NCBI's self-parent root, as holobase_species_map does.
    """
    if not db_path.exists():
        print(f"Warning: {db_path} not found — BUSCO lineages fall back to genus lists")
        return {}
    con = sqlite3.connect(db_path)
    parent, rank, ncbi, name = {}, {}, {}, {}
    for tid, par, rk, nc, nm in con.execute(
            "SELECT taxon_id, parent_id, rank, ncbi_taxid, scientific_name FROM taxon"):
        parent[tid], rank[tid] = par, rk
        if nc is not None:
            ncbi[tid] = nc
        if rk == "phylum":
            name[tid] = nm
    con.close()

    cache, out = {}, {}
    for nc, tid in ((nc, tid) for tid, nc in ncbi.items()):
        start, chain, found = tid, [], None
        while start is not None and start in rank:
            if start in cache:
                found = cache[start]
                break
            if rank[start] == "phylum":
                found = name.get(start)
                break
            chain.append(start)
            nxt = parent.get(start)
            if nxt == start or nxt in chain:
                break
            start = nxt
        for t in chain:
            cache[t] = found
        if found:
            out[nc] = found
    return out

def to_row(a: dict, taxid: int, organism_name: str, kingdom: str,
           source: str, rank: int, reason: str) -> dict:
    ann  = _ann(a)
    stat = _stats(a)
    acc  = a.get("accession", "")
    genes = ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", "") if ann else ""
    n50   = stat.get("scaffold_n50") or stat.get("contig_n50", "")
    total = stat.get("total_sequence_length", "")
    return {
        "taxid":               taxid,
        "organism_name":       organism_name,
        "kingdom":             kingdom,
        "source":              source,
        "accession":           acc,
        "assembly_level":      assembly_level(a),
        "release_date":        release_date(a),
        "has_annotation":      has_annotation(a),
        "country":             get_country(a),
        "protein_coding_genes": str(genes) if genes else "",
        "scaffold_n50_kb":     f"{int(n50)/1000:.0f}" if n50 else "",
        "total_length_mb":     f"{int(total)/1e6:.1f}" if total else "",
        "fasta_type":          "cds" if ann else "genome",
        "busco_lineage":       busco_lineage_for(kingdom, organism_name, taxid),
        "selection_rank":      rank,
        "selection_reason":    reason,
    }


# ── diversity-ordered selection (no cap) ─────────────────────────────────────

def greedy_ordered(assemblies: list) -> list:
    """Order all assemblies by geographic + temporal diversity.
    assemblies must already be sorted by quality_key descending.
    Returns list of (assembly, reason_str) — all assemblies, no cap."""
    if not assemblies:
        return []
    remaining = list(assemblies)
    selected, reasons = [], []
    selected.append(remaining.pop(0))
    reasons.append("best_quality")
    while remaining:
        sel_countries = {get_country(a) for a in selected}
        sel_bins      = {get_year_bin(a) for a in selected}
        best_score, best_idx, best_parts = -1, 0, []
        for i, a in enumerate(remaining):
            country, ybin = get_country(a), get_year_bin(a)
            score, parts = 0, []
            if country and country not in sel_countries:
                score += 2
                parts.append(f"new_country:{country}")
            if ybin and ybin not in sel_bins:
                score += 1
                parts.append(f"new_year_bin:{ybin}")
            if score > best_score:
                best_score, best_idx, best_parts = score, i, parts
        selected.append(remaining.pop(best_idx))
        reasons.append("; ".join(best_parts) if best_parts else "quality_fill")
    return list(zip(selected, reasons))


def select_seed(taxid: int, name: str, kingdom: str, detected: bool,
                assemblies: list, floor: int = DEFAULT_FLOOR) -> list:
    """detected=True → all scaffold-plus+annotated assemblies, diversity-ordered, no
    cap. detected=False → min(floor, available), diversity-ordered. Falls back through
    annotation tiers if needed.

    db_v2 took exactly one assembly for an undetected genus, which is why 170 of its
    309 taxa sat at depth 1 and were the ones losing k-mers to deeper congeners. The
    floor closes that without touching the deep species.
    """
    for pool in [
        [a for a in assemblies if scaffold_plus(a) and has_annotation(a)],
        [a for a in assemblies if scaffold_plus(a)],
        assemblies,
    ]:
        if pool:
            break
    if not pool:
        return []
    pool.sort(key=quality_key, reverse=True)
    pairs = greedy_ordered(pool)
    if not detected and floor > 0:
        pairs = pairs[:floor]
    return [to_row(a, taxid, name, kingdom, "seed", rank + 1, reason)
            for rank, (a, reason) in enumerate(pairs)]


def select_rank_sweep(rank_name: str, kingdom: str, covered_taxids: set,
                      assemblies: list, floor: int = DEFAULT_FLOOR,
                      exclude_broad: bool = False,
                      to_species: dict = None) -> list:
    """Sweep one taxon at --sweep-rank: for each species in it not already covered by
    a PHI-base seed, take min(floor, available) annotated scaffold-plus assemblies,
    diversity-ordered within the species. floor=0 takes every assembly.

    One NCBI query covers the whole clade, so an order sweep costs 20 queries rather
    than one per species. Species with a single assembly are kept at one, not dropped:
    presence is what creates the LCA competition the sweep exists for.

    NCBI labels many assemblies at a node BELOW species — strain, varietas, forma —
    so grouping on the assembly's own tax_id is wrong in two ways that both defeat
    the design. `to_species` maps each NCBI taxid up to its species node and fixes
    both:

      1. A seed's own strains stopped being recognised as covered, because the
         covered check compared taxids for equality and a strain taxid never equals
         its species'. That re-selected assemblies already taken as seeds — 159 of
         them at order scope, e.g. GCF_014117465.1 as both *Aspergillus flavus*
         (5059) and *A. flavus* NRRL3357 (332952).
      2. The floor counted per taxid, so a species whose assemblies sit under
         several strain nodes looked like several depth-1 taxa and the floor passed
         while the species stayed shallow. 87 species were split this way, worst
         *Zymoseptoria tritici* across 6 nodes and *Verticillium dahliae* across 4.

    Rows are emitted under the species taxid, which is the rank the chapter reports
    at and the rank db_v2 used for its seeds.
    """
    to_species = to_species or {}
    by_species = defaultdict(list)
    for a in assemblies:
        organism = a.get("organism") or {}
        tid = organism.get("tax_id")
        if not tid:
            continue
        # Normalise to the species node before either check below.
        sid = to_species.get(tid, tid)
        if sid in covered_taxids or tid in covered_taxids:
            continue
        if exclude_broad and (organism.get("organism_name") or "").split()[:1] \
                and (organism.get("organism_name") or "").split()[0] in BROAD_GENERA:
            continue
        if scaffold_plus(a) and has_annotation(a):
            by_species[sid].append(a)

    rows = []
    for sid, pool in by_species.items():
        pool.sort(key=quality_key, reverse=True)
        # Prefer a name recorded at species rank; strain names carry isolate suffixes
        # that would make the same species look like several.
        org_name = SPECIES_NAME.get(sid) or \
            (pool[0].get("organism") or {}).get("organism_name", str(sid))
        pairs = greedy_ordered(pool)
        if floor > 0:
            pairs = pairs[:floor]
        for rank, (a, reason) in enumerate(pairs):
            rows.append(to_row(a, sid, org_name, kingdom,
                               f"{rank_name}_sweep", rank + 1, reason))
    return rows


def holobase_sweep_targets(seed_taxids, rank: str, db_path: Path = HOLOBASE) -> dict:
    """{clade name: n seed species it contains} for every clade at `rank` holding a
    seed species. Walks holobase's taxonomy tree, which is keyed on an internal
    surrogate taxon_id, so seeds are joined in through the ncbi_taxid column.

    Returns {} if holobase is absent, so the caller can fall back to a genus sweep
    derived from the seed names alone.
    """
    if not db_path.exists():
        print(f"Warning: {db_path} not found — cannot resolve {rank}-level sweep targets")
        return {}

    con = sqlite3.connect(db_path)
    node = {tid: (parent, rk, nm) for tid, parent, rk, nm
            in con.execute("SELECT taxon_id, parent_id, rank, scientific_name FROM taxon")}
    ncbi2internal = dict(con.execute(
        "SELECT ncbi_taxid, taxon_id FROM taxon WHERE ncbi_taxid IS NOT NULL"))
    con.close()

    def ancestor_at(tid):
        seen = set()
        while tid in node and tid not in seen:
            seen.add(tid)
            parent, rk, name = node[tid]
            if rk == rank:
                return name
            tid = parent
        return None

    targets, unresolved = defaultdict(int), 0
    for taxid in seed_taxids:
        internal = ncbi2internal.get(int(taxid))
        name = ancestor_at(internal) if internal else None
        if name:
            targets[name] += 1
        else:
            unresolved += 1
    if unresolved:
        print(f"  {unresolved} seed taxids had no {rank} in holobase — not swept")
    return dict(targets)


def holobase_species_map(db_path: Path = HOLOBASE) -> tuple:
    """Build ({ncbi taxid: species-level ncbi taxid}, {species taxid: species name}).

    NCBI assigns many assemblies to a node below species. Kraken2 would happily build
    from those, but the depth floor and the seed-coverage check both reason per
    species, so both need the species node. Species taxids map to themselves, so the
    caller can apply `to_species.get(tid, tid)` unconditionally.

    Returns ({}, {}) if holobase is absent; the sweep then falls back to grouping on
    the assembly's own taxid, which is db_v2's behaviour and its bug.
    """
    if not db_path.exists():
        print(f"Warning: {db_path} not found — cannot normalise taxids to species rank")
        return {}, {}

    con = sqlite3.connect(db_path)
    # The taxon table is ~3M rows, so hold only what the walk needs: parent and rank
    # for every node, but names for species nodes alone.
    parent, rank, ncbi, name = {}, {}, {}, {}
    for tid, par, rk, nc, nm in con.execute(
            "SELECT taxon_id, parent_id, rank, ncbi_taxid, scientific_name FROM taxon"):
        parent[tid], rank[tid] = par, rk
        if nc is not None:
            ncbi[tid] = nc
        if rk == "species":
            name[tid] = nm
    con.close()

    internal_of = {nc: tid for tid, nc in ncbi.items()}
    to_species, species_name, cache = {}, {}, {}

    def species_internal(start):
        """Nearest ancestor-or-self at species rank, memoised over shared lineages.

        NCBI's root is its own parent, so an unbounded walk never terminates for a
        taxon with no species ancestor (a genus, a family, root itself). `chain`
        doubles as the visited set to stop that.
        """
        tid, chain, found = start, [], None
        while tid is not None and tid in rank:
            if tid in cache:
                found = cache[tid]
                break
            if rank[tid] == "species":
                found = tid
                break
            chain.append(tid)
            nxt = parent.get(tid)
            if nxt == tid or nxt in chain:      # self-parent root, or a cycle
                break
            tid = nxt
        for t in chain:
            cache[t] = found
        return found

    for nc, tid in internal_of.items():
        sp = species_internal(tid)
        if sp is None or sp not in ncbi:
            continue
        to_species[nc] = ncbi[sp]
        species_name[ncbi[sp]] = name[sp]
    return to_species, species_name


# Filled in by run_select once holobase has been read, so select_rank_sweep can name a
# species without carrying the whole map through every call.
SPECIES_NAME = {}


def retain_previous(rows: list, prev_tsv: Path) -> list:
    """Add back any assembly in a previous candidate table that the new selection
    dropped, so a rebuild never shrinks a species' depth.

    The depth policy measured on 2026-09-24 is a floor with NO cap: capping a deep
    species buys a median 0.42 pp of congener specificity and costs a median 16% of
    its own observed pangenome, up to 62%. But the sweep's `min(floor, available)` is
    a cap for any species the sweep owns, and db_v2 took non-seed species through
    genus fill with no limit. Without this, a rebuild silently trimmed 45 assemblies,
    including *Fusarium* cf. *solani* 14 -> 4 and *F. lateritium* 15 -> 4, both deep
    partners in the worst-measured misassignment pairs.

    Retained rows keep their original taxid, source and reason, tagged so the table
    shows they came from the previous build rather than this selection's rules.
    """
    if not prev_tsv.exists():
        return rows
    have = {r["accession"] for r in rows}
    added = []
    with open(prev_tsv) as fh:
        for prev in csv.DictReader(fh, delimiter="\t"):
            acc = prev.get("accession")
            if not acc or acc in have:
                continue
            have.add(acc)
            kept = dict(prev)
            # Rows read from TSV carry every field as a string, but rows built by
            # to_row() carry taxid as an int. Left mixed, the per-taxon counters treat
            # 5507 and "5507" as two taxa and the selection summary overcounts.
            if str(kept.get("taxid", "")).isdigit():
                kept["taxid"] = int(kept["taxid"])
            kept["selection_reason"] = (
                f"retained_previous:{prev.get('selection_reason','') or 'na'}")
            added.append(kept)
    if added:
        by_taxon = defaultdict(int)
        for r in added:
            by_taxon[r.get("organism_name", "?")] += 1
        top = sorted(by_taxon.items(), key=lambda kv: -kv[1])[:5]
        print(f"  retained {len(added)} assemblies the new rules would have dropped "
              f"(no-cap policy), across {len(by_taxon)} taxa")
        for nm, n in top:
            print(f"      +{n:<3} {nm}")
    return rows + added


def load_detected_genera() -> set:
    genera = set()
    try:
        with open(RUNS_TSV) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                for part in (row.get("stat_pathogens") or "").split(";"):
                    name = part.strip().split(":")[0].strip()
                    if name:
                        genera.add(name.split()[0])
    except FileNotFoundError:
        print(f"Warning: {RUNS_TSV} not found — treating all genera as detected")
    return genera


# ── --scope reconnaissance (from legacy scope_db.py) ─────────────────────────

def summarise_assemblies(rows: list) -> dict:
    if not rows:
        return {"n_total": 0, "n_scaffold_plus": 0, "n_annotated": 0,
                "date_min": "", "date_max": "", "n_countries": 0,
                "countries": "", "best_level": ""}
    dates = sorted(d for a in rows if (d := release_date(a)))
    countries = sorted({c for a in rows if (c := get_country(a))})
    scaffold_up = [a for a in rows if scaffold_plus(a)]
    annotated = [a for a in rows if has_annotation(a)]
    levels = [assembly_level(a) for a in rows]
    best = max(levels, key=lambda l: LEVEL_RANK.get(l, 0), default="")
    return {"n_total": len(rows), "n_scaffold_plus": len(scaffold_up),
            "n_annotated": len(annotated), "date_min": dates[0] if dates else "",
            "date_max": dates[-1] if dates else "", "n_countries": len(countries),
            "countries": "; ".join(countries[:20]), "best_level": best}


def run_scope(seeds: dict, t2n: dict, workers: int) -> None:
    print(f"[--scope] Seeds: {len(seeds)}")

    print("\n[1/2] Querying assemblies per seed taxid …")
    pan_rows = []

    def query_seed(tid, kingdom):
        name = t2n.get(str(tid), str(tid))
        return {"taxid": tid, "name": name, "kingdom": kingdom,
                **summarise_assemblies(datasets_query(tid))}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(query_seed, tid, kgd): tid for tid, kgd in seeds.items()}
        for done, fut in enumerate(as_completed(futs), 1):
            pan_rows.append(fut.result())
            if done % 10 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} seeds …", end="\r", flush=True)
    print()

    pan_rows.sort(key=lambda r: (-r["n_total"], r["name"]))
    pan_path = DATA_DIR / "pangenome.tsv"
    pan_cols = ["taxid", "name", "kingdom", "n_total", "n_scaffold_plus", "n_annotated",
                "best_level", "date_min", "date_max", "n_countries", "countries"]
    with open(pan_path, "w") as fh:
        fh.write("\t".join(pan_cols) + "\n")
        for r in pan_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in pan_cols) + "\n")
    print(f"Written: {pan_path}")

    print("\n[2/2] Querying genus-level fill-in …")
    seed_taxids_set = set(seeds.keys())
    genus_to_seeds = defaultdict(list)
    for tid in seeds:
        parts = t2n.get(str(tid), "").split()
        if len(parts) >= 2:
            genus_to_seeds[parts[0]].append(tid)

    fill_rows = []

    def query_genus(genus, seed_tids):
        rows = datasets_query(genus)
        if not rows:
            return None
        new_species = defaultdict(list)
        phibase_n = 0
        for a in rows:
            tid = a.get("organism", {}).get("tax_id")
            if tid in seed_taxids_set:
                phibase_n += 1
            elif tid:
                new_species[tid].append(a)
        new_sp_names = sorted({
            a.get("organism", {}).get("organism_name", "") for a in rows
            if a.get("organism", {}).get("tax_id") not in seed_taxids_set
            and a.get("organism", {}).get("tax_id")
        })
        return {
            "genus": genus, "n_phibase_seeds": len(seed_tids),
            "n_phibase_assemblies": phibase_n, "n_new_species": len(new_species),
            "n_new_assemblies": sum(len(v) for v in new_species.values()),
            "n_new_scaffold": sum(1 for v in new_species.values() for a in v if scaffold_plus(a)),
            "n_new_annotated": sum(1 for v in new_species.values() for a in v if has_annotation(a)),
            "new_species": "; ".join(new_sp_names[:30]),
        }

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(query_genus, g, s): g for g, s in genus_to_seeds.items()}
        for done, fut in enumerate(as_completed(futs), 1):
            row = fut.result()
            if row:
                fill_rows.append(row)
            if done % 5 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} genera …", end="\r", flush=True)
    print()

    fill_rows.sort(key=lambda r: -r["n_new_assemblies"])
    fill_path = DATA_DIR / "genus_fill.tsv"
    fill_cols = ["genus", "n_phibase_seeds", "n_phibase_assemblies", "n_new_species",
                 "n_new_assemblies", "n_new_scaffold", "n_new_annotated", "new_species"]
    with open(fill_path, "w") as fh:
        fh.write("\t".join(fill_cols) + "\n")
        for r in fill_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in fill_cols) + "\n")
    print(f"Written: {fill_path}")


# ── candidate selection ───────────────────────────────────────────────────────

def run_select(seeds: dict, t2n: dict, workers: int, floor: int = DEFAULT_FLOOR,
               sweep_rank: str = DEFAULT_SWEEP_RANK,
               exclude_broad: bool = False, retain: Path = None) -> list:
    global SPECIES_NAME, PHYLUM_OF
    detected_genera = load_detected_genera()
    print(f"STAT-detected genera: {len(detected_genera)}")
    print(f"Depth floor: {floor or 'none'}   sweep rank: {sweep_rank}"
          f"   broad genera: {'excluded' if exclude_broad else 'included'}")

    # Taxid -> species taxid, so the sweep's coverage check and depth floor both reason
    # per species rather than per strain. See select_rank_sweep for what breaks without it.
    to_species, SPECIES_NAME = holobase_species_map()
    PHYLUM_OF = holobase_phylum_map()
    covered = set(seeds.keys())
    if to_species:
        # A seed given at strain or forma-specialis rank still covers its species, and
        # vice versa, so expand coverage both ways before the sweep runs.
        covered |= {to_species[t] for t in seeds if t in to_species}
        print(f"Species-normalised {len(to_species):,} taxids; "
              f"seed coverage expanded {len(seeds)} -> {len(covered)} taxa")

    print(f"\n[1/2] Querying {len(seeds)} seed taxids (all assemblies for detected genera) …")
    all_rows = []

    def work_seed(taxid, kingdom):
        name = t2n.get(str(taxid), str(taxid))
        detected = name.split()[0] in detected_genera
        return select_seed(taxid, name, kingdom, detected, datasets_query(taxid), floor)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work_seed, tid, kgd): tid for tid, kgd in seeds.items()}
        for done, fut in enumerate(as_completed(futs), 1):
            all_rows.extend(fut.result())
            if done % 10 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} seeds … ({len(all_rows)} assemblies so far)",
                      end="\r", flush=True)
    print()
    seed_count = len(all_rows)
    print(f"  Seed assemblies selected: {seed_count}")

    # Sweep targets: the clades at sweep_rank that actually contain a seed species.
    # Above genus that needs a taxonomy tree, which holobase has and NCBI's assembly
    # records do not (they carry no rank above species). A genus sweep needs no tree,
    # so it stays available and reproduces db_v2 with --floor 0 --exclude-broad.
    if sweep_rank == "genus":
        sweep_targets = {}
        for tid, kingdom in seeds.items():
            parts = t2n.get(str(tid), "").split()
            if parts and parts[0] in detected_genera:
                sweep_targets[parts[0]] = kingdom
    else:
        clades = holobase_sweep_targets(covered, sweep_rank)
        if not clades:
            print(f"  No {sweep_rank} targets resolved — falling back to a genus sweep")
            sweep_targets = {parts[0]: kingdom
                             for tid, kingdom in seeds.items()
                             if (parts := t2n.get(str(tid), "").split())
                             and parts[0] in detected_genera}
            sweep_rank = "genus"
        else:
            # Kingdom only drives the BUSCO lineage call, and busco_lineage_for()
            # refines it per species anyway. Default the clade to fungal, then mark
            # the clades that oomycete seeds put in scope — holobase holds no oomycete
            # assemblies, so those clades come back from NCBI or not at all.
            sweep_targets = dict.fromkeys(clades, "fungal")
            oomycete_seeds = {t for t, k in seeds.items() if k == "oomycete"}
            if oomycete_seeds:
                for name in holobase_sweep_targets(oomycete_seeds, sweep_rank):
                    sweep_targets[name] = "oomycete"
            print(f"  {len(clades)} {sweep_rank}-level clades hold the "
                  f"{len(covered)} seed species")

    if exclude_broad and sweep_rank == "genus":
        sweep_targets = {g: k for g, k in sweep_targets.items() if g not in BROAD_GENERA}

    label = "all qualifying assemblies" if not floor else f"up to {floor} per species"
    print(f"\n[2/2] Sweeping {len(sweep_targets)} {sweep_rank} clades ({label}) …")

    def work_clade(name, kingdom):
        return select_rank_sweep(sweep_rank, kingdom, covered, datasets_query(name),
                                 floor, exclude_broad, to_species)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work_clade, n, k): n for n, k in sweep_targets.items()}
        for done, fut in enumerate(as_completed(futs), 1):
            all_rows.extend(fut.result())
            if done % 5 == 0 or done == len(futs):
                print(f"  {done}/{len(sweep_targets)} clades … "
                      f"({len(all_rows)-seed_count} swept so far)", end="\r", flush=True)
    print()
    print(f"  Sweep assemblies selected: {len(all_rows) - seed_count}")

    # Deduplicate before retaining: one accession can qualify under a seed and again
    # under the sweep. Keep the seed row, which carries the curated taxid and reason.
    seen, deduped, dropped = set(), [], 0
    for r in sorted(all_rows, key=lambda r: 0 if r["source"] == "seed" else 1):
        if r["accession"] in seen:
            dropped += 1
            continue
        seen.add(r["accession"])
        deduped.append(r)
    if dropped:
        print(f"  Dropped {dropped} duplicate accessions (same assembly under a seed "
              f"and a sweep taxid)")
    all_rows = deduped

    if retain:
        all_rows = retain_previous(all_rows, retain)

    all_rows.sort(key=lambda r: (r["kingdom"], r["source"], r["organism_name"].lower(),
                                 int(r["selection_rank"] or 0)))
    return all_rows


def write_candidates(rows: list, version: str = DEFAULT_VERSION,
                     floor: int = DEFAULT_FLOOR) -> None:
    """Write data/ref_candidates_<version>.tsv and point ref_candidates.tsv at it.

    The version-stamped name is what keeps the previous build's selection on disk:
    diffing it against the new one is how the download delta is worked out, and
    kraken_db_busco.py reads the plain ref_candidates.tsv name, so the symlink lets
    it pick up the current version without being passed a path.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    versioned = DATA_DIR / f"ref_candidates_{version}.tsv"
    with open(versioned, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CANDIDATE_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # Before db_v3, ref_candidates.tsv was a REGULAR FILE holding the previous
    # selection, and that selection is the only record of what the current database
    # was built from. Replacing it with a symlink would destroy it, so preserve it
    # first and never overwrite an existing preserved copy.
    if CANDIDATES_TSV.is_file() and not CANDIDATES_TSV.is_symlink():
        legacy = DATA_DIR / "ref_candidates_legacy.tsv"
        if legacy.exists():
            raise SystemExit(
                f"{CANDIDATES_TSV} is a real file and {legacy} already exists — "
                f"refusing to overwrite either. Move or remove one of them by hand.")
        CANDIDATES_TSV.rename(legacy)
        print(f"  preserved the previous selection: {CANDIDATES_TSV.name} -> {legacy.name}")
    elif CANDIDATES_TSV.is_symlink():
        CANDIDATES_TSV.unlink()
    CANDIDATES_TSV.symlink_to(versioned.name)

    n_cds    = sum(1 for r in rows if r["fasta_type"] == "cds")
    n_genome = sum(1 for r in rows if r["fasta_type"] == "genome")
    lineage_counts = defaultdict(int)
    per_taxon = defaultdict(int)
    for r in rows:
        lineage_counts[r["busco_lineage"]] += 1
        per_taxon[r["taxid"]] += 1
    source_counts = defaultdict(int)
    for r in rows:
        source_counts[r["source"]] += 1

    depth_hist = defaultdict(int)
    for n in per_taxon.values():
        depth_hist[min(n, 5)] += 1
    at_floor = sum(v for k, v in depth_hist.items() if floor and k >= floor)

    print(f"\n── Selection summary ────────────────────────────────────────────")
    print(f"  Total candidates:  {len(rows):>6,}")
    print(f"  Distinct taxa:     {len(per_taxon):>6,}")
    print(f"  CDS FASTA:         {n_cds:>6,}")
    print(f"  Genomic FASTA:     {n_genome:>6,}  (no annotation — not downloaded)")
    print(f"  By source:")
    for source, n in sorted(source_counts.items()):
        print(f"    {source:<30} {n:>6,}")
    print(f"  Assemblies per taxon:")
    for depth in sorted(depth_hist):
        label = f"{depth}" if depth < 5 else "5+"
        print(f"    {label:>3}: {depth_hist[depth]:>6,} taxa")
    if floor:
        print(f"  Taxa at the floor of {floor}: {at_floor:,}/{len(per_taxon):,} "
              f"({at_floor/max(len(per_taxon),1)*100:.0f}%) — the rest have no deeper "
              f"assemblies in existence and go in shallow, see the docstring")
    print(f"\nOutput: {versioned}")
    print(f"        {CANDIDATES_TSV} -> {versioned.name}")


# ── CDS download (the ONE place in kraken_db_* that downloads) ───────────────

# Genomic assembly downloads (include="genome", used for HAL host references)
# range up to ~22.4Gb (Pinus radiata) across the field/aerial cohort's named
# host candidates, with wheat (Triticum aestivum, ~14-17Gb, ~55% of the whole
# cohort) the single most common case — nowhere near a corner case. A flat
# 300s timeout here is far too short for those and was the real cause of the
# "300s timeout on datasets download genome" failures seen 2026-08-31 and
# reproduced 2026-09-07 (both on the login node AND, when retested via SLURM
# on an actual compute node, at the exact same accession/timeout — ruling out
# the earlier "login-node network throttling" theory entirely; it's simply
# that large genome downloads need more than 300s, full stop). CDS downloads
# (include="cds", used for pathogen references) are typically megabytes, so
# they complete almost immediately regardless of this value — a generous
# shared timeout costs nothing there.
_DOWNLOAD_TIMEOUT = 3600   # 1hr — comfortably covers even the largest (~22Gb) genome
_UNZIP_TIMEOUT    = 1800   # 30min — large zips can genuinely take a while to extract


def download_cds(accession: str, dest_dir: Path, include: str = "cds") -> list:
    """Download sequence FASTA for accession to dest_dir. Resumable: skips if
    .fna files already present. Returns list of .fna paths (empty on failure).
    include: 'cds' (default, used for pathogen references — smaller, keeps the
    Kraken2 DB k-mer-specific) or 'genome' (whole genomic assembly — used for
    host references by kraken_run_select.py, since BBSplit aligns reads rather
    than doing k-mer LCA, so it doesn't need CDS and genomic sequence also
    catches intron/UTR-spanning reads a CDS-only reference would miss; genomic
    assemblies are also far more available for plants than annotated ones —
    NCBI's plant gene-annotation pipeline coverage is much patchier than for
    fungi/vertebrates, so requiring annotation would exclude most hosts)."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    existing = [f for f in dest_dir.glob("**/*.fna")
                if not f.name.endswith(("_clean.fna", "_combined.fna", ".tagged.fna"))]
    if existing:
        return existing

    zip_path = dest_dir / "ncbi_dataset.zip"
    cmd = ["datasets", "download", "genome", "accession", accession,
           "--include", include, "--filename", str(zip_path)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=_DOWNLOAD_TIMEOUT)
    if r.returncode != 0 or not zip_path.exists():
        print(f"  [{_ts()}] {accession}  download FAILED (rc={r.returncode})", flush=True)
        return []

    subprocess.run(["unzip", "-q", "-o", str(zip_path), "-d", str(dest_dir)],
                   capture_output=True, timeout=_UNZIP_TIMEOUT)
    zip_path.unlink(missing_ok=True)

    fnas = [f for f in dest_dir.glob("**/*.fna")
            if not f.name.endswith(("_clean.fna", "_combined.fna", ".tagged.fna"))]
    if not fnas:
        for gz in dest_dir.glob("**/*.fna.gz"):
            out = gz.with_suffix("")
            with gzip.open(gz, "rb") as fin, open(out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            gz.unlink()
        fnas = [f for f in dest_dir.glob("**/*.fna")
                if not f.name.endswith(("_clean.fna", "_combined.fna", ".tagged.fna"))]
    return fnas


def run_download(candidates: list, genomes_dir: Path, workers: int) -> None:
    cds_rows = [r for r in candidates if r["fasta_type"] == "cds"]
    print(f"\nDownloading CDS for {len(cds_rows):,} candidates → {genomes_dir}")
    genomes_dir.mkdir(parents=True, exist_ok=True)

    n_ok, n_fail, n_cached = 0, 0, 0
    t0 = time.time()

    def work(row):
        acc = row["accession"]
        dest = genomes_dir / acc
        already = list(dest.glob("*.fna")) if dest.exists() else []
        fnas = download_cds(acc, dest)
        return acc, fnas, bool(already)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work, row): row["accession"] for row in cds_rows}
        for done, fut in enumerate(as_completed(futs), 1):
            acc, fnas, was_cached = fut.result()
            if fnas:
                n_ok += 1
                if was_cached:
                    n_cached += 1
            else:
                n_fail += 1
            if done % 25 == 0 or done == len(cds_rows):
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                print(f"  [{_ts()}] [{done}/{len(cds_rows)}] ok={n_ok} "
                      f"(cached={n_cached}) fail={n_fail}  rate={rate:.1f}/s", flush=True)

    print(f"\n── Download summary ─────────────────────────────────────────────")
    print(f"  Downloaded/cached OK: {n_ok:,}  (already on disk: {n_cached:,})")
    print(f"  Failed:               {n_fail:,}")
    print(f"\nNext: python kraken/utilities/busco.py")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scope", action="store_true",
                    help="Reconnaissance report only (pangenome.tsv, genus_fill.tsv); "
                         "does not select or download")
    ap.add_argument("--download", action="store_true",
                    help="After selecting candidates, download CDS for all fasta_type=cds rows")
    ap.add_argument("--select-only", action="store_true",
                    help="Select candidates and write ref_candidates.tsv without downloading "
                         "(default if neither --scope nor --download given)")
    ap.add_argument("--genomes-dir", default=str(DEFAULT_GENOMES_DIR),
                    help="CDS download directory (default: kraken/search/data/cds/pathogen)")
    ap.add_argument("--version", default=DEFAULT_VERSION,
                    help=f"Build version, names the output table (default: {DEFAULT_VERSION}). "
                         "Writes ref_candidates_<version>.tsv so the previous version's "
                         "selection stays on disk to diff the download delta against")
    ap.add_argument("--floor", type=int, default=DEFAULT_FLOOR,
                    help=f"Assemblies per species to aim for (default: {DEFAULT_FLOOR}). "
                         "Applies to swept species and to seeds in undetected genera; "
                         "never caps a detected seed. 0 takes every assembly, which is "
                         "db_v2's behaviour. Species with fewer in existence go in shallow")
    ap.add_argument("--sweep-rank", default=DEFAULT_SWEEP_RANK, choices=SWEEP_RANKS,
                    help=f"Rank to sweep for breadth (default: {DEFAULT_SWEEP_RANK}). "
                         "Ranks above genus are resolved through holobase's taxonomy "
                         "tree; genus needs no tree. Order captures all of Pucciniales, "
                         "where 11 of the 18 rust species missing from db_v2 sit outside "
                         "its rust genera")
    ap.add_argument("--exclude-broad", action="store_true",
                    help="Exclude BROAD_GENERA (Aspergillus, Penicillium …) from the "
                         "sweep, as db_v2 did. Off by default: the sweep exists to give "
                         "LCA competitors and a saprophyte competes as well as a pathogen")
    ap.add_argument("--retain", default=str(DATA_DIR / "ref_candidates_legacy.tsv"),
                    help="Previous candidate table whose assemblies must never be "
                         "dropped, enforcing the no-cap policy. The sweep's "
                         "min(floor, available) is a cap for species it owns, and "
                         "db_v2 took non-seed species with no limit, so without this a "
                         "rebuild trims deep species the bias analysis says to keep. "
                         "Pass --retain '' to select from the rules alone")
    ap.add_argument("--from-table", nargs="?", const=str(CANDIDATES_TSV), default=None,
                    help="Skip selection and download exactly what this candidate table "
                         "names (default: the ref_candidates.tsv symlink). Use this on "
                         "Setonix: selection is reviewed locally and committed, and "
                         "re-selecting there would re-query NCBI, whose catalogue moves, "
                         "so the cluster's table would stop matching the repo's. "
                         "Requires --download")
    ap.add_argument("--relabel-lineage", nargs="?", const=str(CANDIDATES_TSV),
                    default=None,
                    help="Recompute the busco_lineage column of an existing candidate "
                         "table in place and exit (default: the ref_candidates.tsv "
                         "symlink). Use after changing busco_lineage_for(); avoids "
                         "re-selecting, which would re-query NCBI and could change the "
                         "assembly set the CDS were downloaded against")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_db_search.log")
    link_latest(logs_base, log_dir / "kraken_db_search.log")
    sys.stdout = log

    try:
        # --from-table downloads what a table already names, so it needs neither the
        # PHI-base seeds nor holobase. Handle it before loading either: phibase_db.json
        # is gitignored and regenerated by stat/build.py, so on Setonix it is
        # routinely absent after a scratch purge and loading it here failed the job in
        # 15 seconds.
        if args.relabel_lineage:
            # Recompute busco_lineage in place, without re-selecting. The lineage is a
            # label on an already-chosen assembly set, so re-running selection to fix it
            # would re-query NCBI and risk changing the set the CDS were downloaded
            # against. Cheap to re-run whenever busco_lineage_for() changes.
            global PHYLUM_OF
            PHYLUM_OF = holobase_phylum_map()
            table = Path(args.relabel_lineage)
            if not table.exists():
                raise SystemExit(f"--relabel-lineage {table} not found")
            with open(table) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            changed = defaultdict(int)
            for r in rows:
                new = busco_lineage_for(r["kingdom"], r["organism_name"], r["taxid"])
                if new != r["busco_lineage"]:
                    changed[(r["busco_lineage"], new)] += 1
                    r["busco_lineage"] = new
            with open(table, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=CANDIDATE_COLS, delimiter="\t",
                                   extrasaction="ignore")
                w.writeheader()
                w.writerows(rows)
            total = sum(changed.values())
            print(f"Relabelled {total} of {len(rows)} rows in {table}")
            for (old, new), n in sorted(changed.items(), key=lambda kv: -kv[1]):
                print(f"  {n:>5}  {old} -> {new}")
            if not total:
                print("  (nothing to change)")
            return

        if args.from_table:
            if not args.download:
                raise SystemExit("--from-table only makes sense with --download")
            table = Path(args.from_table)
            if not table.exists():
                raise SystemExit(f"--from-table {table} not found")
            with open(table) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            n_cds = sum(1 for r in rows if r["fasta_type"] == "cds")
            print(f"Downloading from {table} as selected: {len(rows):,} candidates, "
                  f"{n_cds:,} with CDS")
            run_download(rows, Path(args.genomes_dir), args.workers)
            return

        db  = json.load(open(DB_PATH))
        t2n = db["taxid_to_name"]
        seeds = {}
        for t in set(db["fungal_to_seed"].values()):
            seeds[int(t)] = "fungal"
        for t in set(db["oomycete_to_seed"].values()):
            seeds[int(t)] = "oomycete"

        if args.scope:
            run_scope(seeds, t2n, args.workers)
            return

        rows = run_select(seeds, t2n, args.workers, args.floor,
                          args.sweep_rank, args.exclude_broad,
                          Path(args.retain) if args.retain else None)
        write_candidates(rows, args.version, args.floor)

        if args.download:
            run_download(rows, Path(args.genomes_dir), args.workers)
        else:
            print(f"\nNext: python kraken/search.py --download   "
                  f"(or re-run this with --download)")
    finally:
        log.close()


if __name__ == "__main__":
    main()
