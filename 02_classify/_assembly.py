#!/usr/bin/env python3
"""_assembly.py — the NCBI assembly vocabulary, carried across from 03_kraken/01_search.

This is METHOD, not wiring, and it is lifted verbatim from `kraken_search.py` rather than
rewritten. Every function here encodes something that was measured or that cost a bug to find,
and the comments are the record of that:

  * `datasets_query` warns on a non-zero exit instead of returning []. A taxonomic HOMONYM
    makes the CLI refuse ("exact match for more than one taxid") and the clade then selects
    nothing with no sign anything went wrong. *Cryptococcus* is both a fungus (5206) and a
    scale insect (79213); *Rhizophagus* a fungus (1129544) and a beetle (295919). Between them
    that silently cost 42 assemblies and 17 species.
  * `greedy_ordered` orders by geographic then temporal diversity, not by quality alone.
  * `busco_lineage_for` resolves phylum from holobase FIRST and falls back to the genus sets.
    The genus lists alone covered db_v2 because db_v2 only ever saw pathogens; 160 of db_v3's
    228 Basidiomycota candidates fall in 63 genera outside them and would be scored against
    ascomycota_odb10, giving spuriously low completeness for exactly the breadth that was
    added on purpose.
  * `holobase_species_map` normalises every NCBI taxid to its species node. Without it a
    seed's own strains are not recognised as covered and the depth floor counts per strain,
    so a species split across several nodes passes the floor while staying shallow. That was
    87 species, worst *Zymoseptoria tritici* across 6 nodes.
  * `download_cds` carries the 3600s timeout. The old flat 300s was the real cause of the
    "300s timeout on datasets download genome" failures of 2026-08-31, reproduced on a compute
    node, which ruled out the login-node throttling theory: large genomes simply need longer.

What is deliberately NOT here, because it is policy rather than mechanics, and policy changed
for db_v4: the PHI-base seed source, `load_detected_genera()` (a STAT dependency), and
`retain_previous()`. See `01_refs/select_refs.py`.

`init_taxonomy()` fills the two module globals that `to_row` and the selector read, and returns
the species map. The original carried these as `global` statements inside `run_select`, which
meant importing the module was not enough to make `to_row` correct.
"""

import gzip
import json
import shutil
import sqlite3
import subprocess
import sys
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _layout import HOLOBASE  # noqa: E402

YEAR_BIN_SIZE = 5
LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

# Basidiomycete genera -> basidiomycota_odb10. A FALLBACK only; holobase's phylum is tried
# first. Kept because holobase may be absent and because it is the db_v2-era record.
BASIDIOMYCETE_GENERA = {
    "Puccinia", "Melampsora", "Phakopsora", "Hemileia", "Uromyces",
    "Tranzschelia", "Phragmidium", "Gymnosporangium",        # rusts
    "Ustilago", "Tilletia", "Sporisorium", "Testicularia",   # smuts
    "Mycosarcoma",                                            # Ustilaginomycotina
    "Rhizoctonia", "Heterobasidion", "Moniliophthora", "Crinipellis",
}

# Chytrid/other early-diverging fungi with no phylum-level BUSCO lineage -> fungi_odb10
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

# Filled by datasets_query when a clade returns nothing, so a caller can report "resolved but
# NCBI holds no genome" separately from "the query itself failed". Two different facts: the
# first is data, the second is a bug waiting to be noticed.
NO_DATA = {}

def datasets_query(taxon) -> list:
    r = subprocess.run(
        ["datasets", "summary", "genome", "taxon", str(taxon),
         "--as-json-lines", "--limit", "all"],
        capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        # Do NOT swallow this. A genus that is a taxonomic HOMONYM makes the CLI exit
        # non-zero ("exact match for more than one taxid") and the clade then selects zero
        # assemblies with no sign anything went wrong. Found 2026-10-13: Cryptococcus is
        # both a fungus (5206) and a scale insect (79213); Rhizophagus is both a fungus
        # (1129544) and a beetle (295919). Pass the taxid instead of the name for those.
        err = (r.stderr or "").strip().splitlines()
        msg = err[0] if err else "exit " + str(r.returncode)
        # "valid for X, but no genome data is currently available" is an ANSWER, not a
        # failure: the taxon resolved and NCBI holds nothing for it. Surveying thousands of
        # clades turns that into hundreds of WARNING lines, which trains the reader to skim
        # past the one line that means a homonym just silently returned zero. Counted
        # instead, by the caller, from the empty return.
        if "no genome data is currently available" not in msg:
            print(f"\n  WARNING datasets_query({taxon}) failed: {msg}", flush=True)
            NO_DATA.setdefault("_failures", []).append((str(taxon), msg))
        else:
            NO_DATA.setdefault("_no_genome_data", []).append(str(taxon))
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
    print(f"\nNext: python 02_classify/utilities/busco.py")


# ── taxonomy, read once ───────────────────────────────────────────────────────

# Filled by init_taxonomy(). The selector names a species from here rather than from the
# assembly's own organism_name, because a strain name carries an isolate suffix that makes
# one species look like several.
SPECIES_NAME = {}


def init_taxonomy(db_path: Path = HOLOBASE) -> dict:
    """Read holobase once and fill the module globals to_row/busco_lineage_for depend on.

    Returns {ncbi taxid: species-level ncbi taxid}, with species mapping to themselves so a
    caller may apply `to_species.get(tid, tid)` unconditionally.

    This exists because the original carried PHYLUM_OF and SPECIES_NAME as `global`
    statements inside `run_select`: importing the module was not enough to make `to_row`
    assign the right BUSCO lineage, and any new entry point had to know to call the selector
    first. Both maps are a full pass over a ~3M-row table, so they are read once here.
    """
    global PHYLUM_OF, SPECIES_NAME
    to_species, SPECIES_NAME = holobase_species_map(db_path)
    PHYLUM_OF = holobase_phylum_map(db_path)
    return to_species
