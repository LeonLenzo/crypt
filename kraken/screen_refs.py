#!/usr/bin/env python3
"""
screen_refs.py — select candidate assemblies for BUSCO quality screening.

Two selection passes:

  SEED (pan-genome): for each PHI-base euk seed taxid, collect ALL
  scaffold-plus+annotated assemblies.  Detected-genera seeds get all of them
  ordered by geographic/temporal diversity; undetected-genera seeds get the
  single best-quality assembly only.  Falls back through scaffold-plus, then
  contig-level — always take something rather than nothing.

  GENUS FILL (breadth): for PHI-base genera that appear in STAT detections,
  add ALL annotated scaffold-plus assemblies for each non-PHI-base species
  within that genus.  Broad/saprophytic genera (Aspergillus, Penicillium …)
  are excluded.

Output: kraken/ref_candidates.tsv  (one row per candidate assembly)
The busco_lineage column tells busco_screen.py which BUSCO lineage DB to use.
Run busco_screen.py next, then filter_refs.py, to produce the final ref_screen.tsv.

Run from crypt/:
    python kraken/screen_refs.py [--workers N]
"""

import argparse
import csv
import json
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

DB_PATH   = Path("stat/output/build/data/phibase_db.json")
RUNS_TSV  = Path("stat/output/filter_runs/data/runs.tsv")
OUT_TSV   = Path("kraken/ref_candidates.tsv")

YEAR_BIN_SIZE = 5

LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

# Genera excluded from fill-in: broad/saprophytic, not primarily plant pathogens
BROAD_GENERA = {
    "Aspergillus", "Penicillium", "Trichoderma", "Beauveria", "Metarhizium",
    "Claviceps", "Epichloe", "Ceratocystis", "Leptographium", "Ciboria",
}

# Basidiomycete PHI-base genera → use basidiomycota_odb10 instead of fungi_odb10
BASIDIOMYCETE_GENERA = {
    "Puccinia", "Melampsora", "Phakopsora", "Hemileia", "Uromyces",
    "Tranzschelia", "Phragmidium", "Gymnosporangium",   # rusts
    "Ustilago", "Tilletia", "Sporisorium", "Testicularia",  # smuts
    "Rhizoctonia", "Moniliophthora", "Crinipellis",
}

OUT_COLS = [
    "taxid", "organism_name", "kingdom", "source",
    "accession", "assembly_level", "release_date", "has_annotation",
    "country", "protein_coding_genes", "scaffold_n50_kb", "total_length_mb",
    "fasta_type", "busco_lineage", "selection_rank", "selection_reason",
]


# ── assembly field helpers ────────────────────────────────────────────────────

def _ai(a: dict) -> dict:
    return a.get("assembly_info") or {}

def _ann(a: dict) -> dict:
    return a.get("annotation_info") or {}

def _stats(a: dict) -> dict:
    return a.get("assembly_stats") or {}

def assembly_level(a: dict) -> str:
    return _ai(a).get("assembly_level", "")

def release_date(a: dict) -> str:
    return _ai(a).get("release_date", "") or ""

def has_annotation(a: dict) -> bool:
    return bool(_ann(a))

def get_country(a: dict) -> str:
    for attr in (_ai(a).get("biosample") or {}).get("attributes", []) or []:
        if attr.get("name") == "geo_loc_name":
            v = (attr.get("value") or "").strip()
            if v and v.lower() not in ("not applicable", "missing", "not collected", ""):
                return v.split(":")[0].strip()
    return ""

def get_year_bin(a: dict) -> str:
    y = release_date(a)[:4]
    try:
        return str((int(y) // YEAR_BIN_SIZE) * YEAR_BIN_SIZE)
    except ValueError:
        return ""

def level_rank(a: dict) -> int:
    return LEVEL_RANK.get(assembly_level(a), 0)

def quality_key(a: dict) -> tuple:
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

def scaffold_plus(a: dict) -> bool:
    return level_rank(a) >= 2

def busco_lineage_for(kingdom: str, organism_name: str) -> str:
    if kingdom == "oomycete":
        return "stramenopiles_odb10"
    genus = organism_name.split()[0]
    if genus in BASIDIOMYCETE_GENERA:
        return "basidiomycota_odb10"
    return "fungi_odb10"

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
        "busco_lineage":       busco_lineage_for(kingdom, organism_name),
        "selection_rank":      rank,
        "selection_reason":    reason,
    }


# ── NCBI datasets query ───────────────────────────────────────────────────────

def datasets_query(taxon: str | int) -> list[dict]:
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


# ── diversity-ordered selection (no cap) ─────────────────────────────────────

def greedy_ordered(assemblies: list[dict]) -> list[tuple[dict, str]]:
    """
    Order all assemblies by geographic + temporal diversity.
    assemblies must already be sorted by quality_key descending.
    Returns list of (assembly, reason_str) — all assemblies, no cap.
    """
    if not assemblies:
        return []

    remaining = list(assemblies)
    selected: list[dict] = []
    reasons:  list[str]  = []

    selected.append(remaining.pop(0))
    reasons.append("best_quality")

    while remaining:
        sel_countries = {get_country(a) for a in selected}
        sel_bins      = {get_year_bin(a) for a in selected}

        best_score = -1
        best_idx   = 0
        best_parts: list[str] = []

        for i, a in enumerate(remaining):
            country = get_country(a)
            ybin    = get_year_bin(a)
            score   = 0
            parts:  list[str] = []
            if country and country not in sel_countries:
                score += 2
                parts.append(f"new_country:{country}")
            if ybin and ybin not in sel_bins:
                score += 1
                parts.append(f"new_year_bin:{ybin}")
            if score > best_score:
                best_score = score
                best_idx   = i
                best_parts = parts

        selected.append(remaining.pop(best_idx))
        reasons.append("; ".join(best_parts) if best_parts else "quality_fill")

    return list(zip(selected, reasons))


def select_seed(taxid: int, name: str, kingdom: str,
                detected: bool, assemblies: list[dict]) -> list[dict]:
    """
    Select candidate assemblies for a single seed taxid.
    detected=True  → all scaffold-plus+annotated assemblies, diversity-ordered.
    detected=False → single best-quality assembly only.
    Falls back through annotation tiers if needed.
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

    if not detected:
        pairs = [(pool[0], "best_quality")]
    else:
        pairs = greedy_ordered(pool)

    return [
        to_row(a, taxid, name, kingdom, "seed", rank + 1, reason)
        for rank, (a, reason) in enumerate(pairs)
    ]


# ── genus fill-in ─────────────────────────────────────────────────────────────

def select_genus_fill(kingdom: str,
                      covered_taxids: set[int],
                      assemblies: list[dict]) -> list[dict]:
    """
    For a detected genus: take ALL annotated scaffold-plus assemblies for
    each species not already covered by a PHI-base seed, diversity-ordered
    within each species.
    """
    by_species: dict[int, list[dict]] = defaultdict(list)
    for a in assemblies:
        tid = (a.get("organism") or {}).get("tax_id")
        if not tid or tid in covered_taxids:
            continue
        if scaffold_plus(a) and has_annotation(a):
            by_species[tid].append(a)

    rows = []
    for tid, pool in by_species.items():
        pool.sort(key=quality_key, reverse=True)
        org_name = (pool[0].get("organism") or {}).get("organism_name", str(tid))
        pairs = greedy_ordered(pool)
        for rank, (a, reason) in enumerate(pairs):
            rows.append(to_row(a, tid, org_name, kingdom, "genus_fill", rank + 1, reason))
    return rows


# ── main ──────────────────────────────────────────────────────────────────────

def load_detected_genera() -> set[str]:
    genera: set[str] = set()
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    db  = json.load(open(DB_PATH))
    t2n = db["taxid_to_name"]

    covered: set[int] = set()
    for key in ("fungal_to_seed", "oomycete_to_seed", "nematode_to_seed"):
        covered.update(int(t) for t in db.get(key, {}))

    seeds: dict[int, str] = {}
    for t in set(db["fungal_to_seed"].values()):
        seeds[int(t)] = "fungal"
    for t in set(db["oomycete_to_seed"].values()):
        seeds[int(t)] = "oomycete"

    detected_genera = load_detected_genera()
    print(f"STAT-detected genera: {len(detected_genera)}")

    # ── Pass 1: seeds ─────────────────────────────────────────────────────────
    print(f"\n[1/2] Querying {len(seeds)} seed taxids (all assemblies for detected genera) …")

    all_rows: list[dict] = []

    def work_seed(taxid: int, kingdom: str) -> list[dict]:
        name = t2n.get(str(taxid), str(taxid))
        genus = name.split()[0]
        detected = genus in detected_genera
        assemblies = datasets_query(taxid)
        return select_seed(taxid, name, kingdom, detected, assemblies)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(work_seed, tid, kgd): tid for tid, kgd in seeds.items()}
        done = 0
        for fut in as_completed(futs):
            done += 1
            rows = fut.result()
            all_rows.extend(rows)
            if done % 10 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} seeds … ({len(all_rows)} assemblies so far)",
                      end="\r", flush=True)
    print()

    seed_count = len(all_rows)
    print(f"  Seed assemblies selected: {seed_count}")

    # ── Pass 2: genus fill-in ─────────────────────────────────────────────────
    fill_genera: dict[str, str] = {}
    for tid, kingdom in seeds.items():
        parts = t2n.get(str(tid), "").split()
        if not parts:
            continue
        genus = parts[0]
        if genus and genus in detected_genera and genus not in BROAD_GENERA:
            fill_genera[genus] = kingdom

    print(f"\n[2/2] Querying {len(fill_genera)} genera for fill-in (all qualifying assemblies) …")

    def work_genus(genus: str, kingdom: str) -> list[dict]:
        assemblies = datasets_query(genus)
        return select_genus_fill(kingdom, covered, assemblies)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(work_genus, g, k): g for g, k in fill_genera.items()}
        done = 0
        for fut in as_completed(futs):
            done += 1
            rows = fut.result()
            all_rows.extend(rows)
            if done % 5 == 0 or done == len(futs):
                fill_n = len(all_rows) - seed_count
                print(f"  {done}/{len(fill_genera)} genera … ({fill_n} fill-in assemblies so far)",
                      end="\r", flush=True)
    print()

    fill_count = len(all_rows) - seed_count
    print(f"  Genus fill-in assemblies selected: {fill_count}")

    # ── Write output ──────────────────────────────────────────────────────────
    all_rows.sort(key=lambda r: (r["kingdom"], r["source"], r["organism_name"].lower(), r["selection_rank"]))

    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_cds    = sum(1 for r in all_rows if r["fasta_type"] == "cds")
    n_genome = sum(1 for r in all_rows if r["fasta_type"] == "genome")
    n_seeds_multi = len({r["taxid"] for r in all_rows if r["source"] == "seed" and r["selection_rank"] > 1})
    n_fill_genera = len({r["organism_name"].split()[0] for r in all_rows if r["source"] == "genus_fill"})
    lineage_counts: dict[str, int] = defaultdict(int)
    for r in all_rows:
        lineage_counts[r["busco_lineage"]] += 1

    print(f"\n── Summary ──────────────────────────────────────────────────")
    print(f"  Total candidates:             {len(all_rows):>6,}")
    print(f"    Seed (pan-genome):          {seed_count:>6,}")
    print(f"    Genus fill-in:              {fill_count:>6,}")
    print(f"  CDS FASTA (annotated):        {n_cds:>6,}")
    print(f"  Genomic FASTA (unannotated):  {n_genome:>6,}")
    print(f"  Seeds with >1 assembly:       {n_seeds_multi:>6,}")
    print(f"  Genera contributing fill-in:  {n_fill_genera:>6,}")
    print(f"  BUSCO lineages:")
    for lineage, n in sorted(lineage_counts.items()):
        print(f"    {lineage:<30} {n:>6,}")
    print(f"\nOutput: {OUT_TSV}")
    print(f"Next:   python kraken/busco_screen.py  (on Setonix)")


if __name__ == "__main__":
    main()
