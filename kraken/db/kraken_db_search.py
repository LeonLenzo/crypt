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
              get the single best-quality assembly only. Falls back through
              scaffold-plus, then contig-level — always take something over nothing.

              GENUS FILL (breadth): for PHI-base genera that appear in STAT
              detections, add ALL annotated scaffold-plus assemblies for each
              non-PHI-base species within that genus. Broad/saprophytic genera
              (Aspergillus, Penicillium …) are excluded.

              Output: data/ref_candidates.tsv (one row per candidate assembly)

  --download  Download CDS FASTA for every candidate with fasta_type=cds into
              data/cds/pathogen/{accession}/. This is the ONLY place in the kraken_db_*
              pipeline that downloads CDS — kraken_db_busco.py and
              kraken_db_build.py both read from here, never fetch their own copies.
              Resumable: skips accessions with .fna files already on disk.

Typical run (select + download in one pass):
    python kraken/db/kraken_db_search.py --download

Then: python kraken/db/kraken_db_busco.py   (on Setonix)

Output:
    kraken/output/db/search/data/ref_candidates.tsv   (tracked)
    kraken/output/db/search/data/pangenome.tsv        (tracked, --scope only)
    kraken/output/db/search/data/genus_fill.tsv       (tracked, --scope only)
    kraken/output/db/search/data/cds/pathogen/{accession}/  (gitignored, --download only)
"""

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _util import _Tee, make_log_dir, link_latest

DB_PATH      = Path("stat/output/stat_build/data/phibase_db.json")
RUNS_TSV     = Path("stat/output/stat_filter/data/runs.tsv")
OUT_DIR      = Path("kraken/output/db/search")
DATA_DIR     = OUT_DIR / "data"
CANDIDATES_TSV = DATA_DIR / "ref_candidates.tsv"
DEFAULT_GENOMES_DIR = DATA_DIR / "cds" / "pathogen"

YEAR_BIN_SIZE = 5
LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}

# Genera excluded from fill-in: broad/saprophytic, not primarily plant pathogens
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

def busco_lineage_for(kingdom: str, organism_name: str) -> str:
    if kingdom == "oomycete":
        return "stramenopiles_odb10"
    genus = organism_name.split()[0]
    if genus in BASIDIOMYCETE_GENERA:
        return "basidiomycota_odb10"
    if genus in CHYTRID_GENERA:
        return "fungi_odb10"
    return "ascomycota_odb10"

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


def select_seed(taxid: int, name: str, kingdom: str,
                detected: bool, assemblies: list) -> list:
    """detected=True → all scaffold-plus+annotated assemblies, diversity-ordered.
    detected=False → single best-quality assembly only. Falls back through
    annotation tiers if needed."""
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
    pairs = [(pool[0], "best_quality")] if not detected else greedy_ordered(pool)
    return [to_row(a, taxid, name, kingdom, "seed", rank + 1, reason)
            for rank, (a, reason) in enumerate(pairs)]


def select_genus_fill(kingdom: str, covered_taxids: set, assemblies: list) -> list:
    """For a detected genus: take ALL annotated scaffold-plus assemblies for each
    species not already covered by a PHI-base seed, diversity-ordered within species."""
    by_species = defaultdict(list)
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
        for rank, (a, reason) in enumerate(greedy_ordered(pool)):
            rows.append(to_row(a, tid, org_name, kingdom, "genus_fill", rank + 1, reason))
    return rows


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

def run_select(seeds: dict, t2n: dict, workers: int) -> list:
    covered = set(seeds.keys())
    detected_genera = load_detected_genera()
    print(f"STAT-detected genera: {len(detected_genera)}")

    print(f"\n[1/2] Querying {len(seeds)} seed taxids (all assemblies for detected genera) …")
    all_rows = []

    def work_seed(taxid, kingdom):
        name = t2n.get(str(taxid), str(taxid))
        detected = name.split()[0] in detected_genera
        return select_seed(taxid, name, kingdom, detected, datasets_query(taxid))

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

    fill_genera = {}
    for tid, kingdom in seeds.items():
        parts = t2n.get(str(tid), "").split()
        if not parts:
            continue
        genus = parts[0]
        if genus in detected_genera and genus not in BROAD_GENERA:
            fill_genera[genus] = kingdom

    print(f"\n[2/2] Querying {len(fill_genera)} genera for fill-in (all qualifying assemblies) …")

    def work_genus(genus, kingdom):
        return select_genus_fill(kingdom, covered, datasets_query(genus))

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work_genus, g, k): g for g, k in fill_genera.items()}
        for done, fut in enumerate(as_completed(futs), 1):
            all_rows.extend(fut.result())
            if done % 5 == 0 or done == len(futs):
                print(f"  {done}/{len(fill_genera)} genera … ({len(all_rows)-seed_count} fill-in so far)",
                      end="\r", flush=True)
    print()
    print(f"  Genus fill-in assemblies selected: {len(all_rows) - seed_count}")

    all_rows.sort(key=lambda r: (r["kingdom"], r["source"], r["organism_name"].lower(), r["selection_rank"]))
    return all_rows


def write_candidates(rows: list) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(CANDIDATES_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CANDIDATE_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n_cds    = sum(1 for r in rows if r["fasta_type"] == "cds")
    n_genome = sum(1 for r in rows if r["fasta_type"] == "genome")
    lineage_counts = defaultdict(int)
    for r in rows:
        lineage_counts[r["busco_lineage"]] += 1

    print(f"\n── Selection summary ────────────────────────────────────────────")
    print(f"  Total candidates:  {len(rows):>6,}")
    print(f"  CDS FASTA:         {n_cds:>6,}")
    print(f"  Genomic FASTA:     {n_genome:>6,}  (no annotation — not downloaded)")
    print(f"  BUSCO lineages:")
    for lineage, n in sorted(lineage_counts.items()):
        print(f"    {lineage:<30} {n:>6,}")
    print(f"\nOutput: {CANDIDATES_TSV}")


# ── CDS download (the ONE place in kraken_db_* that downloads) ───────────────

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
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0 or not zip_path.exists():
        print(f"  [{_ts()}] {accession}  download FAILED (rc={r.returncode})", flush=True)
        return []

    subprocess.run(["unzip", "-q", "-o", str(zip_path), "-d", str(dest_dir)],
                   capture_output=True, timeout=300)
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
    print(f"\nNext: python kraken/db/kraken_db_busco.py")


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
                    help="CDS download directory (default: kraken/output/db/search/data/cds/pathogen)")
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

        rows = run_select(seeds, t2n, args.workers)
        write_candidates(rows)

        if args.download:
            run_download(rows, Path(args.genomes_dir), args.workers)
        else:
            print(f"\nNext: python kraken/db/kraken_db_search.py --download   "
                  f"(or re-run this with --download)")
    finally:
        log.close()


if __name__ == "__main__":
    main()
