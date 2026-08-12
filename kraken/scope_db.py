#!/usr/bin/env python3
"""
scope_db.py — scope the Kraken2 DB expansion before rebuilding.

Two output tables:

  kraken/output/scope/pangenome.tsv
    Per PHI-base seed: how many assemblies exist at NCBI, geographic/temporal
    diversity — informs how many assemblies to include per species.

  kraken/output/scope/genus_fill.tsv
    Per PHI-base genus: assemblies and species NOT already in PHI-base seeds —
    informs how much the genus fill-in adds.

Run from crypt/:
    python kraken/scope_db.py [--workers N]
"""

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

DB_PATH  = Path("stat/output/build/data/phibase_db.json")
OUT_DIR  = Path("kraken/output/scope")

KINGDOMS = ["fungal", "oomycete", "nematode"]
LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}


# ── helpers ───────────────────────────────────────────────────────────────────

def datasets_query(taxon: str | int) -> list[dict]:
    """Return list of assembly summary dicts for a taxon (name or taxid)."""
    cmd = [
        "datasets", "summary", "genome", "taxon", str(taxon),
        "--as-json-lines", "--limit", "all",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        return []
    rows = []
    for line in r.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    return rows


def get_geo(assembly: dict) -> str:
    """Extract geo_loc_name from biosample attributes."""
    biosample = (
        assembly.get("assembly_info", {}).get("biosample") or {}
    )
    for attr in biosample.get("attributes", []) or []:
        if attr.get("name") == "geo_loc_name":
            v = attr.get("value", "").strip()
            return v.split(":")[0].strip() if v and v != "not applicable" else ""
    return ""


def assembly_level(assembly: dict) -> str:
    return assembly.get("assembly_info", {}).get("assembly_level", "")


def release_date(assembly: dict) -> str:
    return assembly.get("assembly_info", {}).get("release_date", "") or ""


def has_annotation(assembly: dict) -> bool:
    return bool(assembly.get("annotation_info"))


def scaffold_n50_kb(assembly: dict) -> float:
    v = assembly.get("assembly_stats", {}).get("scaffold_n50")
    try:
        return round(int(v) / 1000, 1)
    except (TypeError, ValueError):
        return 0.0


def level_ge_scaffold(lvl: str) -> bool:
    return LEVEL_RANK.get(lvl, 0) >= 2


def summarise_assemblies(rows: list[dict]) -> dict:
    """Aggregate list of NCBI assembly dicts into summary stats."""
    if not rows:
        return {
            "n_total": 0, "n_scaffold_plus": 0, "n_annotated": 0,
            "date_min": "", "date_max": "", "n_countries": 0,
            "countries": "", "best_level": "",
        }
    dates       = sorted(d for a in rows if (d := release_date(a)))
    countries   = sorted({c for a in rows if (c := get_geo(a))})
    scaffold_up = [a for a in rows if level_ge_scaffold(assembly_level(a))]
    annotated   = [a for a in rows if has_annotation(a)]
    levels      = [assembly_level(a) for a in rows]
    best        = max(levels, key=lambda l: LEVEL_RANK.get(l, 0), default="")
    return {
        "n_total":       len(rows),
        "n_scaffold_plus": len(scaffold_up),
        "n_annotated":   len(annotated),
        "date_min":      dates[0]  if dates else "",
        "date_max":      dates[-1] if dates else "",
        "n_countries":   len(countries),
        "countries":     "; ".join(countries[:20]),
        "best_level":    best,
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    db     = json.load(open(DB_PATH))
    t2n    = db["taxid_to_name"]
    seeds: dict[int, str] = {}   # taxid → kingdom

    for kingdom in KINGDOMS:
        key = f"{kingdom}_to_seed"
        if key not in db:
            continue
        for seed_tid in set(db[key].values()):
            seeds[seed_tid] = kingdom

    print(f"Seeds: {len(seeds)} ({', '.join(f'{k}: {sum(1 for v in seeds.values() if v==k)}' for k in KINGDOMS)})")

    # ── Part 1: pan-genome — assemblies per seed ──────────────────────────────

    print("\n[1/2] Querying assemblies per seed taxid …")
    pan_rows = []

    def query_seed(tid: int, kingdom: str) -> dict:
        name = t2n.get(str(tid), str(tid))
        rows = datasets_query(tid)
        s    = summarise_assemblies(rows)
        return {"taxid": tid, "name": name, "kingdom": kingdom, **s}

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(query_seed, tid, kgd): tid for tid, kgd in seeds.items()}
        done = 0
        for fut in as_completed(futs):
            done += 1
            row = fut.result()
            pan_rows.append(row)
            if done % 10 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} seeds …", end="\r", flush=True)

    print()
    pan_rows.sort(key=lambda r: (-r["n_total"], r["name"]))

    pan_path = OUT_DIR / "pangenome.tsv"
    pan_cols = [
        "taxid", "name", "kingdom",
        "n_total", "n_scaffold_plus", "n_annotated",
        "best_level", "date_min", "date_max",
        "n_countries", "countries",
    ]
    with open(pan_path, "w") as fh:
        fh.write("\t".join(pan_cols) + "\n")
        for r in pan_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in pan_cols) + "\n")
    print(f"Written: {pan_path}")

    # Quick stats
    has_assemblies = [r for r in pan_rows if r["n_total"] > 0]
    print(f"\nPan-genome summary:")
    print(f"  Seeds with ≥1 assembly: {len(has_assemblies)}/{len(pan_rows)}")
    print(f"  Seeds with ≥5 scaffold-plus assemblies: {sum(1 for r in pan_rows if r['n_scaffold_plus'] >= 5)}")
    print(f"  Seeds with ≥2 countries: {sum(1 for r in pan_rows if r['n_countries'] >= 2)}")
    total_assemblies = sum(r["n_total"] for r in pan_rows)
    total_scaffold   = sum(r["n_scaffold_plus"] for r in pan_rows)
    print(f"  Total assemblies available: {total_assemblies:,} ({total_scaffold:,} scaffold-plus)")

    # ── Part 2: genus fill-in — new species within PHI-base genera ────────────

    print("\n[2/2] Querying genus-level fill-in …")

    # Map genus → list of seed taxids in that genus
    seed_taxids_set = set(seeds.keys())
    genus_to_seeds: dict[str, list[int]] = defaultdict(list)
    for tid in seeds:
        name = t2n.get(str(tid), "")
        parts = name.split()
        if len(parts) >= 2:
            genus_to_seeds[parts[0]].append(tid)

    print(f"  Unique genera: {len(genus_to_seeds)}")

    fill_rows = []

    def query_genus(genus: str, seed_tids: list[int]) -> dict:
        rows = datasets_query(genus)
        if not rows:
            return None

        # Separate PHI-base vs new species
        new_species: dict[int, list[dict]] = defaultdict(list)
        phibase_n = 0
        for a in rows:
            tid = a.get("organism", {}).get("tax_id")
            if tid in seed_taxids_set:
                phibase_n += 1
            else:
                if tid:
                    new_species[tid].append(a)

        new_tids = set(new_species.keys())
        new_assemblies = sum(len(v) for v in new_species.values())
        new_scaffold   = sum(
            1 for v in new_species.values() for a in v
            if level_ge_scaffold(assembly_level(a))
        )
        new_annotated  = sum(
            1 for v in new_species.values() for a in v
            if has_annotation(a)
        )
        new_sp_names = sorted(
            {a.get("organism", {}).get("organism_name", "") for a in rows
             if a.get("organism", {}).get("tax_id") not in seed_taxids_set
             and a.get("organism", {}).get("tax_id")}
        )
        return {
            "genus":             genus,
            "n_phibase_seeds":   len(seed_tids),
            "n_phibase_assemblies": phibase_n,
            "n_new_species":     len(new_tids),
            "n_new_assemblies":  new_assemblies,
            "n_new_scaffold":    new_scaffold,
            "n_new_annotated":   new_annotated,
            "new_species":       "; ".join(new_sp_names[:30]),
        }

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {
            pool.submit(query_genus, genus, stids): genus
            for genus, stids in genus_to_seeds.items()
        }
        done = 0
        for fut in as_completed(futs):
            done += 1
            row = fut.result()
            if row:
                fill_rows.append(row)
            if done % 5 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} genera …", end="\r", flush=True)

    print()
    fill_rows.sort(key=lambda r: -r["n_new_assemblies"])

    fill_path = OUT_DIR / "genus_fill.tsv"
    fill_cols = [
        "genus",
        "n_phibase_seeds", "n_phibase_assemblies",
        "n_new_species", "n_new_assemblies", "n_new_scaffold", "n_new_annotated",
        "new_species",
    ]
    with open(fill_path, "w") as fh:
        fh.write("\t".join(fill_cols) + "\n")
        for r in fill_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in fill_cols) + "\n")
    print(f"Written: {fill_path}")

    total_new_sp  = sum(r["n_new_species"]    for r in fill_rows)
    total_new_asm = sum(r["n_new_assemblies"] for r in fill_rows)
    total_new_scf = sum(r["n_new_scaffold"]   for r in fill_rows)
    print(f"\nGenus fill-in summary:")
    print(f"  New species (non-PHI-base, within PHI-base genera): {total_new_sp:,}")
    print(f"  New assemblies total: {total_new_asm:,}")
    print(f"  New assemblies scaffold-plus: {total_new_scf:,}")

    print("\nTop genera by new assemblies:")
    for r in fill_rows[:15]:
        print(f"  {r['genus']:25s}  +{r['n_new_species']:3d} spp  +{r['n_new_assemblies']:5d} asm  ({r['n_new_scaffold']} scaffold+)")


if __name__ == "__main__":
    main()
