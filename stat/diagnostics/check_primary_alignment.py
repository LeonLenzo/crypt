#!/usr/bin/env python3
"""
scripts/check_primary_alignment.py — Three-way alignment: SRA library vs STAT vs BioProject metadata.

Four output columns (per biosample_representative run):
  sra_organism      library_organism from runs.tsv (declared library organism)
  stat_top_pathogen Top PHI-base/ICTV pathogen from stat_pathogens (first entry, highest %)
  meta_pathogen     First PHI-base pathogen name/genus found in BioProject title+description+abstract
  mismatch_flag     all_match | stat_conflict | meta_conflict | both_conflict | stat_undetected

MAL interpretation:
  sra_organism = declared pathogen (library is the pathogen)
  stat_top_pathogen = top STAT detection — may differ if STAT sees a different organism more clearly
  meta_pathogen = what the BioProject says the study was about

HAL interpretation:
  sra_organism = host plant (library is the host)
  stat_top_pathogen = top STAT-detected pathogen (STAT-derived)
  meta_pathogen = pathogen mentioned in BioProject text
  mismatch_flag compares stat_top_pathogen vs meta_pathogen (sra_organism is host, not comparable)

Usage:
  python scripts/check_primary_alignment.py

Output:
  output/analysis/primary_alignment.tsv
"""

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

RUNS_TSV   = Path("output/02_filter_runs/data/runs.tsv")
PHIBASE_DB = Path("output/00_build/data/phibase_db.json")
KW_TSV     = Path("output/04_filter_kw/data/biosample_kw.tsv")
OUT_DIR    = Path("output/analysis")
OUT_TSV    = OUT_DIR / "primary_alignment.tsv"

OUTPUT_FIELDS = [
    "mode", "Run", "BioSample", "BioProject", "biosample_rep",
    "sra_organism",       # library_organism from runs.tsv
    "stat_top_pathogen",  # top PHI-base pathogen in STAT by %
    "stat_top_pct",       # its % in STAT
    "meta_pathogen",      # first PHI-base match in BioProject title+description+abstract
    "mismatch_flag",      # all_match | stat_conflict | meta_conflict | both_conflict | stat_undetected
    "bp_title",
]


def _genus(name: str) -> str:
    return name.split()[0].lower() if name else ""


_COLON_PCT_RE = re.compile(r':([\d.]+)%$')

def _top_stat(stat_pathogens: str) -> tuple[str, str]:
    """Parse top (first) pathogen name and pct from stat_pathogens column."""
    if not stat_pathogens:
        return "", ""
    first = stat_pathogens.split("; ")[0].strip()
    m = _COLON_PCT_RE.search(first)
    pct  = m.group(1) + "%" if m else ""
    name = _COLON_PCT_RE.sub("", first).strip()
    return name, pct


def _build_meta_search_index(db: dict) -> tuple[set[str], list[str]]:
    """
    Build (pathogen_taxids, sorted_pathogen_names) for text search.
    Returns taxids of all pathogen kingdom entries and a list of pathogen
    binomial names sorted longest-first (so species beat genera on match).
    """
    pathogen_taxids: set[int] = set()
    for dk in ("fungal_to_seed", "bacterial_to_seed", "oomycete_to_seed",
               "nematode_to_seed", "virus_to_seed"):
        for t in db.get(dk, {}):
            pathogen_taxids.add(int(t))

    # Collect names that resolve to pathogen taxids; prefer 2-word binomials
    names: list[str] = []
    for name, taxid in db["name_to_taxid"].items():
        if taxid in pathogen_taxids and len(name.split()) >= 2:
            names.append(name)
    names.sort(key=len, reverse=True)  # longest first → species beat genus
    return pathogen_taxids, names


def _find_meta_pathogen(text: str, pathogen_names: list[str]) -> str:
    """Return the first PHI-base pathogen name found in text (longest match first)."""
    if not text:
        return ""
    text_lower = text.lower()
    for name in pathogen_names:
        if name in text_lower:
            return name
    return ""


def _mismatch_flag(mode: str, sra: str, stat: str, meta: str) -> str:
    if not stat:
        return "stat_undetected"
    if mode == "mal":
        sra_genus  = _genus(sra)
        stat_genus = _genus(stat)
        meta_genus = _genus(meta)
        sra_stat_ok  = (sra_genus == stat_genus)
        sra_meta_ok  = (not meta) or (sra_genus == meta_genus)
        if sra_stat_ok and sra_meta_ok:
            return "all_match"
        if not sra_stat_ok and not sra_meta_ok:
            return "both_conflict"
        if not sra_stat_ok:
            return "stat_conflict"
        return "meta_conflict"
    else:  # hal — sra_organism is host, compare stat vs meta
        stat_genus = _genus(stat)
        meta_genus = _genus(meta)
        if not meta:
            return "no_meta"
        if stat_genus == meta_genus:
            return "all_match"
        return "stat_vs_meta_conflict"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(PHIBASE_DB) as f:
        db = json.load(f)
    _, pathogen_names = _build_meta_search_index(db)
    print(f"PHI-base search index: {len(pathogen_names):,} pathogen binomials")

    with open(RUNS_TSV, newline="") as f:
        all_runs = list(csv.DictReader(f, delimiter="\t"))
    runs = [r for r in all_runs if r.get("biosample_representative") == "True"]
    print(f"Loaded {len(all_runs):,} runs → {len(runs):,} biosample_representative")

    # Load BioProject metadata from biosample_kw (first row per BP)
    bp_meta: dict[str, dict] = {}
    if KW_TSV.exists():
        with open(KW_TSV, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                bp_meta.setdefault(row["BioProject"], row)

    # Pre-compute meta_pathogen per BioProject (expensive search done once per BP)
    bp_meta_pathogen: dict[str, str] = {}
    for bp, meta in bp_meta.items():
        text = (meta.get("title", "") + " " +
                meta.get("description", "") + " " +
                meta.get("abstract", ""))
        bp_meta_pathogen[bp] = _find_meta_pathogen(text, pathogen_names)

    output_rows = []
    flag_counts: Counter = Counter()

    mal_stat_conflict_ex: list[str] = []
    hal_stat_meta_ex: list[str] = []

    for r in runs:
        mode      = r["mode"]
        sra_org   = r["library_organism"]
        stat_top, stat_pct = _top_stat(r.get("stat_pathogens", ""))
        bp        = r["BioProject"]
        meta_path = bp_meta_pathogen.get(bp, "")
        bp_title  = bp_meta.get(bp, {}).get("title", "")

        flag = _mismatch_flag(mode, sra_org, stat_top, meta_path)
        flag_counts[f"{mode}:{flag}"] += 1

        if mode == "mal" and flag == "stat_conflict" and len(mal_stat_conflict_ex) < 5:
            mal_stat_conflict_ex.append(
                f"sra={sra_org!r} stat={stat_top!r} ({stat_pct}%) | {bp_title[:55]}")
        if mode == "hal" and flag == "stat_vs_meta_conflict" and len(hal_stat_meta_ex) < 5:
            hal_stat_meta_ex.append(
                f"stat={stat_top!r} meta={meta_path!r} | {bp_title[:55]}")

        output_rows.append({
            "mode":             mode,
            "Run":              r["Run"],
            "BioSample":        r["BioSample"],
            "BioProject":       bp,
            "biosample_rep":    r.get("biosample_representative", ""),
            "sra_organism":     sra_org,
            "stat_top_pathogen": stat_top,
            "stat_top_pct":     stat_pct,
            "meta_pathogen":    meta_path,
            "mismatch_flag":    flag,
            "bp_title":         bp_title,
        })

    # ── Summary ───────────────────────────────────────────────────────────────
    mal_total = sum(v for k, v in flag_counts.items() if k.startswith("mal:"))
    hal_total = sum(v for k, v in flag_counts.items() if k.startswith("hal:"))

    print(f"\n── MAL alignment ({mal_total:,} biosample_representative runs) ──")
    for flag in ("all_match", "stat_conflict", "meta_conflict",
                 "both_conflict", "stat_undetected"):
        n = flag_counts[f"mal:{flag}"]
        pct = n / max(mal_total, 1) * 100
        print(f"  {flag:<22}: {n:5d}  ({pct:.1f}%)")
    if mal_stat_conflict_ex:
        print("  stat_conflict examples:")
        for ex in mal_stat_conflict_ex:
            print(f"    {ex}")

    print(f"\n── HAL alignment ({hal_total:,} biosample_representative runs) ──")
    for flag in ("all_match", "stat_vs_meta_conflict", "no_meta", "stat_undetected"):
        n = flag_counts[f"hal:{flag}"]
        pct = n / max(hal_total, 1) * 100
        print(f"  {flag:<22}: {n:5d}  ({pct:.1f}%)")
    if hal_stat_meta_ex:
        print("  stat_vs_meta_conflict examples (STAT detects pathogen not in BP text):")
        for ex in hal_stat_meta_ex:
            print(f"    {ex}")

    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(output_rows)

    print(f"\nWritten: {OUT_TSV}  ({len(output_rows):,} rows)")


if __name__ == "__main__":
    main()
