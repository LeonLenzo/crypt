#!/usr/bin/env python3
"""
02_filter.py — validate thresholds, apply retention gate, and detect cryptic
co-infections.  Produces a unified crypt.tsv combining MAL and HAL results.

Replaces three separate scripts (02_stat gate pass, 03_validate, 04_crypt):
  Phase 1  retention gate   — filter stat_cache.jsonl to confirmed runs
  Phase 2  validate         — empirical kingdom distribution tables (default on)
  Phase 3  crypt            — detect cryptic co-infecting pathogens

Unified output (output/02_filter/data/crypt.tsv) includes a `mode` column
(mal/hal) so MAL and HAL results can be analysed together or filtered
independently.  Run with --mode both (default) to process and merge both;
--mode mal or --mode hal for a single mode.

Cache compatibility
  Reads stat_cache.jsonl and {mode}_runs.json from output/01_fetch/data/ if
  present, falling back to output/02_stat/data/ and output/01_sra/data/ for
  backward compatibility with the legacy pipeline.

Usage:
  python 02_filter.py                     # both modes (default)
  python 02_filter.py --mode mal
  python 02_filter.py --mode hal
  python 02_filter.py --mode mal --skip-validate
"""

import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

from _util import _Tee, load_json, save_json

# ── Paths ──────────────────────────────────────────────────────────────────────

DB_PATH = Path("output/00_build/data/phibase_db.json")
OUT_DIR = Path("output/02_filter")


def _find_stat_cache() -> Path:
    for p in [Path("output/01_fetch/data/stat_cache.jsonl"),
              Path("output/02_stat/data/stat_cache.jsonl"),
              Path("output/02_stat/stat_cache.jsonl")]:
        if p.exists():
            return p
    raise FileNotFoundError(
        "stat_cache.jsonl not found in output/01_fetch/data/ or output/02_stat/\n"
        "Run: python 01_fetch.py --mode mal  (then hal)")


def _find_runs_path(mode: str) -> Path:
    for p in [Path(f"output/01_fetch/data/{mode}_runs.json"),
              Path(f"output/01_sra/data/{mode}_runs.json"),
              Path(f"output/01_sra/{mode}_runs.json")]:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"{mode}_runs.json not found in output/01_fetch/data/ or output/01_sra/\n"
        f"Run: python 01_fetch.py --mode {mode}")


# ── Shared constants ───────────────────────────────────────────────────────────

KINGDOMS = ["Fungi", "Viruses", "Bacteria", "Oomycota", "Nematoda"]

# Kingdom detection thresholds used by both validate (as reference) and crypt.
# Edit these if validate output suggests different values, then re-run with
# --skip-validate.
KINGDOM_THRESHOLDS: dict[str, float] = {
    "Fungi":    1.0,
    "Viruses":  10.0,
    "Bacteria": 5.0,
    "Oomycota": 0.5,
    "Nematoda": 1.0,
}

MAL_MIN_HOST_PCT     = 1.0   # Viridiplantae % gate for in-planta confirmation
HAL_MIN_PATHOGEN_PCT = 1.0   # min % to count a PHI-base pathogen as gate signal
ABS_MIN_PCT          = 0.1   # floor below which no organism is reported
LEAF_FRAC            = 0.75  # child must be >= this fraction of parent count
HOST_NODE            = "Viridiplantae"

STAT_KINGDOM_TO_DB_KEY: dict[str, str] = {
    "Fungi":    "fungal_to_seed",
    "Viruses":  "virus_to_seed",
    "Bacteria": "bacterial_to_seed",
    "Oomycota": "oomycete_to_seed",
    "Nematoda": "nematode_to_seed",
}
_ALL_KINGDOM_KEYS = tuple(STAT_KINGDOM_TO_DB_KEY.values())

MYCOVIRUS_KEYWORDS = {
    "mycovirus", "mitovirus", "hypovirus", "chrysovirus", "partitivirus",
    "totivirus", "endornavirus", "victorivirus", "botybirnavirus",
    "mitoviridae", "narnaviridae", "hypoviridae", "mycovirales",
}

# Validate breakpoints and absolute count floor candidates
BREAKPOINTS     = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
ABS_COUNT_FLOORS = [100, 250, 500, 1_000, 2_000, 5_000]


# ── STAT parsing utilities ────────────────────────────────────────────────────

def _analyzed(stat_data: list) -> int:
    t = stat_data[0].get("tax_totals", {})
    return t.get("analysed", 0) or t.get("analyzed", 0)


def _table(stat_data: list) -> list:
    return stat_data[0].get("tax_table", [])


def _node_count(table: list, node: str) -> int:
    for e in table:
        if e.get("org") == node:
            return e.get("total_count", 0)
    return 0


def _is_mycovirus(name: str) -> bool:
    n = name.lower()
    return any(kw in n for kw in MYCOVIRUS_KEYWORDS)


def _genus(name: str) -> str:
    parts = name.strip().split()
    return parts[0].lower() if parts else ""


def specific_hits(table: list, node: str, analyzed: int,
                  min_pct: float = ABS_MIN_PCT) -> list[tuple[str, float]]:
    """
    Return most-specific (leaf-level) organisms under a taxonomy node as
    (name, pct) tuples, sorted by pct descending.

    Leaf detection: a count is a leaf if no other count falls in the interval
    [LEAF_FRAC * count, count) — meaning no finer child captures most of the
    signal. Prefers 2-word binomials over clade names or strain strings.
    """
    top = _node_count(table, node)
    if not top:
        return []

    min_count = max(1, int(analyzed * min_pct / 100))
    under = [
        (e["org"], e["total_count"]) for e in table
        if 0 < e.get("total_count", 0) <= top
        and e.get("org") != node
        and e.get("total_count", 0) >= min_count
    ]
    if not under:
        return [(node, top / analyzed * 100)]

    count_vals  = sorted({c for _, c in under}, reverse=True)
    leaf_counts = {
        c for c in count_vals
        if not any(c2 < c and c2 >= LEAF_FRAC * c for c2 in count_vals)
    }

    best: dict[int, tuple[str, float]] = {}
    for name, count in under:
        if count not in leaf_counts:
            continue
        pct = count / analyzed * 100
        nw  = len(name.split())
        if count not in best:
            best[count] = (name, pct)
        elif nw == 2 and len(best[count][0].split()) != 2:
            best[count] = (name, pct)

    return sorted(best.values(), key=lambda x: -x[1])


def _parse_kingdom_pcts(stat_data: list) -> dict:
    """Return kingdom % dict + analyzed count from a STAT response."""
    if not stat_data or not isinstance(stat_data, list):
        return {}
    an  = _analyzed(stat_data)
    if not an:
        return {}
    tbl = _table(stat_data)
    pct = {e["org"]: e.get("total_count", 0) / an * 100
           for e in tbl if e.get("org")}
    return {
        "analyzed":     an,
        "host_pct":     pct.get(HOST_NODE,    0.0),
        "fungi_pct":    pct.get("Fungi",      0.0),
        "virus_pct":    pct.get("Viruses",    0.0),
        "bacteria_pct": pct.get("Bacteria",   0.0),
        "oomycete_pct": pct.get("Oomycota",   0.0),
        "nematode_pct": pct.get("Nematoda",   0.0),
        "_table":       tbl,
    }


# ── Phase 1: retention gate ───────────────────────────────────────────────────

def _passes_mal_gate(stat: dict) -> bool:
    return stat["host_pct"] >= MAL_MIN_HOST_PCT


def _passes_hal_gate(stat: dict, name_to_taxid: dict,
                     pathogen_taxids: set[int]) -> bool:
    an = stat["analyzed"]
    for entry in stat["_table"]:
        name = entry.get("org", "")
        pct  = entry.get("total_count", 0) / an * 100
        if pct < HAL_MIN_PATHOGEN_PCT:
            continue
        taxid = name_to_taxid.get(name.lower())
        if taxid and taxid in pathogen_taxids:
            return True
    return False


def apply_gate(jsonl_path: Path, runs: dict, mode: str,
               name_to_taxid: dict | None,
               pathogen_taxids: set | None) -> tuple[dict, int, int]:
    """
    Stream stat_cache.jsonl, apply the mode gate, return (confirmed, n_no_stat, n_fail).
    Only processes lines whose accession is in `runs`.
    """
    confirmed = {}
    n_no_stat = n_fail = 0

    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                acc, data_str = line.split("\t", 1)
                data = json.loads(data_str)
            except (ValueError, json.JSONDecodeError):
                continue

            run_row = runs.get(acc)
            if run_row is None:
                continue

            stat = _parse_kingdom_pcts(data)
            if not stat:
                n_no_stat += 1
                continue

            passes = (_passes_mal_gate(stat) if mode == "mal"
                      else _passes_hal_gate(stat, name_to_taxid, pathogen_taxids))
            if passes:
                row = dict(run_row)
                row.update({
                    "host_pct":     round(stat["host_pct"],     2),
                    "fungi_pct":    round(stat["fungi_pct"],    2),
                    "virus_pct":    round(stat["virus_pct"],    2),
                    "bacteria_pct": round(stat["bacteria_pct"], 2),
                    "oomycete_pct": round(stat["oomycete_pct"], 2),
                    "nematode_pct": round(stat["nematode_pct"], 2),
                    "analyzed":     stat["analyzed"],
                })
                confirmed[acc] = row
            else:
                n_fail += 1

    return confirmed, n_no_stat, n_fail


# ── Phase 2: validate distributions ──────────────────────────────────────────

def _percentile(sorted_vals: list[float], p: float) -> float:
    if not sorted_vals:
        return 0.0
    idx = (len(sorted_vals) - 1) * p / 100
    lo, hi = int(idx), min(int(idx) + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (idx - lo)


def _ascii_histogram(values: list[float], title: str, width: int = 50) -> str:
    pos = [v for v in values if v > 0]
    if not pos:
        return f"{title}: no values > 0\n"
    lo = math.floor(math.log2(min(pos))) if min(pos) >= 1 else -4
    hi = math.ceil(math.log2(max(pos))) + 1 if max(pos) >= 1 else 1
    edges  = [2 ** i for i in range(lo, hi + 1)]
    edges  = [e for e in edges if e <= max(pos) * 2]
    counts = [0] * (len(edges) - 1)
    for v in pos:
        for i in range(len(edges) - 1):
            if edges[i] <= v < edges[i + 1]:
                counts[i] += 1
                break
        else:
            counts[-1] += 1
    max_count = max(counts) if counts else 1
    lines = [f"\n{title} (n={len(pos):,} runs with signal, total={len(values):,})"]
    for i, c in enumerate(counts):
        label = f"{edges[i]:>7.2f}–{edges[i+1]:<7.2f}%"
        bar   = "█" * int(c / max_count * width)
        lines.append(f"  {label} | {bar} {c:,}")
    return "\n".join(lines) + "\n"


def _gap_threshold(sorted_vals: list[float]) -> float | None:
    pos = [v for v in sorted_vals if v > 0]
    if len(pos) < 10:
        return None
    log_vals = sorted(set(round(math.log10(v), 3) for v in pos))
    if len(log_vals) < 3:
        return None
    max_gap = 0.0
    gap_lo  = None
    for i in range(len(log_vals) - 1):
        gap = log_vals[i + 1] - log_vals[i]
        if gap > max_gap:
            max_gap = gap
            gap_lo  = 10 ** log_vals[i]
    return gap_lo


def _breakpoint_table(kingdom: str, kpcts: list[float],
                      analyzed_vals: list[int]) -> str:
    n_total    = len(kpcts)
    cur_thresh = KINGDOM_THRESHOLDS[kingdom]
    lines = [f"\n{kingdom} (n={n_total:,} confirmed runs, current threshold={cur_thresh}%)"]
    lines.append(
        f"  {'Thresh%':>8}  {'Runs≥thresh':>12}  {'% of total':>10}  "
        f"{'Median abs count':>17}  {'P10 abs count':>14}"
    )
    lines.append("  " + "-" * 68)
    for bp in BREAKPOINTS:
        passing   = [(p, a) for p, a in zip(kpcts, analyzed_vals) if p >= bp]
        n_passing = len(passing)
        if n_passing == 0:
            lines.append(
                f"  {bp:>7.2f}%  {0:>12,}  {0:>9.1f}%  {'—':>17}  {'—':>14}"
            )
            continue
        abs_counts = sorted(p * a / 100 for p, a in passing)
        med_abs    = _percentile(abs_counts, 50)
        p10_abs    = _percentile(abs_counts, 10)
        marker     = " ← current" if bp == cur_thresh else ""
        lines.append(
            f"  {bp:>7.2f}%  {n_passing:>12,}  {n_passing/n_total*100:>9.1f}%  "
            f"{med_abs:>17,.0f}  {p10_abs:>14,.0f}{marker}"
        )
    return "\n".join(lines)


def _abs_count_table(kingdom: str, species_rows: list[dict]) -> str:
    cur_thresh = KINGDOM_THRESHOLDS[kingdom]
    at_thresh  = [r for r in species_rows
                  if r["kingdom"] == kingdom and r["pct"] >= cur_thresh]
    if not at_thresh:
        return f"\n{kingdom}: no species detections at current threshold"
    abs_counts = sorted(r["abs_count"] for r in at_thresh)
    lines = [
        f"\n{kingdom} species-level detections at ≥{cur_thresh}% "
        f"(n={len(at_thresh):,}): abs count floor impact"
    ]
    lines.append(f"  {'Min abs count':>14}  {'Survives':>10}  {'% kept':>8}")
    lines.append("  " + "-" * 38)
    for floor in ABS_COUNT_FLOORS:
        n_survive = sum(1 for c in abs_counts if c >= floor)
        lines.append(
            f"  {floor:>14,}  {n_survive:>10,}  {n_survive/len(at_thresh)*100:>7.1f}%"
        )
    return "\n".join(lines)


def validate_phase(mode: str, confirmed: dict,
                   jsonl_path: Path) -> None:
    """
    Stream stat_cache for confirmed runs; print kingdom distribution tables
    and write {mode}_kingdom_dist.tsv + {mode}_species_dist.tsv.
    """
    needed = set(confirmed.keys())
    n_confirmed = len(needed)
    print(f"\n── Phase 2: validate ({mode.upper()}, {n_confirmed:,} confirmed runs) ──",
          flush=True)

    kingdom_rows: list[dict] = []
    species_rows: list[dict] = []
    n_found = n_nodata = 0

    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                acc, data_str = line.split("\t", 1)
            except ValueError:
                continue
            if acc not in needed:
                continue
            try:
                stat_data = json.loads(data_str)
            except json.JSONDecodeError:
                continue

            an = _analyzed(stat_data) if stat_data and isinstance(stat_data, list) else 0
            if not an:
                n_nodata += 1
                continue
            tbl = _table(stat_data)
            n_found += 1

            row = {"run": acc, "analyzed": an}
            for k in KINGDOMS:
                kc   = _node_count(tbl, k)
                kpct = kc / an * 100 if an else 0.0
                row[f"{k.lower()}_pct"] = round(kpct, 4)
                for name, pct in specific_hits(tbl, k, an):
                    if k == "Viruses" and _is_mycovirus(name):
                        continue
                    species_rows.append({
                        "run":       acc,
                        "kingdom":   k,
                        "name":      name,
                        "pct":       round(pct, 4),
                        "abs_count": round(pct * an / 100, 1),
                        "analyzed":  an,
                    })
            kingdom_rows.append(row)

    # Write TSVs
    kd_path = OUT_DIR / "data" / f"{mode}_kingdom_dist.tsv"
    sp_path = OUT_DIR / "data" / f"{mode}_species_dist.tsv"

    kd_fields = ["run", "analyzed"] + [f"{k.lower()}_pct" for k in KINGDOMS]
    with open(kd_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=kd_fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader(); w.writerows(kingdom_rows)

    sp_fields = ["run", "kingdom", "name", "pct", "abs_count", "analyzed"]
    with open(sp_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sp_fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader(); w.writerows(species_rows)

    print(f"  Written: {kd_path}  ({len(kingdom_rows):,} rows)", flush=True)
    print(f"  Written: {sp_path}  ({len(species_rows):,} rows)", flush=True)

    # Build and print summary
    analyzed_vals = [r["analyzed"] for r in kingdom_rows]
    lines = [
        f"\n── validate {mode.upper()} ─────────────────────────────",
        f"Confirmed runs processed: {n_found:,}  |  no STAT data: {n_nodata:,}",
        "", "=" * 72, "KINGDOM BREAKPOINT TABLES", "=" * 72,
    ]
    for k in KINGDOMS:
        kpcts = [r[f"{k.lower()}_pct"] for r in kingdom_rows]
        lines.append(_breakpoint_table(k, kpcts, analyzed_vals))

    lines += ["", "=" * 72, "KINGDOM DISTRIBUTIONS (log2 bins)", "=" * 72]
    for k in KINGDOMS:
        kpcts     = [r[f"{k.lower()}_pct"] for r in kingdom_rows]
        lines.append(_ascii_histogram(kpcts, k))
        gap_val   = _gap_threshold(sorted(v for v in kpcts if v > 0))
        if gap_val is not None:
            lines.append(
                f"  → Largest log-space gap below: {gap_val:.3f}%  "
                f"(suggested empirical threshold)"
            )
        else:
            lines.append("  → Distribution too sparse for gap detection")

    lines += ["", "=" * 72, "ABSOLUTE COUNT FLOOR (species-level)", "=" * 72]
    for k in KINGDOMS:
        lines.append(_abs_count_table(k, species_rows))

    if analyzed_vals:
        sa = sorted(analyzed_vals)
        lines += ["", "=" * 72, "READ DEPTH (all confirmed runs)", "=" * 72,
                  f"  p10={_percentile(sa,10):,.0f}  "
                  f"median={_percentile(sa,50):,.0f}  "
                  f"p90={_percentile(sa,90):,.0f}  "
                  f"max={max(sa):,.0f}"]

    summary = "\n".join(lines) + "\n"
    (OUT_DIR / "logs" / f"{mode}_validate_summary.txt").write_text(summary)
    print(summary, flush=True)
    print(f"  If thresholds need adjustment, edit KINGDOM_THRESHOLDS in "
          f"02_filter.py and re-run with --skip-validate.", flush=True)


# ── Phase 3: crypt — pathogen detection ───────────────────────────────────────

def _load_phibase() -> dict:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"PHI-base DB not found: {DB_PATH}\nRun:  python3 00_build.py")
    with open(DB_PATH) as f:
        raw = json.load(f)

    p_to_seed: dict[int, int] = {}
    kingdom_dicts: dict[str, dict[int, int]] = {}
    for key in _ALL_KINGDOM_KEYS:
        m = {int(k): v for k, v in raw[key].items()}
        kingdom_dicts[key] = m
        p_to_seed.update(m)

    db = {
        "name_to_taxid":  raw["name_to_taxid"],
        "taxid_to_name":  {int(k): v for k, v in raw["taxid_to_name"].items()},
        "pathogen_taxids": set(p_to_seed.keys()),
        "p_to_seed":      p_to_seed,
        "h_to_seed":      {int(k): v for k, v in raw["host_to_seed"].items()},
        "p_to_hosts":     {int(k): set(v) for k, v in raw["pathogen_to_hosts"].items()},
        **kingdom_dicts,
    }
    print(f"PHI-base: {len(db['pathogen_taxids']):,} pathogen taxids, "
          f"{len(db['name_to_taxid']):,} names", flush=True)
    return db


def detect_pathogens(stat_data: list, db: dict,
                     exclude_seed: int | None = None
                     ) -> list[tuple[int, str, float, str]]:
    """
    Return PHI-base pathogens detected in STAT above KINGDOM_THRESHOLDS.
    Each entry: (taxid, name, pct, kingdom).  Excludes same primary seed (MAL).
    """
    an = _analyzed(stat_data)
    if not an:
        return []
    tbl  = _table(stat_data)
    seen: dict[int, tuple[str, float, str]] = {}

    for kingdom, kthresh in KINGDOM_THRESHOLDS.items():
        kpct = _node_count(tbl, kingdom) / an * 100
        if kpct < kthresh:
            continue
        db_key      = STAT_KINGDOM_TO_DB_KEY[kingdom]
        kingdom_map = db[db_key]
        for name, pct in specific_hits(tbl, kingdom, an):
            if kingdom == "Viruses" and _is_mycovirus(name):
                continue
            taxid = db["name_to_taxid"].get(name.lower())
            if not taxid or taxid not in kingdom_map:
                continue
            seed = kingdom_map.get(taxid, taxid)
            if exclude_seed is not None and seed == exclude_seed:
                continue
            if taxid not in seen or pct > seen[taxid][1]:
                seen[taxid] = (name, pct, kingdom)

    return sorted(
        [(tid, nm, pct, kg) for tid, (nm, pct, kg) in seen.items()],
        key=lambda x: -x[2]
    )


def detect_host_species(stat_data: list) -> tuple[str, float]:
    an  = _analyzed(stat_data)
    if not an:
        return "", 0.0
    tbl  = _table(stat_data)
    hits = specific_hits(tbl, HOST_NODE, an)
    if not hits:
        return HOST_NODE, 0.0
    for name, pct in hits:
        if len(name.split()) == 2:
            return name, pct
    return hits[0]


def _interaction_status(p_taxid: int | None, h_taxid: int | None, db: dict) -> str:
    if p_taxid is None or h_taxid is None:
        return "unresolved"
    p_seed      = db["p_to_seed"].get(p_taxid, p_taxid)
    h_seed      = db["h_to_seed"].get(h_taxid, h_taxid)
    known_hosts = db["p_to_hosts"].get(p_seed, set())
    if h_seed in known_hosts:
        return "known"
    if known_hosts:
        return "novel_host_range"
    return "novel_combination"


def _classify_mal(run_row: dict, stat_data: list, db: dict) -> dict | None:
    an = _analyzed(stat_data)
    if not an:
        return None
    sra_name      = run_row.get("ScientificName", "") or run_row.get("Organism", "")
    primary_taxid = db["name_to_taxid"].get(sra_name.lower())
    primary_seed  = (db["p_to_seed"].get(primary_taxid, primary_taxid)
                     if primary_taxid else None)
    host_name, host_species_pct = detect_host_species(stat_data)
    host_taxid  = db["name_to_taxid"].get(host_name.lower())
    secondaries = detect_pathogens(stat_data, db, exclude_seed=primary_seed)
    kingdoms    = {kg for _, _, _, kg in secondaries}
    flag = ("multi_kingdom" if len(kingdoms) > 1
            else "multi_species" if secondaries else "single")
    pri_genus  = _genus(sra_name)
    same_genus = bool(pri_genus and any(_genus(nm) == pri_genus for _, nm, _, _ in secondaries))
    return {
        "host":                 host_name,
        "host_pct":             round(run_row.get("host_pct", 0), 2),
        "host_species_pct":     round(host_species_pct, 2),
        "primary_pathogen":     sra_name,
        "primary_taxid":        primary_taxid or "",
        "primary_pct":          "",
        "secondary_pathogens":  "; ".join(f"{nm} ({pct:.1f}%)" for _, nm, pct, _ in secondaries),
        "secondary_kingdoms":   "; ".join(sorted(kingdoms)),
        "n_secondary":          len(secondaries),
        "interaction_status":   _interaction_status(primary_taxid, host_taxid, db),
        "co_infection_flag":    flag,
        "same_genus_secondary": str(same_genus),
    }


def _classify_hal(run_row: dict, stat_data: list, db: dict) -> dict | None:
    an = _analyzed(stat_data)
    if not an:
        return None
    host_name      = run_row.get("ScientificName", "") or run_row.get("Organism", "")
    host_taxid_str = run_row.get("TaxID", "")
    try:
        host_taxid = int(host_taxid_str)
    except (ValueError, TypeError):
        host_taxid = db["name_to_taxid"].get(host_name.lower())
    pathogens = detect_pathogens(stat_data, db)
    if not pathogens:
        return None
    primary_taxid, primary_name, primary_pct, _ = pathogens[0]
    secondaries = pathogens[1:]
    kingdoms    = {pathogens[0][3]} | {kg for _, _, _, kg in secondaries}
    flag = ("multi_kingdom" if len(kingdoms) > 1
            else "multi_species" if secondaries else "single")
    pri_genus  = _genus(primary_name)
    same_genus = bool(pri_genus and any(_genus(nm) == pri_genus for _, nm, _, _ in secondaries))
    return {
        "host":                 host_name,
        "host_pct":             round(run_row.get("host_pct", 0), 2),
        "host_species_pct":     "",
        "primary_pathogen":     primary_name,
        "primary_taxid":        primary_taxid,
        "primary_pct":          round(primary_pct, 2),
        "secondary_pathogens":  "; ".join(f"{nm} ({pct:.1f}%)" for _, nm, pct, _ in secondaries),
        "secondary_kingdoms":   "; ".join(sorted({kg for _, _, _, kg in secondaries})),
        "n_secondary":          len(secondaries),
        "interaction_status":   _interaction_status(primary_taxid, host_taxid, db),
        "co_infection_flag":    flag,
        "same_genus_secondary": str(same_genus),
    }


def _load_stat_for_confirmed(jsonl_path: Path, needed: set[str]) -> dict:
    cache = {}
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                acc, data_str = line.split("\t", 1)
                if acc in needed:
                    cache[acc] = json.loads(data_str)
            except (ValueError, json.JSONDecodeError):
                continue
    return cache


# ── Output helpers ────────────────────────────────────────────────────────────

OUTPUT_FIELDS = [
    "Run", "mode", "BioSample", "BioProject", "SRAStudy", "Platform",
    "ScientificName",
    "host", "host_pct", "host_species_pct",
    "primary_pathogen", "primary_taxid", "primary_pct",
    "secondary_pathogens", "secondary_kingdoms", "n_secondary",
    "interaction_status", "co_infection_flag", "same_genus_secondary",
    "biosample_n_runs", "biosample_representative",
    "fungi_pct", "virus_pct", "bacteria_pct", "oomycete_pct", "nematode_pct",
    "analyzed",
]


def _annotate_biosamples(rows: list[dict]) -> None:
    """Add biosample_n_runs and biosample_representative in-place (across modes)."""
    groups: dict[str, list[int]] = {}
    for i, row in enumerate(rows):
        bs  = row.get("BioSample", "").strip()
        key = bs if bs else f"__noBS_{row['Run']}"
        groups.setdefault(key, []).append(i)
    for indices in groups.values():
        best = max(indices, key=lambda i: int(rows[i].get("analyzed") or 0))
        for i in indices:
            rows[i]["biosample_n_runs"]        = len(indices)
            rows[i]["biosample_representative"] = (i == best)


def write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _mode_summary(rows: list[dict], mode: str) -> str:
    flags    = Counter(r["co_infection_flag"] for r in rows)
    statuses = Counter(r["interaction_status"] for r in rows)
    repr_rows  = [r for r in rows if r.get("biosample_representative")]
    bs_flags   = Counter(r["co_infection_flag"] for r in repr_rows)
    bs_statuses = Counter(r["interaction_status"] for r in repr_rows)
    top_secondary: Counter = Counter()
    for r in repr_rows:
        for entry in r["secondary_pathogens"].split(";"):
            name = entry.strip().split("(")[0].strip()
            if name:
                top_secondary[name] += 1
    top_lines = "\n".join(
        f"  {name:<45} {n:>5,}" for name, n in top_secondary.most_common(15)
    )
    s = (
        f"{mode.upper()}:\n"
        f"  Runs classified:   {len(rows):>8,}\n"
        f"    single:          {flags['single']:>8,}\n"
        f"    multi_species:   {flags['multi_species']:>8,}\n"
        f"    multi_kingdom:   {flags['multi_kingdom']:>8,}\n"
        f"  Unique BioSamples: {len(repr_rows):>8,}\n"
        f"    single:          {bs_flags['single']:>8,}\n"
        f"    multi_species:   {bs_flags['multi_species']:>8,}\n"
        f"    multi_kingdom:   {bs_flags['multi_kingdom']:>8,}\n"
        f"  Interaction status (BioSample):\n"
        f"    known:           {bs_statuses['known']:>6,}\n"
        f"    novel_host_range:{bs_statuses['novel_host_range']:>6,}\n"
        f"    novel_combo:     {bs_statuses['novel_combination']:>6,}\n"
        f"    unresolved:      {bs_statuses['unresolved']:>6,}\n"
    )
    if top_lines:
        s += f"  Top secondary pathogens (BioSample):\n{top_lines}\n"
    return s


# ── Single-mode pipeline ──────────────────────────────────────────────────────

def run_mode(mode: str, db: dict, jsonl_path: Path,
             skip_validate: bool) -> tuple[list[dict], dict]:
    """
    Run gate + validate + crypt for one mode.
    Returns (classified_rows, gate_stats).
    """
    runs_path = _find_runs_path(mode)
    runs      = load_json(runs_path)
    print(f"\n── {mode.upper()}: {len(runs):,} runs from {runs_path} ──", flush=True)

    # Phase 1: gate
    print(f"\n── Phase 1: retention gate ({mode.upper()}) ──", flush=True)
    name_to_taxid = pathogen_taxids = None
    if mode == "hal":
        name_to_taxid   = db["name_to_taxid"]
        pathogen_taxids = db["pathogen_taxids"]

    confirmed, n_no_stat, n_fail = apply_gate(
        jsonl_path, runs, mode, name_to_taxid, pathogen_taxids)

    gate_desc = (f"≥{MAL_MIN_HOST_PCT}% Viridiplantae" if mode == "mal"
                 else f"PHI-base pathogen ≥{HAL_MIN_PATHOGEN_PCT}%")
    n_total     = len(runs)
    n_stat      = n_total - n_no_stat
    n_confirmed = len(confirmed)
    print(f"  Runs:        {n_total:,}", flush=True)
    print(f"  With STAT:   {n_stat:,}  ({n_stat/max(n_total,1)*100:.1f}%)", flush=True)
    print(f"  Failed gate: {n_fail:,}", flush=True)
    print(f"  Confirmed:   {n_confirmed:,}  "
          f"({n_confirmed/max(n_stat,1)*100:.1f}% of screened)  [{gate_desc}]",
          flush=True)

    confirmed_path = OUT_DIR / "data" / f"{mode}_confirmed.json"
    save_json(confirmed, confirmed_path)

    gate_stats = {
        "n_total": n_total, "n_stat": n_stat,
        "n_fail": n_fail, "n_confirmed": n_confirmed,
        "gate_desc": gate_desc,
    }

    # Phase 2: validate
    if not skip_validate:
        validate_phase(mode, confirmed, jsonl_path)

    # Phase 3: crypt
    print(f"\n── Phase 3: crypt ({mode.upper()}) ──", flush=True)
    stat_cache = _load_stat_for_confirmed(jsonl_path, set(confirmed.keys()))
    classify   = _classify_mal if mode == "mal" else _classify_hal
    rows: list[dict] = []
    skipped = 0

    for acc, run_row in confirmed.items():
        stat_data = stat_cache.get(acc)
        if not stat_data or not isinstance(stat_data, list):
            skipped += 1
            continue
        result = classify(run_row, stat_data, db)
        if result is None:
            skipped += 1
            continue
        row = {"Run": acc, "mode": mode}
        row.update(run_row)
        row.update(result)
        rows.append(row)

    order = {"multi_kingdom": 0, "multi_species": 1, "single": 2}
    rows.sort(key=lambda r: (order.get(r["co_infection_flag"], 3),
                             -int(r["n_secondary"] or 0)))

    print(f"  Classified: {len(rows):,}  |  skipped (no STAT): {skipped:,}",
          flush=True)
    return rows, gate_stats


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal", "both"], default="both",
                    help="mode(s) to process (default: both)")
    ap.add_argument("--skip-validate", action="store_true",
                    help="skip phase 2 distribution validation")
    args = ap.parse_args()

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "logs").mkdir(parents=True, exist_ok=True)

    log = _Tee(OUT_DIR / "logs" / "filter.log")
    sys.stdout = log

    try:
        target_modes = ["mal", "hal"] if args.mode == "both" else [args.mode]
        jsonl_path   = _find_stat_cache()
        print(f"STAT cache: {jsonl_path}", flush=True)

        db = _load_phibase()

        all_rows:   list[dict] = []
        gate_stats: dict[str, dict] = {}

        for mode in target_modes:
            rows, gs = run_mode(mode, db, jsonl_path,
                                skip_validate=args.skip_validate)
            all_rows.extend(rows)
            gate_stats[mode] = gs

        # Dedup: same Run accession in both modes (rare edge case)
        seen: set[str] = set()
        deduped: list[dict] = []
        for row in all_rows:
            acc = row["Run"]
            if acc not in seen:
                seen.add(acc)
                deduped.append(row)
        n_dupes = len(all_rows) - len(deduped)
        if n_dupes:
            print(f"\nDeduplicated {n_dupes} Run(s) present in multiple modes.",
                  flush=True)

        _annotate_biosamples(deduped)

        out_path = OUT_DIR / "data" / "crypt.tsv"
        write_tsv(deduped, out_path)

        # ── Summary ───────────────────────────────────────────────────────────
        mode_label = " + ".join(m.upper() for m in target_modes)
        summary_lines = [
            f"── 02_filter summary ────────────────────────────────",
            f"Modes: {mode_label}",
            "",
        ]

        for mode in target_modes:
            gs = gate_stats[mode]
            summary_lines.append(
                f"Gate ({mode.upper()}):  {gs['n_confirmed']:,} / {gs['n_total']:,} "
                f"confirmed  ({gs['n_confirmed']/max(gs['n_stat'],1)*100:.1f}%)  "
                f"[{gs['gate_desc']}]"
            )

        summary_lines += [""]

        for mode in target_modes:
            mode_rows = [r for r in deduped if r["mode"] == mode]
            summary_lines.append(_mode_summary(mode_rows, mode))

        if len(target_modes) > 1:
            summary_lines += [
                f"Combined:",
                f"  Total rows:          {len(deduped):,}",
                f"  Duplicates removed:  {n_dupes:,}",
            ]

        summary_lines += [f"\nOutput: {out_path}"]
        summary = "\n".join(summary_lines) + "\n"
        (OUT_DIR / "logs" / "filter_summary.txt").write_text(summary)
        print(f"\n{summary}")

    finally:
        log.close()


if __name__ == "__main__":
    main()
