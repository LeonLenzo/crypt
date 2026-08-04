#!/usr/bin/env python3
"""
02_filter_runs.py — gate, detect eukaryotic co-pathogens, annotate.

Scope: eukaryotic plant pathogens only (Fungi, Oomycota, Nematoda).
Bacteria and viruses are excluded from co-infection detection; their
kingdom-level percentages are retained in runs.tsv as informational
columns. Rationale: plant RNA-seq uses polyA+ selection, which depletes
bacterial mRNA (no polyA tails) and most plant virus families lacking 3'
polyA sequences (tospoviruses, tobamoviruses, cucumoviruses, luteoviruses),
making kingdom-level thresholds for those groups biologically
uninterpretable. Eukaryotic mRNA is faithfully captured by polyA
selection. See methods for full justification.

Three phases in a single STAT-cache pass per mode:
  Gate   — MAL: ≥1% Viridiplantae reads; HAL: any eukaryotic PHI-base
            pathogen ≥1%.
  Detect — leaf-level species via specific_hits(); PHI-base cross-ref;
            co-infection flag; interaction status; same-genus flag.
  Annotate — biosample dedup (highest-analysed run per BioSample).

Validate (skipped with --skip-validate):
  Kingdom distribution tables, breakpoint analysis, R-ready TSVs.
"""
import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

from _util import _Tee, link_latest, load_json, make_log_dir, save_json

# ── Detection scope ───────────────────────────────────────────────────────────

KINGDOM_THRESHOLDS: dict[str, float] = {
    "Fungi":    0.5,
    "Oomycota": 0.5,
    "Nematoda": 1.0,
}
# Maps STAT kingdom node name → phibase_db JSON key
EUK_KINGDOM_KEYS: dict[str, str] = {
    "Fungi":    "fungal_to_seed",
    "Oomycota": "oomycete_to_seed",
    "Nematoda": "nematode_to_seed",
}

MAL_MIN_HOST_PCT     = 1.0   # Viridiplantae % gate for MAL in-planta confirmation
HAL_MIN_PATHOGEN_PCT = 1.0   # eukaryotic pathogen % gate for HAL
ABS_MIN_PCT          = 0.1   # floor below which no organism is reported
LEAF_FRAC            = 0.75  # child must be ≥ this fraction of parent count
HOST_NODE            = "Viridiplantae"

_RNA_SOURCES = {
    "TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "METATRANSCRIPTOMIC", "VIRAL RNA"
}

# ── Paths ─────────────────────────────────────────────────────────────────────

OUT_DIR   = Path("output/02_filter_runs")
DB_PATH   = Path("output/00_build/data/phibase_db.json")
JSONL_PATH = Path("output/01_fetch_runs/data/stat_cache.jsonl")

BREAKPOINTS      = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
ABS_COUNT_FLOORS = [100, 250, 500, 1_000, 2_000, 5_000]


# ── DB loading ────────────────────────────────────────────────────────────────

def _load_db() -> dict:
    if not DB_PATH.exists():
        sys.exit(f"PHI-base DB not found: {DB_PATH}\nRun: python3 00_build.py")
    raw = json.loads(DB_PATH.read_text())

    euk_p_to_seed: dict[int, int] = {}
    kingdom_maps:  dict[str, dict[int, int]] = {}
    for kingdom, db_key in EUK_KINGDOM_KEYS.items():
        m = {int(k): v for k, v in raw[db_key].items()}
        kingdom_maps[db_key] = m
        euk_p_to_seed.update(m)

    # p_to_seed includes non-euk kingdoms so that bacterial/viral MAL library
    # organisms still resolve to their PHI-base seed for interaction_status.
    # Detection (detect_pathogens, HAL gate) uses euk_taxids only.
    all_p_to_seed = dict(euk_p_to_seed)
    for extra_key in ("bacterial_to_seed", "virus_to_seed"):
        if extra_key in raw:
            all_p_to_seed.update({int(k): v for k, v in raw[extra_key].items()})

    db = {
        "name_to_taxid":       raw["name_to_taxid"],
        "taxid_to_name":       {int(k): v for k, v in raw["taxid_to_name"].items()},
        "euk_taxids":          set(euk_p_to_seed),
        "p_to_seed":           all_p_to_seed,
        "h_to_seed":           {int(k): v for k, v in raw["host_to_seed"].items()},
        "p_to_hosts":          {int(k): set(v) for k, v in raw["pathogen_to_hosts"].items()},
        "viridiplantae_names": set(raw.get("viridiplantae_names", [])),
        **kingdom_maps,
    }
    print(
        f"PHI-base (euk scope): {len(euk_p_to_seed):,} pathogen taxids, "
        f"{len(raw['name_to_taxid']):,} names, "
        f"{len(db['viridiplantae_names']):,} Viridiplantae names",
        flush=True,
    )
    return db


# ── STAT parsing ──────────────────────────────────────────────────────────────

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


def _parse_kingdom_pcts(stat_data: list) -> dict | None:
    if not stat_data or not isinstance(stat_data, list):
        return None
    an = _analyzed(stat_data)
    if not an:
        return None
    tbl = _table(stat_data)
    pct = {e["org"]: e.get("total_count", 0) / an * 100 for e in tbl if e.get("org")}
    return {
        "analyzed":     an,
        "host_pct":     pct.get(HOST_NODE,  0.0),
        "fungi_pct":    pct.get("Fungi",    0.0),
        "virus_pct":    pct.get("Viruses",  0.0),  # informational only
        "bacteria_pct": pct.get("Bacteria", 0.0),  # informational only
        "oomycete_pct": pct.get("Oomycota", 0.0),
        "nematode_pct": pct.get("Nematoda", 0.0),
        "_table":       tbl,
    }


# ── specific_hits — core leaf-level detection algorithm ──────────────────────

def specific_hits(table: list, node: str, analyzed: int,
                  min_pct: float = ABS_MIN_PCT) -> list[tuple[str, float]]:
    """
    Return most-specific (leaf-level) organisms under a STAT taxonomy node as
    (name, pct) tuples sorted by pct descending.

    A count is a leaf if no other count falls in [LEAF_FRAC * count, count),
    meaning no finer child captures most of the signal. 2-word binomials are
    preferred over clade names or strain strings at the same count level.
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
        if count not in best or (len(name.split()) == 2
                                  and len(best[count][0].split()) != 2):
            best[count] = (name, pct)
    return sorted(best.values(), key=lambda x: -x[1])


# ── Gate functions ────────────────────────────────────────────────────────────

def _passes_mal_gate(kp: dict) -> bool:
    return kp["host_pct"] >= MAL_MIN_HOST_PCT


def _passes_hal_gate(kp: dict, db: dict) -> bool:
    """HAL gate: any eukaryotic PHI-base pathogen ≥ HAL_MIN_PATHOGEN_PCT."""
    an = kp["analyzed"]
    n2t = db["name_to_taxid"]
    euk = db["euk_taxids"]
    for entry in kp["_table"]:
        if entry.get("total_count", 0) / an * 100 < HAL_MIN_PATHOGEN_PCT:
            continue
        taxid = n2t.get(entry.get("org", "").lower())
        if taxid and taxid in euk:
            return True
    return False


# ── Pathogen and host detection ───────────────────────────────────────────────

def detect_pathogens(stat_data: list, db: dict,
                     exclude_seed: int | None = None
                     ) -> list[tuple[int, str, float, str]]:
    """
    Return eukaryotic PHI-base pathogens detected above KINGDOM_THRESHOLDS.
    Each entry: (taxid, name, pct, kingdom).
    exclude_seed: MAL library organism seed taxid — excluded from results.
    """
    an  = _analyzed(stat_data)
    if not an:
        return []
    tbl  = _table(stat_data)
    seen: dict[int, tuple[str, float, str]] = {}

    for kingdom, kthresh in KINGDOM_THRESHOLDS.items():
        if _node_count(tbl, kingdom) / an * 100 < kthresh:
            continue
        db_key = EUK_KINGDOM_KEYS[kingdom]
        kmap   = db[db_key]
        for name, pct in specific_hits(tbl, kingdom, an):
            taxid = db["name_to_taxid"].get(name.lower())
            if not taxid or taxid not in kmap:
                continue
            seed = kmap.get(taxid, taxid)
            if exclude_seed is not None and seed == exclude_seed:
                continue
            if taxid not in seen or pct > seen[taxid][1]:
                seen[taxid] = (name, pct, kingdom)

    return sorted(
        [(tid, nm, pct, kg) for tid, (nm, pct, kg) in seen.items()],
        key=lambda x: -x[2],
    )


def _plant_hits(stat_data: list, db: dict) -> list[tuple[str, float]]:
    """All Viridiplantae species detected, filtered by allowlist, pct desc."""
    an = _analyzed(stat_data)
    if not an:
        return []
    tbl   = _table(stat_data)
    hits  = specific_hits(tbl, HOST_NODE, an)
    if not hits:
        return []
    vnames  = db["viridiplantae_names"]
    n2t     = db["name_to_taxid"]
    p_seeds = db["euk_taxids"]
    h_seeds = set(db["h_to_seed"])
    result  = []
    for name, pct in hits:
        lower = name.lower()
        if vnames:
            if lower in vnames:
                result.append((name, pct))
        else:
            taxid = n2t.get(lower)
            if taxid is None or (taxid not in p_seeds and taxid in h_seeds):
                result.append((name, pct))
    return result


# ── Interaction status ────────────────────────────────────────────────────────

def _interaction_status(p_taxid: int | None, h_taxid: int | None, db: dict) -> str:
    if p_taxid is None or h_taxid is None:
        return "unresolved"
    p_seed      = db["p_to_seed"].get(p_taxid, p_taxid)
    h_seed      = db["h_to_seed"].get(h_taxid, h_taxid)
    known_hosts = db["p_to_hosts"].get(p_seed, set())
    if h_seed in known_hosts:
        return "known"
    return "novel_host_range" if known_hosts else "novel_combination"


# ── Per-run classification ────────────────────────────────────────────────────

def _genus(name: str) -> str:
    parts = name.strip().split()
    return parts[0].lower() if parts else ""


def _classify(run_row: dict, stat_data: list, kp: dict,
              db: dict, mode: str) -> dict | None:
    """Return co-infection columns for one gate-pass run, or None if unusable."""
    library_organism = run_row.get("ScientificName", "") or run_row.get("Organism", "")
    lib_lower = library_organism.lower()
    lib_taxid = db["name_to_taxid"].get(lib_lower)
    lib_seed  = db["p_to_seed"].get(lib_taxid) if lib_taxid else None

    all_pathogens = detect_pathogens(stat_data, db)
    all_hosts     = _plant_hits(stat_data, db)

    # HAL requires at least one eukaryotic pathogen detected above threshold
    if mode == "hal" and not all_pathogens:
        return None

    stat_pathogens_str = "; ".join(f"{nm}:{pct:.1f}%" for _, nm, pct, _ in all_pathogens)
    stat_hosts_str = (
        "; ".join(f"{nm}:{pct:.1f}%" for nm, pct in all_hosts)
        if all_hosts else
        (f"{HOST_NODE}:{kp['host_pct']:.1f}%" if kp["host_pct"] else "")
    )

    # host column: best species-level Viridiplantae hit (MAL) or library organism (HAL)
    if mode == "mal":
        best = next(((nm, pct) for nm, pct in all_hosts if len(nm.split()) == 2), None)
        host_name = best[0] if best else (all_hosts[0][0] if all_hosts else HOST_NODE)
    else:
        host_name = library_organism

    # library_detected
    if mode == "mal":
        lib_detected = lib_seed is not None and any(
            db["p_to_seed"].get(tid, tid) == lib_seed for tid, *_ in all_pathogens
        )
    else:
        lib_detected = bool(all_hosts and any(
            lib_lower in nm.lower() or nm.lower() == lib_lower for nm, _ in all_hosts
        ))

    # co-infection secondaries (MAL: exclude library organism seed; HAL: exclude top)
    if mode == "mal":
        co_pats   = [(tid, nm, pct, kg) for tid, nm, pct, kg in all_pathogens
                     if lib_seed is None or db["p_to_seed"].get(tid, tid) != lib_seed]
        pri_genus = _genus(library_organism)
    else:
        co_pats   = all_pathogens[1:]
        pri_genus = _genus(all_pathogens[0][1]) if all_pathogens else ""

    same_genus = bool(pri_genus and any(_genus(nm) == pri_genus for _, nm, *_ in co_pats))

    chk_kingdoms = {kg for _, _, _, kg in co_pats}
    if mode == "hal" and all_pathogens:
        chk_kingdoms.add(all_pathogens[0][3])
    flag = ("multi_kingdom" if len(chk_kingdoms) > 1
            else "multi_species" if co_pats else "single")

    # interaction_status: MAL → pathogen vs plant host; HAL → primary vs library host
    if mode == "mal":
        host_taxid = db["name_to_taxid"].get(host_name.lower())
        istatus    = _interaction_status(lib_taxid, host_taxid, db)
    else:
        try:
            host_taxid = int(run_row.get("TaxID", ""))
        except (ValueError, TypeError):
            host_taxid = db["name_to_taxid"].get(lib_lower)
        istatus = _interaction_status(
            all_pathogens[0][0] if all_pathogens else None, host_taxid, db
        )

    return {
        "host":                 host_name,
        "host_pct":             round(kp["host_pct"], 2),
        "library_organism":     library_organism,
        "library_detected":     str(lib_detected),
        "stat_pathogens":       stat_pathogens_str,
        "stat_hosts":           stat_hosts_str,
        "n_pathogens":          len(all_pathogens),
        "interaction_status":   istatus,
        "co_infection_flag":    flag,
        "same_genus_secondary": str(same_genus),
    }


# ── Biosample deduplication ───────────────────────────────────────────────────

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
            rows[i]["biosample_n_runs"]         = len(indices)
            rows[i]["biosample_representative"]  = (i == best)


# ── Validate output ───────────────────────────────────────────────────────────

def _percentile(sv: list[float], p: float) -> float:
    if not sv:
        return 0.0
    idx = (len(sv) - 1) * p / 100
    lo, hi = int(idx), min(int(idx) + 1, len(sv) - 1)
    return sv[lo] + (sv[hi] - sv[lo]) * (idx - lo)


def _ascii_histogram(values: list[float], title: str, width: int = 50) -> str:
    pos = [v for v in values if v > 0]
    if not pos:
        return f"{title}: no values > 0\n"
    lo    = math.floor(math.log2(min(pos))) if min(pos) >= 1 else -4
    hi    = math.ceil(math.log2(max(pos))) + 1 if max(pos) >= 1 else 1
    edges = [2 ** i for i in range(lo, hi + 1)]
    edges = [e for e in edges if e <= max(pos) * 2]
    counts    = [sum(1 for v in pos if edges[i] <= v < edges[i + 1])
                 for i in range(len(edges) - 1)]
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
    log_vals = sorted({round(math.log10(v), 3) for v in pos})
    if len(log_vals) < 3:
        return None
    max_gap, gap_lo = 0.0, None
    for i in range(len(log_vals) - 1):
        gap = log_vals[i + 1] - log_vals[i]
        if gap > max_gap:
            max_gap, gap_lo = gap, 10 ** log_vals[i]
    return gap_lo


def _breakpoint_table(kingdom: str, kpcts: list[float],
                      analyzed_vals: list[int]) -> str:
    n_total    = len(kpcts)
    cur_thresh = KINGDOM_THRESHOLDS[kingdom]
    lines = [f"\n{kingdom} (n={n_total:,} confirmed runs, current threshold={cur_thresh}%)"]
    lines.append(f"  {'Thresh%':>8}  {'Runs≥thresh':>12}  {'% of total':>10}  "
                 f"{'Median abs count':>17}  {'P10 abs count':>14}")
    lines.append("  " + "-" * 68)
    for bp in BREAKPOINTS:
        passing = [(p, a) for p, a in zip(kpcts, analyzed_vals) if p >= bp]
        if not passing:
            lines.append(f"  {bp:>7.2f}%  {'0':>12}  {'0.0%':>10}  {'—':>17}  {'—':>14}")
            continue
        abs_counts = sorted(p * a / 100 for p, a in passing)
        med_abs    = _percentile(abs_counts, 50)
        p10_abs    = _percentile(abs_counts, 10)
        marker     = " ← current" if bp == cur_thresh else ""
        lines.append(
            f"  {bp:>7.2f}%  {len(passing):>12,}  {len(passing)/n_total*100:>9.1f}%  "
            f"{med_abs:>17,.0f}  {p10_abs:>14,.0f}{marker}"
        )
    return "\n".join(lines)


def _write_validate(mode: str, kpct_rows: list[dict],
                    species_rows: list[dict], log_dir: Path) -> None:
    """Write kingdom_dist.tsv, species_dist.tsv, and validate_summary.txt."""
    print(f"\n── Phase 2: validate ({mode.upper()}, {len(kpct_rows):,} confirmed runs) ──",
          flush=True)
    data_dir = OUT_DIR / "data"

    kfields = ["Run", "fungi_pct", "oomycota_pct", "nematoda_pct", "analyzed"]
    kdist_path = data_dir / f"{mode}_kingdom_dist.tsv"
    with open(kdist_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=kfields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(kpct_rows)

    sdist_path = data_dir / f"{mode}_species_dist.tsv"
    with open(sdist_path, "w", newline="") as f:
        w = csv.DictWriter(f,
            fieldnames=["Run", "kingdom", "name", "pct", "abs_count"],
            delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(species_rows)

    lines = []
    for kingdom in KINGDOM_THRESHOLDS:
        col    = f"{kingdom.lower()}_pct"
        kpcts  = [r[col] for r in kpct_rows]
        a_vals = [r["analyzed"] for r in kpct_rows]
        lines.append(_ascii_histogram(kpcts, f"{kingdom} kingdom %"))
        lines.append(_breakpoint_table(kingdom, kpcts, a_vals))
        gap_val = _gap_threshold(sorted(v for v in kpcts if v > 0))
        lines.append(
            f"  → Largest log-space gap below: {gap_val:.3f}%"
            if gap_val is not None else
            "  → Distribution too sparse for gap detection"
        )
        lines.append("")
    summary = "\n".join(lines)
    summary_path = log_dir / f"{mode}_validate_summary.txt"
    summary_path.write_text(summary)
    print(summary, flush=True)
    print(f"  If thresholds need adjustment, edit KINGDOM_THRESHOLDS in "
          f"this script, then re-run with --skip-validate.", flush=True)


# ── Output ────────────────────────────────────────────────────────────────────

OUTPUT_FIELDS = [
    "Run", "mode", "BioSample", "BioProject", "SRAStudy", "Platform",
    "host", "host_pct",
    "library_organism", "library_detected",
    "stat_pathogens", "stat_hosts", "n_pathogens",
    "interaction_status", "co_infection_flag", "same_genus_secondary",
    "biosample_n_runs", "biosample_representative",
    "fungi_pct", "virus_pct", "bacteria_pct", "oomycete_pct", "nematode_pct",
    "analyzed",
]


def _write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _mode_summary(rows: list[dict], mode: str) -> str:
    flags     = Counter(r["co_infection_flag"] for r in rows)
    repr_rows = [r for r in rows if r.get("biosample_representative")]
    bs_flags  = Counter(r["co_infection_flag"] for r in repr_rows)
    bs_status = Counter(r["interaction_status"] for r in repr_rows)
    lib_det   = sum(1 for r in repr_rows if r.get("library_detected") == "True")
    top: Counter = Counter()
    for r in repr_rows:
        for ent in r.get("stat_pathogens", "").split(";"):
            nm = ent.strip().split(":")[0].strip()
            if nm:
                top[nm] += 1
    top_lines = "\n".join(
        f"  {nm:<45} {n:>5,}" for nm, n in top.most_common(15)
    )
    s = (
        f"{mode.upper()}:\n"
        f"  Runs classified:      {len(rows):>8,}\n"
        f"    single:             {flags['single']:>8,}\n"
        f"    multi_species:      {flags['multi_species']:>8,}\n"
        f"    multi_kingdom:      {flags['multi_kingdom']:>8,}\n"
        f"  Unique BioSamples:    {len(repr_rows):>8,}\n"
        f"    single:             {bs_flags['single']:>8,}\n"
        f"    multi_species:      {bs_flags['multi_species']:>8,}\n"
        f"    multi_kingdom:      {bs_flags['multi_kingdom']:>8,}\n"
        f"  Library organism detected by STAT: {lib_det:,} / {len(repr_rows):,}\n"
        f"  Interaction status (BioSample):\n"
        f"    known:            {bs_status['known']:>6,}\n"
        f"    novel_host_range: {bs_status['novel_host_range']:>6,}\n"
        f"    novel_combo:      {bs_status['novel_combination']:>6,}\n"
        f"    unresolved:       {bs_status['unresolved']:>6,}\n"
    )
    if top_lines:
        s += f"  Top detected pathogens (BioSample):\n{top_lines}\n"
    return s


# ── Single-pass mode processing ───────────────────────────────────────────────

def _find_runs_path(mode: str) -> Path:
    for base in (Path("output/01_fetch_runs/data"), Path("output/01_fetch_runs")):
        p = base / f"{mode}_runs.json"
        if p.exists():
            return p
    sys.exit(f"{mode}_runs.json not found under output/01_fetch_runs/")


def _process_mode(mode: str, db: dict, skip_validate: bool,
                  log_dir: Path) -> tuple[list[dict], dict]:
    """
    Single-pass through stat_cache.jsonl: gate + detect + classify.
    Accumulates validate data during the pass if not skip_validate.
    """
    runs_path = _find_runs_path(mode)
    runs      = load_json(runs_path)
    needed    = set(runs.keys())
    print(f"\n── {mode.upper()}: {len(needed):,} runs from {runs_path} ──", flush=True)
    gate_desc = (f"≥{MAL_MIN_HOST_PCT}% Viridiplantae" if mode == "mal"
                 else f"euk PHI-base pathogen ≥{HAL_MIN_PATHOGEN_PCT}%")
    print(f"\n── Phase 1+3: gate + detect ({mode.upper()}, [{gate_desc}]) ──", flush=True)

    rows:         list[dict] = []
    kpct_rows:    list[dict] = []   # validate: one row per confirmed run
    species_rows: list[dict] = []   # validate: one row per species detection
    n_wrong_source = n_no_stat = n_fail = n_skip = 0

    with open(JSONL_PATH) as f:
        for line in f:
            if "\t" not in line:
                continue
            acc, _, rest = line.partition("\t")
            if acc not in needed:
                continue

            run_row = runs[acc]
            if run_row.get("LibrarySource") not in _RNA_SOURCES:
                n_wrong_source += 1
                continue

            try:
                stat_data = json.loads(rest)
            except json.JSONDecodeError:
                n_no_stat += 1
                continue

            kp = _parse_kingdom_pcts(stat_data)
            if kp is None:
                n_no_stat += 1
                continue

            passes = (_passes_mal_gate(kp) if mode == "mal"
                      else _passes_hal_gate(kp, db))
            if not passes:
                n_fail += 1
                continue

            # Accumulate validate data (before classify, so even skip-classify runs count)
            if not skip_validate:
                kpct_rows.append({
                    "Run":          acc,
                    "fungi_pct":    kp["fungi_pct"],
                    "oomycota_pct": kp["oomycete_pct"],
                    "nematoda_pct": kp["nematode_pct"],
                    "analyzed":     kp["analyzed"],
                })
                for kingdom in KINGDOM_THRESHOLDS:
                    for name, pct in specific_hits(kp["_table"], kingdom, kp["analyzed"]):
                        species_rows.append({
                            "Run":       acc,
                            "kingdom":   kingdom,
                            "name":      name,
                            "pct":       round(pct, 4),
                            "abs_count": round(pct * kp["analyzed"] / 100),
                        })

            result = _classify(run_row, stat_data, kp, db, mode)
            if result is None:
                n_skip += 1
                continue

            row = {
                "Run":          acc,
                "mode":         mode,
                "BioSample":    run_row.get("BioSample",  ""),
                "BioProject":   run_row.get("BioProject", ""),
                "SRAStudy":     run_row.get("SRAStudy",   ""),
                "Platform":     run_row.get("Platform",   ""),
                "fungi_pct":    round(kp["fungi_pct"],    2),
                "virus_pct":    round(kp["virus_pct"],    2),
                "bacteria_pct": round(kp["bacteria_pct"], 2),
                "oomycete_pct": round(kp["oomycete_pct"], 2),
                "nematode_pct": round(kp["nematode_pct"], 2),
                "analyzed":     kp["analyzed"],
            }
            row.update(result)
            rows.append(row)

    n_total = len(needed)
    n_stat  = n_total - n_no_stat - n_wrong_source
    n_conf  = len(rows) + n_skip
    print(f"  Total:          {n_total:,}", flush=True)
    print(f"  Non-RNA source: {n_wrong_source:,} (excluded — GENOMIC/METAGENOMIC/OTHER)",
          flush=True)
    print(f"  With STAT:      {n_stat:,}  ({n_stat/max(n_total,1)*100:.1f}%)", flush=True)
    print(f"  Failed gate:    {n_fail:,}", flush=True)
    print(f"  Confirmed:      {n_conf:,}  ({n_conf/max(n_stat,1)*100:.1f}%)  [{gate_desc}]",
          flush=True)
    print(f"  Classified:     {len(rows):,}  |  skipped (no classify): {n_skip:,}",
          flush=True)

    rows.sort(key=lambda r: (
        {"multi_kingdom": 0, "multi_species": 1, "single": 2}.get(r["co_infection_flag"], 3),
        -int(r.get("n_pathogens") or 0),
    ))

    if not skip_validate:
        _write_validate(mode, kpct_rows, species_rows, log_dir)

    save_json({r["Run"]: r for r in rows},
              OUT_DIR / "data" / f"{mode}_confirmed.json")

    return rows, {
        "n_total": n_total, "n_wrong_source": n_wrong_source, "n_stat": n_stat,
        "n_fail": n_fail, "n_confirmed": n_conf, "gate_desc": gate_desc,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal", "both"], default="both",
                    help="mode(s) to process (default: both)")
    ap.add_argument("--skip-validate", action="store_true",
                    help="skip phase 2 distribution validation (faster re-runs)")
    args = ap.parse_args()

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "filter.log")
    link_latest(logs_base, log_dir / "filter.log")
    sys.stdout = log

    try:
        if not JSONL_PATH.exists():
            sys.exit(f"stat_cache.jsonl not found: {JSONL_PATH}")
        print(f"STAT cache: {JSONL_PATH}", flush=True)

        db           = _load_db()
        target_modes = ["mal", "hal"] if args.mode == "both" else [args.mode]
        all_rows:    list[dict] = []
        gate_stats:  dict[str, dict] = {}

        for mode in target_modes:
            rows, gs = _process_mode(mode, db, args.skip_validate, log_dir)
            all_rows.extend(rows)
            gate_stats[mode] = gs

        # Dedup cross-mode Run collisions (rare edge case)
        seen: set[str] = set()
        deduped: list[dict] = []
        for row in all_rows:
            if row["Run"] not in seen:
                seen.add(row["Run"])
                deduped.append(row)
        n_dupes = len(all_rows) - len(deduped)
        if n_dupes:
            print(f"\nDeduplicated {n_dupes} Run(s) present in multiple modes.", flush=True)

        _annotate_biosamples(deduped)
        out_path = OUT_DIR / "data" / "runs.tsv"
        _write_tsv(deduped, out_path)

        # Summary
        mode_label = " + ".join(m.upper() for m in target_modes)
        lines = [f"── 02_filter_runs summary ────────────────────────────────",
                 f"Modes: {mode_label}", ""]
        for mode in target_modes:
            gs = gate_stats[mode]
            lines.append(
                f"Gate ({mode.upper()}):  {gs['n_confirmed']:,} / {gs['n_total']:,} "
                f"confirmed  ({gs['n_confirmed']/max(gs['n_stat'],1)*100:.1f}%)  "
                f"[{gs['gate_desc']}]"
            )
        lines.append("")
        for mode in target_modes:
            lines.append(_mode_summary([r for r in deduped if r["mode"] == mode], mode))
        if len(target_modes) > 1:
            lines += [
                "Combined:",
                f"  Total rows:          {len(deduped):,}",
                f"  Duplicates removed:  {n_dupes:,}",
            ]
        lines.append(f"\nOutput: {out_path}  ({len(deduped):,} rows)")
        summary = "\n".join(lines) + "\n"
        summary_path = log_dir / "filter_summary.txt"
        summary_path.write_text(summary)
        link_latest(logs_base, summary_path)
        print(f"\n{summary}")

    finally:
        log.close()


if __name__ == "__main__":
    main()
