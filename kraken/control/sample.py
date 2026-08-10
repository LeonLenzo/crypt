#!/usr/bin/env python3
"""
control/sample.py — build stratified validation set for Kraken2 pathogen-only DB.

Stratified sample of RNA-seq runs with STAT-confirmed ground truth, covering
fungi + oomycetes only (bacteria/viruses excluded from this DB).

Strata:
  A   True negative — non-plant host (mammal/bird/fish).  [STUB: provide --noplant-tsv]
  B   Plant host-only — gate-fail HAL runs, STAT 0% fungi + 0% oomycete.
  B2  Pathogen-only — MAL runs, host_pct=0, single STAT fungal/oomycete pathogen.
  C   Single pathogen, lab, low burden   (STAT euk 0.5–2%)
  D   Single pathogen, lab, medium burden(STAT euk 2–10%)
  E   Single pathogen, lab, high burden  (STAT euk >10%)
  F   Single pathogen, field             (STAT euk any %)
  G   Intentional co-infection           (llm_treatment=coinf_experiment, all)
  H   HC field co-infected, diff-genus   (co_infection_flag=T, same_genus_secondary=F, field)
  I   Same-genus pair                    (same_genus_secondary=T, all available)

Outputs:
  output/control/data/control_runs.tsv   — one row per run, with stratum + ground-truth cols
  output/control/data/run_ids.txt        — flat Run list for classify.py --run-list
  output/control/data/manifest.tsv       — stratum counts + target vs actual

Usage:
  python control/sample.py [--noplant-tsv PATH]
"""

import argparse
import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

SEED = 42
random.seed(SEED)

# ── sample targets ─────────────────────────────────────────────────────────────
N_TARGET = {
    "A":  200,   # non-plant true negative (external; requires --noplant-tsv)
    "B":  500,   # plant host-only (gate-fail HAL)
    "B2": 200,   # pathogen-only (MAL, host_pct=0.0)
    "C":  200,   # single pathogen, lab, low 0.5–2%
    "D":  200,   # single pathogen, lab, medium 2–10%
    "E":  200,   # single pathogen, lab, high >10%
    "F":  400,   # single pathogen, field (any %)
    "G":  None,  # all coinf_experiment
    "H":  400,   # HC field co-infected, diff-genus
    "I":  None,  # all same-genus pairs
}

MAX_PER_HOST = 50  # cap per named_host within each stratum

RNA_SOURCES = {
    "TRANSCRIPTOMIC",
    "TRANSCRIPTOMIC SINGLE CELL",
    "METATRANSCRIPTOMIC",
    "VIRAL RNA",
}

# NCBI taxids for kingdom-level detection in STAT tax_table
FUNGI_TAXID        = 4751
OOMYCOTA_TAXID     = 4762
VIRIDIPLANTAE_TAXID = 33090

# ── paths ──────────────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parents[2]
RUNS_TSV   = ROOT / "stat/output/data/runs.tsv"
KW_TSV     = ROOT / "metadata/output/data/biosample_kw.tsv"
LLM_TSV    = ROOT / "metadata/output/data/bioproject_llm.tsv"
STAT_CACHE = ROOT / "stat/output/data/stat_cache.jsonl"
HAL_RUNS   = ROOT / "stat/output/data/hal_runs.json"
MAL_RUNS   = ROOT / "stat/output/data/mal_runs.json"
OUT_DIR    = ROOT / "kraken/control/output/data"


# ── data loading ───────────────────────────────────────────────────────────────

def load_llm() -> dict[str, dict]:
    """Return {BioProject: {llm_treatment, llm_study_setting, ...}}."""
    out = {}
    with open(LLM_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["BioProject"]] = row
    return out


def load_named_hosts() -> dict[str, str]:
    """Return {BioSample: named_host} from biosample_kw.tsv."""
    out = {}
    with open(KW_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("named_host"):
                out[row["BioSample"]] = row["named_host"]
    return out


def load_enriched_runs(llm: dict, named_hosts: dict) -> list[dict]:
    """
    Load runs.tsv and enrich each row with:
      - named_host (from biosample_kw, falling back to host column)
      - euk_pct    (fungi_pct + oomycete_pct)
      - llm_treatment, llm_study_setting (from bioproject_llm)
    """
    rows = []
    with open(RUNS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            bp   = row.get("BioProject", "")
            bs   = row.get("BioSample", "")
            lrow = llm.get(bp, {})
            row["llm_treatment"]    = lrow.get("llm_treatment", "")
            row["llm_study_setting"] = lrow.get("llm_study_setting", "")
            row["named_host"]       = named_hosts.get(bs) or row.get("host", "")
            fungi   = float(row.get("fungi_pct")    or 0)
            oomyc   = float(row.get("oomycete_pct") or 0)
            row["euk_pct"] = round(fungi + oomyc, 4)
            rows.append(row)
    return rows


def _stat_pcts(entry: dict) -> tuple[float, float, float]:
    """
    From a parsed STAT cache entry dict, return (fungi_pct, oomycete_pct, host_pct).
    Uses NCBI kingdom/clade taxids in tax_table.
    """
    if not entry:
        return 0.0, 0.0, 0.0
    totals   = entry.get("tax_totals", {})
    analysed = totals.get("analysed", 0)
    if not analysed:
        return 0.0, 0.0, 0.0
    fungi_count = oomyc_count = host_count = 0
    for t in entry.get("tax_table", []):
        tid = t.get("tax_id")
        if tid == FUNGI_TAXID:
            fungi_count = t.get("total_count", 0)
        elif tid == OOMYCOTA_TAXID:
            oomyc_count = t.get("total_count", 0)
        elif tid == VIRIDIPLANTAE_TAXID:
            host_count  = t.get("total_count", 0)
    return (
        round(fungi_count / analysed * 100, 4),
        round(oomyc_count / analysed * 100, 4),
        round(host_count  / analysed * 100, 4),
    )


def load_gate_fail_plant_runs() -> list[dict]:
    """
    Stratum B source: HAL gate-fail runs where STAT fungi=0 AND oomycete=0.

    Strategy:
      1. Load hal_runs.json for RNA-source + ScientificName (host).
      2. Scan stat_cache.jsonl for those run IDs.
      3. Keep runs not in gate-pass set where both euk pcts are 0.
    """
    print("Loading HAL run metadata...", flush=True)
    with open(HAL_RUNS) as f:
        hal_meta = json.load(f)  # {Run: {RunInfo fields}}

    hal_rna = {
        run: meta
        for run, meta in hal_meta.items()
        if meta.get("LibrarySource", "").upper() in RNA_SOURCES
    }
    print(f"  HAL RNA runs: {len(hal_rna):,}", flush=True)

    # gate-pass run IDs from runs.tsv
    gate_pass = set()
    with open(RUNS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("mode") == "hal":
                gate_pass.add(row["Run"])
    print(f"  HAL gate-pass (exclude): {len(gate_pass):,}", flush=True)

    candidates = []
    hal_set = set(hal_rna)
    n_scanned = 0

    print("Scanning stat_cache for HAL gate-fail host-only runs...", flush=True)
    with open(STAT_CACHE) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            run, _, rest = line.partition("\t")
            if run not in hal_set or run in gate_pass:
                continue
            n_scanned += 1
            try:
                data = json.loads(rest)
                entry = data[0] if isinstance(data, list) else data
            except (json.JSONDecodeError, IndexError):
                continue
            fungi_pct, oomyc_pct, _ = _stat_pcts(entry)
            if fungi_pct == 0.0 and oomyc_pct == 0.0:
                meta = hal_rna[run]
                candidates.append({
                    "Run":        run,
                    "mode":       "hal",
                    "BioSample":  meta.get("BioSample", ""),
                    "BioProject": meta.get("BioProject", ""),
                    "named_host": meta.get("ScientificName", ""),
                    "host_pct":   "",
                    "euk_pct":    0.0,
                    "fungi_pct":  0.0,
                    "oomycete_pct": 0.0,
                    "co_infection_flag":    "False",
                    "same_genus_secondary": "False",
                    "stat_pathogens":       "",
                    "llm_treatment":        "",
                    "llm_study_setting":    "",
                    "stratum":    "B",
                    "expected":   "0% fungi/oomycete",
                })

    print(f"  Scanned {n_scanned:,} HAL gate-fail STAT entries; "
          f"{len(candidates):,} host-only candidates", flush=True)
    return candidates


def load_pathogen_only_runs() -> list[dict]:
    """
    Stratum B2 source: MAL gate-fail runs where STAT host_pct=0 and euk_pct > 5%.
    These are pure-culture pathogen libraries that fell below the MAL host gate.
    """
    print("Loading MAL run metadata...", flush=True)
    with open(MAL_RUNS) as f:
        mal_meta = json.load(f)

    mal_rna = {
        run: meta
        for run, meta in mal_meta.items()
        if meta.get("LibrarySource", "").upper() in RNA_SOURCES
    }
    print(f"  MAL RNA runs: {len(mal_rna):,}", flush=True)

    gate_pass = set()
    with open(RUNS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("mode") == "mal":
                gate_pass.add(row["Run"])
    print(f"  MAL gate-pass (exclude): {len(gate_pass):,}", flush=True)

    candidates = []
    mal_set = set(mal_rna)
    n_scanned = 0

    print("Scanning stat_cache for MAL gate-fail pathogen-only runs...", flush=True)
    with open(STAT_CACHE) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            run, _, rest = line.partition("\t")
            if run not in mal_set or run in gate_pass:
                continue
            n_scanned += 1
            try:
                data = json.loads(rest)
                entry = data[0] if isinstance(data, list) else data
            except (json.JSONDecodeError, IndexError):
                continue
            fungi_pct, oomyc_pct, host_pct = _stat_pcts(entry)
            euk_pct = round(fungi_pct + oomyc_pct, 4)
            if host_pct != 0.0 or euk_pct < 5.0:
                continue
            meta = mal_rna[run]
            candidates.append({
                "Run":        run,
                "mode":       "mal",
                "BioSample":  meta.get("BioSample", ""),
                "BioProject": meta.get("BioProject", ""),
                "named_host": meta.get("ScientificName", ""),
                "host_pct":   "0.0",
                "euk_pct":    euk_pct,
                "fungi_pct":  fungi_pct,
                "oomycete_pct": oomyc_pct,
                "co_infection_flag":    "single",
                "same_genus_secondary": "False",
                "stat_pathogens":       "",
                "llm_treatment":        "",
                "llm_study_setting":    "",
                "stratum":    "B2",
                "expected":   "single pathogen, 0% host reads",
            })

    print(f"  Scanned {n_scanned:,} MAL gate-fail STAT entries; "
          f"{len(candidates):,} pathogen-only candidates", flush=True)
    return candidates


# ── stratified sampler ─────────────────────────────────────────────────────────

def stratified_sample(rows: list[dict], n: int | None) -> list[dict]:
    """
    Sample up to n rows, capped at MAX_PER_HOST per named_host.
    If n is None, return all rows (still applying the per-host cap).
    """
    by_host: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_host[r.get("named_host") or "Unknown"].append(r)

    pool: list[dict] = []
    for host_rows in by_host.values():
        random.shuffle(host_rows)
        pool.extend(host_rows[:MAX_PER_HOST])

    random.shuffle(pool)
    return pool if n is None else pool[:n]


# ── stratum filters ────────────────────────────────────────────────────────────

def is_single_pathogen(row: dict) -> bool:
    return row.get("co_infection_flag") == "single"


def stratum_single_lab(rows: list[dict], lo: float, hi: float) -> list[dict]:
    """Single pathogen, lab setting, euk_pct in [lo, hi)."""
    out = []
    for r in rows:
        if not is_single_pathogen(r):
            continue
        if r.get("llm_study_setting") != "lab":
            continue
        if not (lo <= r["euk_pct"] < hi):
            continue
        out.append(r)
    return out


def stratum_single_field(rows: list[dict]) -> list[dict]:
    """Single pathogen, field setting, any euk_pct > 0."""
    out = []
    for r in rows:
        if not is_single_pathogen(r):
            continue
        if r.get("llm_study_setting") != "field":
            continue
        if r["euk_pct"] == 0.0:
            continue
        out.append(r)
    return out


def stratum_coinf_experiment(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r.get("llm_treatment") == "coinf_experiment"]


def stratum_hc_field_diffgenus(rows: list[dict]) -> list[dict]:
    return [
        r for r in rows
        if r.get("co_infection_flag") in ("multi_species", "multi_kingdom")
        and r.get("same_genus_secondary") == "False"
        and r.get("llm_study_setting") == "field"
    ]


def stratum_same_genus(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r.get("same_genus_secondary") == "True"]


# ── output columns ─────────────────────────────────────────────────────────────

OUT_COLS = [
    "stratum", "Run", "mode", "BioSample", "BioProject",
    "named_host", "host_pct", "euk_pct", "fungi_pct", "oomycete_pct",
    "stat_pathogens", "co_infection_flag", "same_genus_secondary",
    "llm_treatment", "llm_study_setting", "expected",
]


def tag(rows: list[dict], stratum: str, expected: str) -> list[dict]:
    for r in rows:
        r["stratum"]  = stratum
        r["expected"] = expected
    return rows


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--noplant-tsv", metavar="PATH",
                    help="TSV with Run column for non-plant true negatives (stratum A). "
                         "Must have Run and named_host columns.")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading LLM classifications...", flush=True)
    llm = load_llm()

    print("Loading named hosts...", flush=True)
    named_hosts = load_named_hosts()

    print("Loading enriched runs...", flush=True)
    all_runs = load_enriched_runs(llm, named_hosts)
    print(f"  Loaded {len(all_runs):,} runs", flush=True)

    # ── stratum A: non-plant true negative ────────────────────────────────────
    strat_a: list[dict] = []
    if args.noplant_tsv:
        with open(args.noplant_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                row.setdefault("euk_pct", 0.0)
                row.setdefault("fungi_pct", 0.0)
                row.setdefault("oomycete_pct", 0.0)
                row.setdefault("host_pct", "")
                row.setdefault("stat_pathogens", "")
                row.setdefault("co_infection_flag", "False")
                row.setdefault("same_genus_secondary", "False")
                row.setdefault("llm_treatment", "")
                row.setdefault("llm_study_setting", "")
                strat_a.append(row)
        strat_a = tag(stratified_sample(strat_a, N_TARGET["A"]), "A",
                      "0% fungi/oomycete (non-plant host)")
        print(f"  Stratum A: {len(strat_a)} runs", flush=True)
    else:
        print("  Stratum A: skipped (no --noplant-tsv provided)", flush=True)

    # ── stratum B: plant host-only (gate-fail HAL) ────────────────────────────
    print("Building stratum B (plant host-only)...", flush=True)
    b_candidates = load_gate_fail_plant_runs()
    strat_b = tag(stratified_sample(b_candidates, N_TARGET["B"]), "B",
                  "0% fungi/oomycete (plant host, no pathogen)")
    print(f"  Stratum B: {len(strat_b)} runs", flush=True)

    # ── stratum B2: pathogen-only (gate-fail MAL, host_pct=0) ────────────────
    print("Building stratum B2 (pathogen-only)...", flush=True)
    b2_cands = load_pathogen_only_runs()
    strat_b2 = tag(stratified_sample(b2_cands, N_TARGET["B2"]), "B2",
                   "single pathogen, 0% host reads")
    print(f"  Stratum B2: {len(strat_b2)} (from {len(b2_cands)} candidates)", flush=True)

    # ── strata C/D/E: single pathogen lab, abundance tiers ───────────────────
    strat_c = tag(stratified_sample(stratum_single_lab(all_runs, 0.5,  2.0), N_TARGET["C"]),
                  "C", "single pathogen, lab, low (0.5–2%)")
    strat_d = tag(stratified_sample(stratum_single_lab(all_runs, 2.0, 10.0), N_TARGET["D"]),
                  "D", "single pathogen, lab, medium (2–10%)")
    strat_e = tag(stratified_sample(stratum_single_lab(all_runs, 10.0, 1e9), N_TARGET["E"]),
                  "E", "single pathogen, lab, high (>10%)")
    print(f"  Stratum C: {len(strat_c)}, D: {len(strat_d)}, E: {len(strat_e)}", flush=True)

    # ── stratum F: single pathogen field ─────────────────────────────────────
    strat_f = tag(stratified_sample(stratum_single_field(all_runs), N_TARGET["F"]),
                  "F", "single pathogen, field")
    print(f"  Stratum F: {len(strat_f)}", flush=True)

    # ── stratum G: intentional co-infection (all) ─────────────────────────────
    g_cands = stratum_coinf_experiment(all_runs)
    strat_g = tag(stratified_sample(g_cands, N_TARGET["G"]), "G",
                  "intentional co-infection experiment")
    print(f"  Stratum G: {len(strat_g)} (from {len(g_cands)} candidates)", flush=True)

    # ── stratum H: HC field co-infected, diff-genus ───────────────────────────
    h_cands = stratum_hc_field_diffgenus(all_runs)
    strat_h = tag(stratified_sample(h_cands, N_TARGET["H"]), "H",
                  "HC field co-infected, diff-genus")
    print(f"  Stratum H: {len(strat_h)} (from {len(h_cands)} candidates)", flush=True)

    # ── stratum I: same-genus pair (all) ─────────────────────────────────────
    i_cands = stratum_same_genus(all_runs)
    strat_i = tag(stratified_sample(i_cands, N_TARGET["I"]), "I",
                  "same-genus secondary (LCA collapse test)")
    print(f"  Stratum I: {len(strat_i)} (from {len(i_cands)} candidates)", flush=True)

    # ── combine + deduplicate (Run appears in only one stratum) ───────────────
    strata_order = ["A", "B", "B2", "C", "D", "E", "F", "G", "H", "I"]
    all_strata: list[dict] = (
        strat_a + strat_b + strat_b2 +
        strat_c + strat_d + strat_e + strat_f +
        strat_g + strat_h + strat_i
    )

    seen_runs: set[str] = set()
    deduped: list[dict] = []
    for r in all_strata:
        if r["Run"] not in seen_runs:
            seen_runs.add(r["Run"])
            deduped.append(r)

    print(f"\nTotal unique runs: {len(deduped):,}", flush=True)

    # ── write control_runs.tsv ────────────────────────────────────────────────
    out_tsv = OUT_DIR / "control_runs.tsv"
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(deduped)
    print(f"Wrote {out_tsv}", flush=True)

    # ── write run_ids.txt ─────────────────────────────────────────────────────
    out_ids = OUT_DIR / "run_ids.txt"
    with open(out_ids, "w") as f:
        for r in deduped:
            f.write(r["Run"] + "\n")
    print(f"Wrote {out_ids}", flush=True)

    # ── write manifest.tsv ────────────────────────────────────────────────────
    counts: dict[str, int] = defaultdict(int)
    for r in deduped:
        counts[r["stratum"]] += 1

    out_manifest = OUT_DIR / "manifest.tsv"
    with open(out_manifest, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["stratum", "target", "actual", "description"])
        for s in strata_order:
            tgt = N_TARGET[s]
            w.writerow([s, tgt if tgt else "all", counts.get(s, 0),
                        next((r["expected"] for r in deduped if r["stratum"] == s), "")])
    print(f"Wrote {out_manifest}", flush=True)


if __name__ == "__main__":
    main()
