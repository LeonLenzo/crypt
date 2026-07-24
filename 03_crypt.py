#!/usr/bin/env python3
"""
03_crypt.py — identify cryptic co-infections in confirmed mixed libraries.

Reads confirmed runs from cache/{mode}_confirmed.json (written by 02_stat.py)
and the shared STAT cache, then scans each run's full STAT taxonomy for
PHI-base plant pathogens to classify co-infection status.

MAL (microbe-as-library):
  host              = most abundant Viridiplantae leaf species in STAT
  primary_pathogen  = library organism (from SRA ScientificName)
  secondary         = other PHI-base pathogens detected in STAT

HAL (host-as-library):
  host              = library organism (from SRA ScientificName / TaxID)
  primary_pathogen  = most abundant PHI-base pathogen detected in STAT
  secondary         = remaining PHI-base pathogens, ordered by abundance

co_infection_flag:
  single        — one PHI-base pathogen detected (primary only)
  multi_species — 2+ pathogens, same kingdom
  multi_kingdom — 2+ pathogens from different kingdoms

Output: candidates/{mode}_crypt.tsv

Usage:
  python 03_crypt.py --mode mal
  python 03_crypt.py --mode hal
"""

import argparse
import csv
import json
import sys
from pathlib import Path

# ── Settings ──────────────────────────────────────────────────────────────────

DB_PATH         = Path("output/00_build/phibase_db.json")
STAT_CACHE_PATH = Path("output/02_stat/stat_cache.json")
OUT_DIR         = Path("output/03_crypt")


# ── Logging ───────────────────────────────────────────────────────────────────

class _Tee:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._log    = open(path, "w", buffering=1)
        self._stdout = sys.__stdout__

    def write(self, text: str) -> None:
        self._stdout.write(text)
        self._log.write(text)

    def flush(self) -> None:
        self._stdout.flush()
        self._log.flush()

    def close(self) -> None:
        sys.stdout = self._stdout
        self._log.close()

# Kingdom detection thresholds (% of analysed reads)
KINGDOM_THRESHOLDS = {
    "Fungi":    1.0,
    "Viruses":  5.0,
    "Bacteria": 5.0,
    "Oomycota": 0.5,
}

ABS_MIN_PCT = 0.1   # floor below which no organism is reported (% of analysed)
LEAF_FRAC   = 0.75  # child must be >= this fraction of parent to suppress parent

MYCOVIRUS_KEYWORDS = {
    "mycovirus", "mitovirus", "hypovirus", "chrysovirus", "partitivirus",
    "totivirus", "endornavirus", "victorivirus", "botybirnavirus",
    "mitoviridae", "narnaviridae", "hypoviridae", "mycovirales",
}

HOST_NODE = "Viridiplantae"


# ── PHI-base loader ───────────────────────────────────────────────────────────

def _load_phibase(db_path: Path) -> dict:
    if not db_path.exists():
        raise FileNotFoundError(
            f"PHI-base DB not found: {db_path}\nRun:  python3 00_build.py")
    with open(db_path) as f:
        raw = json.load(f)
    db = {
        "name_to_taxid":    raw["name_to_taxid"],
        "taxid_to_name":    {int(k): v for k, v in raw["taxid_to_name"].items()},
        "pathogen_taxids":  set(raw["pathogen_taxids"]),
        "virus_taxids":     set(raw.get("virus_taxids", [])),
        "p_to_seed":        {int(k): v for k, v in
                             raw["pathogen_taxid_to_seed"].items()},
        "h_to_seed":        {int(k): v for k, v in
                             raw["host_taxid_to_seed"].items()},
        "p_to_hosts":       {int(k): set(v) for k, v in
                             raw["pathogen_to_hosts"].items()},
    }
    n_virus = len(db["virus_taxids"])
    print(f"PHI-base: {len(db['pathogen_taxids']):,} pathogen taxids "
          f"(incl. {n_virus:,} virus), {len(db['name_to_taxid']):,} names",
          flush=True)
    return db


# ── STAT helpers ──────────────────────────────────────────────────────────────

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


def specific_hits(table: list, node: str, analyzed: int,
                  min_pct: float = ABS_MIN_PCT) -> list[tuple[str, float]]:
    """
    Return the most specific (leaf-level) organisms under a taxonomy node
    as (name, pct) tuples, sorted by pct descending.

    Leaf detection: a count is a leaf if no other count falls in the interval
    [LEAF_FRAC * count, count) — meaning there is no finer-grained child
    that captures most of its signal. Prefers 2-word binomials over
    clade names or strain-level strings.
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


# ── Pathogen detection ────────────────────────────────────────────────────────

def detect_pathogens(stat_data: list, db: dict,
                     exclude_seed: int | None = None
                     ) -> list[tuple[int, str, float, str]]:
    """
    Return PHI-base pathogens detected in STAT above kingdom thresholds.
    Each entry: (taxid, name, pct, kingdom).
    Excludes any pathogen whose seed taxid == exclude_seed (for MAL,
    to avoid re-listing the primary/library organism as a secondary).
    Sorted by pct descending.
    """
    analyzed = _analyzed(stat_data)
    if not analyzed:
        return []
    table = _table(stat_data)

    seen: dict[int, tuple[str, float, str]] = {}

    for kingdom, kthresh in KINGDOM_THRESHOLDS.items():
        kpct = _node_count(table, kingdom) / analyzed * 100
        if kpct < kthresh:
            continue

        for name, pct in specific_hits(table, kingdom, analyzed):
            if kingdom == "Viruses" and _is_mycovirus(name):
                continue
            taxid = db["name_to_taxid"].get(name.lower())
            if not taxid or taxid not in db["pathogen_taxids"]:
                continue
            seed = db["p_to_seed"].get(taxid, taxid)
            if exclude_seed is not None and seed == exclude_seed:
                continue
            if taxid not in seen or pct > seen[taxid][1]:
                seen[taxid] = (name, pct, kingdom)

    return sorted(
        [(tid, nm, pct, kg) for tid, (nm, pct, kg) in seen.items()],
        key=lambda x: -x[2]
    )


def detect_host_species(stat_data: list) -> tuple[str, float]:
    """
    Return the most abundant Viridiplantae leaf species in STAT as (name, pct).
    Falls back to ("Viridiplantae", host_pct) if no species-level hit found.
    """
    analyzed = _analyzed(stat_data)
    if not analyzed:
        return "", 0.0
    table = _table(stat_data)
    hits  = specific_hits(table, HOST_NODE, analyzed)
    if not hits:
        return HOST_NODE, 0.0
    # Prefer a 2-word binomial
    for name, pct in hits:
        if len(name.split()) == 2:
            return name, pct
    return hits[0]


# ── PHI-base interaction check ────────────────────────────────────────────────

def known_interaction(p_taxid: int, h_taxid: int, db: dict) -> bool:
    p_seed = db["p_to_seed"].get(p_taxid, p_taxid)
    h_seed = db["h_to_seed"].get(h_taxid, h_taxid)
    return h_seed in db["p_to_hosts"].get(p_seed, set())


# ── Per-run classification ────────────────────────────────────────────────────

def classify_mal(run_row: dict, stat_data: list, db: dict) -> dict | None:
    """
    MAL: library organism is the primary pathogen.
    Look for additional PHI-base pathogens as secondary/cryptic co-infectors.
    """
    analyzed = _analyzed(stat_data)
    if not analyzed:
        return None

    sra_name   = run_row.get("ScientificName", "") or run_row.get("Organism", "")
    primary_taxid = db["name_to_taxid"].get(sra_name.lower())
    primary_seed  = (db["p_to_seed"].get(primary_taxid, primary_taxid)
                     if primary_taxid else None)

    # Host: most abundant Viridiplantae leaf in STAT
    host_name, host_species_pct = detect_host_species(stat_data)
    host_taxid = db["name_to_taxid"].get(host_name.lower())

    # Secondary pathogens: PHI-base pathogens detected in STAT excluding the primary
    secondaries = detect_pathogens(stat_data, db, exclude_seed=primary_seed)

    ki = False
    if primary_taxid and host_taxid:
        ki = known_interaction(primary_taxid, host_taxid, db)
    elif primary_taxid is None:
        ki = None   # primary not in PHI-base name map

    kingdoms_detected = {kg for _, _, _, kg in secondaries}
    if len(secondaries) == 0:
        flag = "single"
    elif len(kingdoms_detected) > 1:
        flag = "multi_kingdom"
    else:
        flag = "multi_species"

    return {
        "host":                host_name,
        "host_pct":            round(run_row.get("host_pct", 0), 2),
        "host_species_pct":    round(host_species_pct, 2),
        "primary_pathogen":    sra_name,
        "primary_taxid":       primary_taxid or "",
        "primary_pct":         "",   # library organism dominates; STAT pct uninformative
        "secondary_pathogens": "; ".join(
            f"{nm} ({pct:.1f}%)" for _, nm, pct, _ in secondaries),
        "secondary_kingdoms":  "; ".join(sorted(kingdoms_detected)),
        "n_secondary":         len(secondaries),
        "known_interaction":   ("" if ki is None else str(ki)),
        "co_infection_flag":   flag,
    }


def classify_hal(run_row: dict, stat_data: list, db: dict) -> dict | None:
    """
    HAL: library organism is the host. Primary pathogen is the most abundant
    PHI-base pathogen detected in STAT.
    """
    analyzed = _analyzed(stat_data)
    if not analyzed:
        return None

    # Host: library organism from SRA metadata
    host_name  = run_row.get("ScientificName", "") or run_row.get("Organism", "")
    host_taxid_str = run_row.get("TaxID", "")
    try:
        host_taxid = int(host_taxid_str)
    except (ValueError, TypeError):
        host_taxid = db["name_to_taxid"].get(host_name.lower())

    pathogens = detect_pathogens(stat_data, db)
    if not pathogens:
        return None   # no PHI-base pathogens detected — skip run

    primary_taxid, primary_name, primary_pct, primary_kg = pathogens[0]
    secondaries = pathogens[1:]

    ki = False
    if host_taxid:
        ki = known_interaction(primary_taxid, host_taxid, db)

    kingdoms_detected = {pathogens[0][3]} | {kg for _, _, _, kg in secondaries}
    if len(pathogens) == 1:
        flag = "single"
    elif len(kingdoms_detected) > 1:
        flag = "multi_kingdom"
    else:
        flag = "multi_species"

    return {
        "host":                host_name,
        "host_pct":            round(run_row.get("host_pct", 0), 2),
        "host_species_pct":    "",
        "primary_pathogen":    primary_name,
        "primary_taxid":       primary_taxid,
        "primary_pct":         round(primary_pct, 2),
        "secondary_pathogens": "; ".join(
            f"{nm} ({pct:.1f}%)" for _, nm, pct, _ in secondaries),
        "secondary_kingdoms":  "; ".join(
            sorted({kg for _, _, _, kg in secondaries})),
        "n_secondary":         len(secondaries),
        "known_interaction":   str(ki),
        "co_infection_flag":   flag,
    }


# ── TSV output ────────────────────────────────────────────────────────────────

OUTPUT_FIELDS = [
    "Run", "BioProject", "SRAStudy", "Platform", "ScientificName",
    "host", "host_pct", "host_species_pct",
    "primary_pathogen", "primary_taxid", "primary_pct",
    "secondary_pathogens", "secondary_kingdoms", "n_secondary",
    "known_interaction", "co_infection_flag",
    "fungi_pct", "virus_pct", "bacteria_pct", "oomycete_pct", "analyzed",
]


def write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def print_summary(rows: list[dict], mode: str) -> None:
    from collections import Counter
    flags  = Counter(r["co_infection_flag"] for r in rows)
    ki_yes = sum(1 for r in rows if r["known_interaction"] == "True")
    ki_no  = sum(1 for r in rows if r["known_interaction"] == "False")

    print(f"\n── {mode.upper()} crypt summary ──")
    print(f"  Runs classified:     {len(rows):>8,}")
    print(f"  single:              {flags['single']:>8,}")
    print(f"  multi_species:       {flags['multi_species']:>8,}")
    print(f"  multi_kingdom:       {flags['multi_kingdom']:>8,}")
    print(f"  known interaction:   {ki_yes:>8,}  (PHI-base confirmed host×pathogen)")
    print(f"  novel interaction:   {ki_no:>8,}  (not in PHI-base)")

    if any(r["n_secondary"] for r in rows):
        top_secondary: Counter = Counter()
        for r in rows:
            for entry in r["secondary_pathogens"].split(";"):
                name = entry.strip().split("(")[0].strip()
                if name:
                    top_secondary[name] += 1
        print(f"\n  Top secondary pathogens:")
        for name, n in top_secondary.most_common(15):
            print(f"    {name:<45} {n:>5,}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal"], required=True,
                    help="mal = microbe-as-library  |  hal = host-as-library")
    args = ap.parse_args()

    log = _Tee(OUT_DIR / f"{args.mode}.log")
    sys.stdout = log

    try:
        confirmed_path = Path(f"output/02_stat/{args.mode}_confirmed.json")
        out_path       = OUT_DIR / f"{args.mode}_crypt.tsv"

        print(f"Mode: {args.mode.upper()}", flush=True)

        if not confirmed_path.exists():
            raise FileNotFoundError(
                f"{confirmed_path} not found — run 02_stat.py --mode {args.mode} first")

        confirmed  = json.loads(confirmed_path.read_text())
        stat_cache = json.loads(STAT_CACHE_PATH.read_text()) if STAT_CACHE_PATH.exists() else {}
        db         = _load_phibase(DB_PATH)

        print(f"Confirmed runs:  {len(confirmed):,}", flush=True)

        classify = classify_mal if args.mode == "mal" else classify_hal
        rows     = []
        skipped  = 0

        for acc, run_row in confirmed.items():
            stat_data = stat_cache.get(acc)
            if not stat_data or not isinstance(stat_data, list):
                skipped += 1
                continue
            result = classify(run_row, stat_data, db)
            if result is None:
                skipped += 1
                continue
            row = {"Run": acc}
            row.update(run_row)
            row.update(result)
            rows.append(row)

        order = {"multi_kingdom": 0, "multi_species": 1, "single": 2}
        rows.sort(key=lambda r: (order.get(r["co_infection_flag"], 3),
                                 -int(r["n_secondary"] or 0)))

        write_tsv(rows, out_path)
        print_summary(rows, args.mode)
        if skipped:
            print(f"\n  Skipped (no STAT / no pathogens detected): {skipped:,}")

        from collections import Counter
        flags  = Counter(r["co_infection_flag"] for r in rows)
        ki_yes = sum(1 for r in rows if r["known_interaction"] == "True")
        ki_no  = sum(1 for r in rows if r["known_interaction"] == "False")

        top_secondary: Counter = Counter()
        for r in rows:
            for entry in r["secondary_pathogens"].split(";"):
                name = entry.strip().split("(")[0].strip()
                if name:
                    top_secondary[name] += 1

        top_lines = "\n".join(
            f"  {name:<45} {n:>5,}"
            for name, n in top_secondary.most_common(15)
        )

        summary = (
            f"── 03_crypt {args.mode.upper()} summary ────────────────────────\n"
            f"Mode:              {args.mode.upper()}\n"
            f"\n"
            f"Runs classified:   {len(rows):>8,}\n"
            f"  single:          {flags['single']:>8,}\n"
            f"  multi_species:   {flags['multi_species']:>8,}\n"
            f"  multi_kingdom:   {flags['multi_kingdom']:>8,}\n"
            f"Skipped:           {skipped:>8,}\n"
            f"\n"
            f"Known interaction (PHI-base): {ki_yes:>6,}\n"
            f"Novel interaction:            {ki_no:>6,}\n"
        )
        if top_lines:
            summary += f"\nTop secondary pathogens:\n{top_lines}\n"
        summary += (
            f"\nOutput: {out_path}\n"
            f"Log:    {OUT_DIR / f'{args.mode}.log'}\n"
        )
        (OUT_DIR / f"{args.mode}_summary.txt").write_text(summary)
        print(f"\n{summary}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
