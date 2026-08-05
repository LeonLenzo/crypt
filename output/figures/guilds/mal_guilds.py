#!/usr/bin/env python3
"""
figure/guilds/mal_guilds.py — build co-infection guild network from runs.tsv (MAL + HAL).

Treats all pathogens (primary + secondary) symmetrically. Nodes = pathogens,
edges = co-occurrence in the same confirmed run, edge weight = run count.

Normalises strain-level names to species (or f. sp.) level.
Kingdom is resolved via the phibase_db name_to_taxid lookup.

Outputs (figure/guilds/):
  mal_guild_nodes.tsv    name, kingdom, n_primary, n_secondary, total
  mal_guild_edges.tsv    node1, node2, weight

Run from crypt/:
  python figure/guilds/mal_guilds.py
"""

import csv
import json
import re
from collections import defaultdict
from pathlib import Path

RUNS_TSV        = Path("output/02_filter_runs/data/runs.tsv")
DB_PATH         = Path("output/00_build/data/phibase_db.json")
NODES_TSV       = Path("output/figures/guilds/mal_guild_nodes.tsv")
EDGES_TSV       = Path("output/figures/guilds/mal_guild_edges.tsv")

MIN_EDGE_WEIGHT = 2    # drop edges seen in fewer than this many runs
_SKIP           = {"environmental samples"}


# ── Name normalisation ────────────────────────────────────────────────────────

_FSP_RE    = re.compile(r'^(\S+\s+\S+\s+f\.\s+sp\.\s+\S+)')
_PCT_RE    = re.compile(r':[\d.]+%$')


def _normalize(name: str) -> str:
    """Collapse strain identifiers; keep virus common names and f. sp. epithets intact."""
    name = _PCT_RE.sub('', name).strip()
    # Virus/viroid common names are multi-word; don't truncate them
    if re.search(r'vir(?:us|oid)', name, re.IGNORECASE):
        return name
    # Keep f. sp. designation
    m = _FSP_RE.match(name)
    if m:
        return m.group(1)
    # Truncate strain-appended binomials to genus + species epithet only
    parts = name.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else name


# ── DB kingdom lookup ─────────────────────────────────────────────────────────

def _build_kingdom_map(db: dict) -> tuple[dict[int, str], dict[str, int]]:
    """Return (taxid→kingdom, lowercase_name→taxid) from phibase_db."""
    taxid_to_kg: dict[int, str] = {}
    for key, kg in (("fungal_to_seed",    "Fungi"),
                    ("bacterial_to_seed",  "Bacteria"),
                    ("oomycete_to_seed",   "Oomycota"),
                    ("nematode_to_seed",   "Nematoda"),
                    ("virus_to_seed",      "Viruses")):
        for tid in db.get(key, {}):
            taxid_to_kg[int(tid)] = kg
    name_to_tid = {k: int(v) for k, v in db.get("name_to_taxid", {}).items()}
    return taxid_to_kg, name_to_tid


def _kingdom(name: str, taxid_str: str,
             taxid_to_kg: dict, name_to_tid: dict) -> str:
    """Resolve kingdom for a pathogen name/taxid pair."""
    if taxid_str:
        try:
            kg = taxid_to_kg.get(int(taxid_str))
            if kg:
                return kg
        except ValueError:
            pass
    tid = name_to_tid.get(name.lower()) or name_to_tid.get(_normalize(name).lower())
    return taxid_to_kg.get(tid, "Unknown") if tid else "Unknown"


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    with open(DB_PATH) as f:
        db = json.load(f)
    taxid_to_kg, name_to_tid = _build_kingdom_map(db)

    node_kingdom:     dict[str, str] = {}
    node_n_primary:   dict[str, int] = defaultdict(int)
    node_n_secondary: dict[str, int] = defaultdict(int)
    edge_counts:      dict[tuple[str, str], int] = defaultdict(int)

    with open(RUNS_TSV) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if not row.get('stat_pathogens'):
                continue

            # Primary node
            prim_norm = _normalize(row['library_organism'])
            prim_kg   = _kingdom(row['library_organism'], '',
                                 taxid_to_kg, name_to_tid)
            node_kingdom.setdefault(prim_norm, prim_kg)
            node_n_primary[prim_norm] += 1

            # Secondary nodes
            sec_norms: list[str] = []
            for entry in row['stat_pathogens'].split('; '):
                raw = _PCT_RE.sub('', entry).strip()
                if raw.lower() in _SKIP or not raw:
                    continue
                norm = _normalize(raw)
                kg   = _kingdom(raw, '', taxid_to_kg, name_to_tid)
                node_kingdom.setdefault(norm, kg)
                node_n_secondary[norm] += 1
                sec_norms.append(norm)

            # Co-occurrence pairs for this run
            run_pathogens = sorted({prim_norm} | set(sec_norms))
            for i, a in enumerate(run_pathogens):
                for b in run_pathogens[i + 1:]:
                    edge_counts[(a, b)] += 1

    # Write nodes
    with open(NODES_TSV, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['name', 'kingdom', 'n_primary', 'n_secondary', 'total'])
        for name in sorted(node_kingdom):
            n_p = node_n_primary.get(name, 0)
            n_s = node_n_secondary.get(name, 0)
            w.writerow([name, node_kingdom[name], n_p, n_s, n_p + n_s])

    # Write edges (filtered)
    with open(EDGES_TSV, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['node1', 'node2', 'weight'])
        for (a, b), wt in sorted(edge_counts.items(), key=lambda x: -x[1]):
            if wt >= MIN_EDGE_WEIGHT:
                w.writerow([a, b, wt])

    n_edges = sum(1 for wt in edge_counts.values() if wt >= MIN_EDGE_WEIGHT)
    print(f"Nodes: {len(node_kingdom)}")
    print(f"Edges (weight >= {MIN_EDGE_WEIGHT}): {n_edges}")
    print(f"Written: {NODES_TSV}, {EDGES_TSV}")


if __name__ == '__main__':
    main()
