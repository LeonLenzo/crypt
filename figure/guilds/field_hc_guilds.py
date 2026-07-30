#!/usr/bin/env python3
"""
figure/guilds/field_hc_guilds.py — co-infection guild network for high-confidence
field samples only.

Filters runs.tsv to:
  - biosample_representative == True     (one row per biological sample)
  - same_genus_secondary == False        (high-confidence, diff-genus co-infections)
  - BioProject llm_study_setting == field  (LLM-classified field studies)

Run from crypt/:
  python figure/guilds/field_hc_guilds.py
"""

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

RUNS_TSV   = Path("output/02_filter_runs/data/runs.tsv")
LLM_TSV    = Path("output/05_llm_classify/data/bioproject_llm.tsv")
KW_TSV     = Path("output/04_filter_kw/data/biosample_kw.tsv")
DB_PATH    = Path("output/00_build/data/phibase_db.json")

MIN_EDGE_WEIGHT = 1
_SKIP = {"environmental samples"}

# Uninformative broad-clade names excluded when --all-hosts is NOT set.
BROAD_CLADE_NAMES = {
    "Viridiplantae", "Mesangiospermae", "eudicotyledons", "rosids",
    "Euphyllophyta", "Pentapetalae", "Poales", "asterids", "Gunneridae",
    "Magnoliopsida", "Lamiales", "BOP clade", "PACMAD clade",
    "Streptophyta",
}

_FSP_RE = re.compile(r'^(\S+\s+\S+\s+f\.\s+sp\.\s+\S+)')
_PCT_RE = re.compile(r':[\d.]+%$')


def _normalize(name: str) -> str:
    name = _PCT_RE.sub('', name).strip()
    if re.search(r'vir(?:us|oid)', name, re.IGNORECASE):
        return name
    m = _FSP_RE.match(name)
    if m:
        return m.group(1)
    parts = name.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else name


def _build_kingdom_map(db: dict) -> tuple[dict[int, str], dict[str, int]]:
    taxid_to_kg: dict[int, str] = {}
    for key, kg in (("fungal_to_seed",   "Fungi"),
                    ("bacterial_to_seed", "Bacteria"),
                    ("oomycete_to_seed",  "Oomycota"),
                    ("nematode_to_seed",  "Nematoda"),
                    ("virus_to_seed",     "Viruses")):
        for tid in db.get(key, {}):
            taxid_to_kg[int(tid)] = kg
    name_to_tid = {k: int(v) for k, v in db.get("name_to_taxid", {}).items()}
    return taxid_to_kg, name_to_tid


def _kingdom(name: str, taxid_str: str, taxid_to_kg: dict, name_to_tid: dict) -> str:
    if taxid_str:
        try:
            kg = taxid_to_kg.get(int(taxid_str))
            if kg:
                return kg
        except ValueError:
            pass
    tid = name_to_tid.get(name.lower()) or name_to_tid.get(_normalize(name).lower())
    return taxid_to_kg.get(tid, "Unknown") if tid else "Unknown"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-hosts", action="store_true",
                    help="Include broad-clade/unresolved hosts (writes *_all_* files)")
    args = ap.parse_args()

    if args.all_hosts:
        nodes_tsv      = Path("figure/guilds/field_hc_all_nodes.tsv")
        edges_tsv      = Path("figure/guilds/field_hc_all_edges.tsv")
        node_hosts_tsv = Path("figure/guilds/field_hc_all_node_hosts.tsv")
    else:
        nodes_tsv      = Path("figure/guilds/field_hc_nodes.tsv")
        edges_tsv      = Path("figure/guilds/field_hc_edges.tsv")
        node_hosts_tsv = Path("figure/guilds/field_hc_node_hosts.tsv")

    db = json.loads(DB_PATH.read_text())
    taxid_to_kg, name_to_tid = _build_kingdom_map(db)

    # Field BPs + treatment from LLM classifications
    field_bps:    set[str] = set()
    bp_treatment: dict[str, str] = {}
    with open(LLM_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("llm_study_setting") == "field":
                field_bps.add(r["BioProject"])
            bp_treatment[r["BioProject"]] = r.get("llm_treatment", "unclear")
    print(f"Field BPs (LLM): {len(field_bps)}")

    # named_host fallback from biosample_kw (first non-empty per BioProject)
    bp_named_host: dict[str, str] = {}
    with open(KW_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            bp = r["BioProject"]
            if bp not in bp_named_host and r.get("named_host"):
                bp_named_host[bp] = r["named_host"]
    print(f"BPs with named_host: {len(bp_named_host)}")

    node_kingdom:      dict[str, str] = {}
    node_n_primary:    dict[str, int] = defaultdict(int)
    node_n_secondary:  dict[str, int] = defaultdict(int)
    node_hosts:        dict[str, Counter] = defaultdict(Counter)
    edge_counts:       dict[tuple[str, str], int] = defaultdict(int)
    edge_treatments:   dict[tuple[str, str], Counter] = defaultdict(Counter)
    n_rows = 0
    n_host_excluded = 0

    with open(RUNS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("biosample_representative") != "True":
                continue
            if row.get("BioProject") not in field_bps:
                continue
            if row.get("co_infection_flag") == "single":
                continue
            if not row.get("stat_pathogens"):
                continue
            host = row.get("host", "").strip()
            if not args.all_hosts and (host in BROAD_CLADE_NAMES or not host):
                fallback = bp_named_host.get(row["BioProject"], "")
                if fallback:
                    host = fallback
                else:
                    n_host_excluded += 1
                    continue
            n_rows += 1

            treat = bp_treatment.get(row["BioProject"], "unclear")

            prim_norm  = _normalize(row["library_organism"])
            prim_genus = prim_norm.split()[0].lower()
            prim_kg    = _kingdom(row["library_organism"], "", taxid_to_kg, name_to_tid)
            node_kingdom.setdefault(prim_norm, prim_kg)
            node_n_primary[prim_norm] += 1
            if host:
                node_hosts[prim_norm][host] += 1

            sec_norms: list[str] = []
            for entry in row["stat_pathogens"].split("; "):
                raw = _PCT_RE.sub("", entry).strip()
                if raw.lower() in _SKIP or not raw:
                    continue
                norm = _normalize(raw)
                # Skip same-genus secondaries (low k-mer specificity within genus).
                # Diff-genus pairs are retained even when same-genus secondaries
                # co-occur in the same run.
                if norm.split()[0].lower() == prim_genus:
                    continue
                kg   = _kingdom(raw, "", taxid_to_kg, name_to_tid)
                node_kingdom.setdefault(norm, kg)
                node_n_secondary[norm] += 1
                if host:
                    node_hosts[norm][host] += 1
                sec_norms.append(norm)

            run_pathogens = sorted({prim_norm} | set(sec_norms))
            for i, a in enumerate(run_pathogens):
                for b in run_pathogens[i + 1:]:
                    edge_counts[(a, b)] += 1
                    edge_treatments[(a, b)][treat] += 1

    with open(nodes_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["name", "kingdom", "n_primary", "n_secondary", "total", "top_host"])
        for name in sorted(node_kingdom):
            n_p = node_n_primary.get(name, 0)
            n_s = node_n_secondary.get(name, 0)
            top_h = node_hosts[name].most_common(1)[0][0] if node_hosts.get(name) else ""
            w.writerow([name, node_kingdom[name], n_p, n_s, n_p + n_s, top_h])

    n_edges_kept = 0
    with open(edges_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["node1", "node2", "weight", "dominant_treatment", "treatment_breakdown"])
        for (a, b), wt in sorted(edge_counts.items(), key=lambda x: -x[1]):
            if wt >= MIN_EDGE_WEIGHT:
                ct = edge_treatments[(a, b)]
                dominant = ct.most_common(1)[0][0]
                breakdown = "|".join(f"{t}:{n}" for t, n in ct.most_common())
                w.writerow([a, b, wt, dominant, breakdown])
                n_edges_kept += 1

    with open(node_hosts_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["name", "host", "n"])
        for name in sorted(node_hosts):
            for host, n in node_hosts[name].most_common():
                w.writerow([name, host, n])

    n_single = sum(1 for wt in edge_counts.values() if wt == 1)
    print(f"Runs processed: {n_rows}  (broad-clade excluded: {n_host_excluded})")
    print(f"Nodes: {len(node_kingdom)}")
    print(f"Edges (weight >= {MIN_EDGE_WEIGHT}): {n_edges_kept}  "
          f"(singleton edges dropped: {n_single})")
    print(f"Written: {nodes_tsv}, {edges_tsv}, {node_hosts_tsv}")


if __name__ == "__main__":
    main()
