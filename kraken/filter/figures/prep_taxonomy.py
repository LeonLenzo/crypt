#!/usr/bin/env python3
"""Lineage tables for the composite grid's marginal trees.

Pathogen taxids come from the PHI-base reference (canonical species identity, preferred
over db_v3 names); host taxids from the read-based host calls. Lineages are walked from
the NCBI taxonomy dump. Output: two small rank tables the R side turns into trees.
"""
import json, csv
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
DUMP = ROOT.parents[1] / "_refdata/databases/ncbi_taxonomy/2026-09-17"
RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

# --- target taxids ---
phi = json.load(open(ROOT / "stat/build/data/phibase_db.json"))
n2t = phi["name_to_taxid"]
cl = list(csv.DictReader(open(ROOT / "kraken/filter/data/_classified.tsv"), delimiter="\t"))
classifiable = {r["species"] for r in cl if r["cls"] != "unclassified"}
patho = {n2t[s.lower()]: s for s in classifiable if s.lower() in n2t}     # taxid -> name

host_tid = {}
for r in csv.DictReader(open(ROOT / "kraken/assign/host/data/host_calls.tsv"), delimiter="\t"):
    t = (r.get("host_taxid") or "").strip()
    if t and (r.get("host") or "").strip():
        host_tid[int(t)] = r["host"].strip()

# --- taxonomy dump: parent + rank, then scientific names for the nodes we touch ---
parent, rank = {}, {}
for line in open(DUMP / "nodes.dmp"):
    f = [x.strip() for x in line.split("|")]
    parent[int(f[0])] = int(f[1]); rank[int(f[0])] = f[2]

def lineage(tid):
    out, seen = {}, set()
    while tid and tid not in seen:
        seen.add(tid)
        if rank.get(tid) in RANKS:
            out[rank[tid]] = tid
        if tid == 1: break
        tid = parent.get(tid, 1) if tid != 1 else 0
    return out

targets = list(patho) + list(host_tid)
lins = {t: lineage(t) for t in targets}
need = {tid for lin in lins.values() for tid in lin.values()} | set(targets)

name = {}
for line in open(DUMP / "names.dmp"):
    f = [x.strip() for x in line.split("|")]
    if f[3] == "scientific name" and int(f[0]) in need:
        name[int(f[0])] = f[1]

def write(path, tid_label):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["taxid", "label"] + RANKS)
        for tid, lbl in sorted(tid_label.items(), key=lambda kv: kv[1]):
            lin = lins[tid]
            row = [tid, lbl]
            last = "root"
            for rk in RANKS:                      # forward-fill missing ranks with the last seen name
                if rk in lin:
                    last = name.get(lin[rk], f"tx{lin[rk]}")
                row.append(last)
            w.writerow(row)
    print(f"  wrote {path.name}: {len(tid_label)} taxa")

write(ROOT / "kraken/filter/data/patho_lineage.tsv", patho)
write(ROOT / "kraken/filter/data/host_lineage.tsv", host_tid)
