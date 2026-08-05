#!/usr/bin/env python3
"""
screen_references.py — screen PHI-base euk pathogen seed taxids for available
transcriptome assemblies, ranked by currency and completeness.

For each seed taxid in phibase_db.json (fungi + oomycetes):
  1. Query NCBI datasets summary for all assemblies under the taxon
  2. Select best assembly: RefSeq+annotation > GenBank+annotation (most recent)
     > RefSeq (no annotation) > GenBank (no annotation, most recent, best level)
  3. Report: accession, release date, assembly level, has_annotation, N50,
     gene count, total size — for manual review before building Kraken2 DB

Output: kraken/ref_screen.tsv
Run:    python scripts/screen_references.py [--workers N]
"""

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

DB_PATH = Path("output/00_build/data/phibase_db.json")
OUT_TSV = Path("kraken/ref_screen.tsv")

LEVEL_RANK = {
    "Complete Genome": 5,
    "Chromosome": 4,
    "Scaffold": 3,
    "Contig": 2,
}


def level_rank(lvl: str) -> int:
    return LEVEL_RANK.get(lvl, 1)


def sort_key(a: dict) -> tuple:
    acc = a.get("accession", "")
    ai = a.get("assembly_info", {})
    ann = a.get("annotation_info") or {}
    stats = a.get("assembly_stats", {})
    genes = ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", 0) or 0
    # Priority: annotation > gene count tier (5k buckets) > assembly level >
    # RefSeq (tie-breaker) > date
    # Gene count outweighs RefSeq: a GenBank with 36k genes beats RefSeq with 15k genes
    return (
        1 if ann else 0,                            # annotation = RNA FASTA available
        min(int(genes), 50000) // 5000,             # gene count tier (capped)
        level_rank(ai.get("assembly_level", "")),   # Chromosome > Scaffold > Contig
        1 if acc.startswith("GCF_") else 0,         # RefSeq as tie-breaker only
        ai.get("release_date", "") or "",           # most recent
    )


def query_taxon(taxid: int, name: str, kingdom: str) -> dict:
    result = subprocess.run(
        ["datasets", "summary", "genome", "taxon", str(taxid),
         "--as-json-lines", "--limit", "all"],
        capture_output=True, text=True, timeout=120
    )

    row = {
        "taxid": taxid,
        "name": name,
        "kingdom": kingdom,
        "best_accession": "",
        "assembly_name": "",
        "assembly_level": "",
        "release_date": "",
        "has_annotation": False,
        "protein_coding_genes": "",
        "scaffold_n50_kb": "",
        "total_length_mb": "",
        "n_assemblies": 0,
        "fasta_type": "",
        "notes": "",
    }

    if result.returncode != 0 or not result.stdout.strip():
        row["notes"] = "NCBI query failed"
        return row

    assemblies = []
    for line in result.stdout.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if "reports" in d:
            assemblies.extend(d["reports"])
        elif "accession" in d:
            assemblies.append(d)

    row["n_assemblies"] = len(assemblies)
    if not assemblies:
        row["notes"] = "no assembly"
        return row

    assemblies.sort(key=sort_key, reverse=True)
    best = assemblies[0]

    acc = best.get("accession", "")
    ai = best.get("assembly_info", {})
    ann = best.get("annotation_info") or {}
    stats = best.get("assembly_stats", {})

    genes = ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", "") if ann else ""
    n50 = stats.get("scaffold_n50") or stats.get("contig_n50", "")
    total = stats.get("total_sequence_length", "")

    row["best_accession"] = acc
    row["assembly_name"] = ai.get("assembly_name", "")
    row["assembly_level"] = ai.get("assembly_level", "")
    row["release_date"] = ai.get("release_date", "")
    row["has_annotation"] = bool(ann)
    row["protein_coding_genes"] = str(genes) if genes else ""
    row["scaffold_n50_kb"] = f"{int(n50)/1000:.0f}" if n50 else ""
    row["total_length_mb"] = f"{int(total)/1e6:.1f}" if total else ""
    row["fasta_type"] = "rna" if ann else "genome"

    date_year = (ai.get("release_date", "") or "")[:4]
    notes = []
    if not ann:
        notes.append("no_annotation→use_genome")
    if date_year and int(date_year) < 2015:
        notes.append(f"OLD_assembly({date_year})")
    if not acc.startswith("GCF_"):
        notes.append("GenBank_only(no_RefSeq)")
    if level_rank(ai.get("assembly_level", "")) <= 2:
        notes.append("contig_level")
    row["notes"] = "; ".join(notes)

    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    db = json.load(open(DB_PATH))
    t2n = db["taxid_to_name"]

    seeds: dict[int, tuple[str, str]] = {}
    for t in set(db["fungal_to_seed"].values()):
        seeds[int(t)] = ("fungi", t2n.get(str(t), str(t)))
    for t in set(db["oomycete_to_seed"].values()):
        seeds[int(t)] = ("oomycete", t2n.get(str(t), str(t)))

    print(f"Screening {len(seeds)} seed taxids "
          f"({sum(1 for k,_ in seeds.values() if k=='fungi')} fungi, "
          f"{sum(1 for k,_ in seeds.values() if k=='oomycete')} oomycetes) ...", flush=True)

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {
            ex.submit(query_taxon, taxid, name, kingdom): taxid
            for taxid, (kingdom, name) in seeds.items()
        }
        done = 0
        for fut in as_completed(futures):
            done += 1
            row = fut.result()
            rows.append(row)
            flag = ("✓rna" if (row["has_annotation"] and row["best_accession"]) else
                    ("✓gen" if row["best_accession"] else "✗"))
            print(f"  [{done:3d}/{len(seeds)}] {flag} {row['taxid']:8d}  "
                  f"{row['name'][:45]:45s}  {row['best_accession']:20s}  "
                  f"{(row['release_date'] or '')[:7]:8s}  "
                  f"genes={row['protein_coding_genes']:6s}  "
                  f"N50={row['scaffold_n50_kb']:6s}kb",
                  flush=True)

    rows.sort(key=lambda r: (r["kingdom"], r["name"].lower()))

    cols = [
        "taxid", "name", "kingdom", "best_accession", "assembly_name",
        "assembly_level", "release_date", "has_annotation", "protein_coding_genes",
        "scaffold_n50_kb", "total_length_mb", "n_assemblies", "fasta_type", "notes",
    ]

    with open(OUT_TSV, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")

    n_rna    = sum(1 for r in rows if r["has_annotation"] and r["best_accession"])
    n_genome = sum(1 for r in rows if not r["has_annotation"] and r["best_accession"])
    n_none   = sum(1 for r in rows if not r["best_accession"])
    n_old    = sum(1 for r in rows if "OLD_assembly" in r["notes"])
    n_no_refseq = sum(1 for r in rows if r["best_accession"] and "GenBank_only" in r["notes"])

    print(f"\n── Summary ({len(rows)} seeds) ──")
    print(f"  Will use RNA FASTA (has annotation):  {n_rna}")
    print(f"  Will use genomic FASTA (no annotation): {n_genome}")
    print(f"  No assembly at all:                   {n_none}")
    print(f"  Old assemblies (<2015):               {n_old}")
    print(f"  No RefSeq (GenBank only):             {n_no_refseq}")
    print(f"\nOutput: {OUT_TSV}")


if __name__ == "__main__":
    main()
