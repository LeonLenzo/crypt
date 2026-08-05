#!/usr/bin/env python3
"""
Check which PHI-base eukaryotic pathogen seed species have RefSeq assemblies
with RNA annotation available (needed for kallisto index build).

Run from crypt/: python scripts/check_refseq_coverage.py
"""
import json
import subprocess
import sys
from pathlib import Path

DB = Path("output/00_build/data/phibase_db.json")
OUT = Path("scripts/refseq_coverage.tsv")

db  = json.load(open(DB))
t2n = db["taxid_to_name"]

# eukaryotic kingdoms only (bacteria excluded — polyA pipeline)
EUK_KINGDOMS = {
    "fungi":    db["fungal_to_seed"],
    "oomycete": db["oomycete_to_seed"],
    "nematode": db["nematode_to_seed"],
}

# collect unique seed taxids per kingdom
seeds = {}   # taxid -> kingdom
for kingdom, d in EUK_KINGDOMS.items():
    for taxid in set(d.values()):
        seeds[taxid] = kingdom

print(f"Checking {len(seeds)} eukaryotic seed taxids against NCBI RefSeq ...")
print(f"{'Taxid':<10} {'Kingdom':<10} {'RefSeq':<8} {'Annotated':<10} {'Accession':<20} Name")
print("-" * 100)

rows = []
n_refseq = n_annotated = 0

for i, (taxid, kingdom) in enumerate(sorted(seeds.items()), 1):
    name = t2n.get(str(taxid), str(taxid))
    print(f"  [{i:3d}/{len(seeds)}] {taxid} {name[:50]}", end="\r", flush=True)

    result = subprocess.run(
        ["datasets", "summary", "genome", "taxon", str(taxid)],
        capture_output=True, text=True
    )

    has_refseq = False
    has_annotation = False
    accession = ""
    asm_name = ""

    if result.returncode == 0 and result.stdout.strip():
        try:
            data = json.loads(result.stdout)
            reports = data.get("reports", [])
            if reports:
                # prefer chromosome/complete with annotation; take first otherwise
                best = reports[0]
                for r in reports:
                    if r.get("annotation_info"):
                        best = r
                        break
                has_refseq = True
                accession = best.get("accession", "")
                asm_name  = best.get("assembly_info", {}).get("assembly_name", "")
                has_annotation = bool(best.get("annotation_info"))
                n_refseq += 1
                if has_annotation:
                    n_annotated += 1
        except json.JSONDecodeError:
            pass

    rows.append({
        "taxid": taxid, "kingdom": kingdom, "name": name,
        "has_refseq": has_refseq, "has_annotation": has_annotation,
        "accession": accession, "asm_name": asm_name,
    })
    status = ("✓ann" if has_annotation else "✓asm" if has_refseq else "✗")
    print(f"  {taxid:<10} {kingdom:<10} {str(has_refseq):<8} {str(has_annotation):<10} {accession:<20} {name[:40]}  [{status}]")

print()
print(f"Summary: {len(seeds)} euk seed species | RefSeq: {n_refseq} ({100*n_refseq/len(seeds):.0f}%) | Annotated: {n_annotated} ({100*n_annotated/len(seeds):.0f}%)")

# missing species
missing = [r for r in rows if not r["has_refseq"]]
no_annot = [r for r in rows if r["has_refseq"] and not r["has_annotation"]]
print(f"\nNo RefSeq assembly ({len(missing)}):")
for r in missing:
    print(f"  {r['taxid']:<10} {r['kingdom']:<10} {r['name']}")
print(f"\nRefSeq but no annotation ({len(no_annot)}):")
for r in no_annot:
    print(f"  {r['taxid']:<10} {r['kingdom']:<10} {r['name']}  [{r['accession']}]")

# write TSV
with open(OUT, "w") as f:
    f.write("taxid\tkingdom\tname\thas_refseq\thas_annotation\taccession\tasm_name\n")
    for r in rows:
        f.write(f"{r['taxid']}\t{r['kingdom']}\t{r['name']}\t"
                f"{r['has_refseq']}\t{r['has_annotation']}\t"
                f"{r['accession']}\t{r['asm_name']}\n")
print(f"\nWritten: {OUT}")
