#!/usr/bin/env python3
"""Decide which plant hosts go into the database, and which Ensembl genomes supply their CDS.

    python release/02_reference/scripts/screen_hosts.py
    python release/02_reference/scripts/screen_hosts.py --min-run-frac 0.02

Reads  data/inputs/scope_rule.tsv                   source=host, the criterion
       ../01_cohort/data/outputs/supp_data_s2_cohort.tsv   the cohort's own host column
       data/raw/raw_ensembl_plants_species.txt      what exists to build from
       data/raw/raw_taxdmp.zip                      names.dmp + nodes.dmp, the name authority
Writes data/curation/screen_hosts.tsv               ONE ROW PER COHORT HOST, kept or dropped

A SCREEN OVER THE COHORT, NOT AN AUTHORED LIST. The research tree's `seeds/hosts.tsv` was
hand-written for a four-cereal cohort and went stale the moment the chapter became wheat-only:
it still names Oryza, Zea and Sorghum, none of which has a run left. Deriving the list from S2
means it cannot drift from the cohort it is supposed to serve.

## The rule

    analysed   Triticum aestivum, 1,450 runs. Always in; it is the sink.
    contrast   a cohort host OUTSIDE the analysed genus whose GENUS holds at least
               `contrast_min_run_frac` of the cohort. Hordeum is 61/1,560 = 3.9% and is in.
    outgroup   Arabidopsis thaliana, DECLARED rather than derived. Mislabel QC only.

**ONE SPECIES PER GENUS, and it is the load-bearing clause**, because it drops hosts the
cohort really contains: T. turgidum (8 runs), Triticum sp. (8), T. dicoccon (1), H. jubatum
(4). Adding them would be actively harmful. db_v3 held both hexaploid T. aestivum and
tetraploid T. turgidum, and 27.7% of the Triticum clade's 51.8M minimizers migrated to the
GENUS node, so a correct wheat call scored low coverage; db_v4 fixed it by consolidating.
Re-adding a congeneric host would re-create a defect already paid for. Note the floor is a
property of the GENUS (Hordeum is 61 runs across vulgare and jubatum) while a database entry
is a SPECIES, so the floor qualifies a genus and one representative is picked from it.

**The sink survives every exclusion** because Triticeae CDS are near-identical: a barley, rye
or triticale read lands on wheat or on a Triticeae ancestor, and either way NOT on a fungus,
which is the whole job of the sink.

**The outgroup is not derived from the cohort on purpose.** Host assignment RANKS reads across
hosts rather than thresholding, so with only cohort hosts present every library is called as
one of them by construction and a mislabelled library cannot be seen. A distant dicot gives
that signal somewhere to land. This is how Fragaria's absence was found: 30 strawberry runs
ranking onto Populus, Arachis and Glycine.

## BREADTH SCALES WITH COHORT REPRESENTATION

The analysed host takes EVERY cultivar Ensembl publishes (19 for T. aestivum), because that is
where the reads are and a cultivar-specific transcript missing from the sink leaks onto a
fungus. A contrast or outgroup host takes its REFERENCE assembly only: barley is 57 runs and
Arabidopsis is pure QC, so the 77 barley cultivars Ensembl carries would be weight with no
question behind it.

Taking every cultivar is safe in a way that one-species-per-genus is not, and for the same
reason: all 19 wheat entries share taxid 4565, so they consolidate on ONE node and no
genus-level LCA bleed is possible. Kraken2 stores each distinct minimizer once, so the
redundancy adds no weight to the taxon either.

## TWO TAXONOMY TRAPS, BOTH LIVE HERE

**Ensembl files barley under 112509 (H. vulgare subsp. vulgare), not the cohort's 4513.** So a
host cannot be matched on taxid equality; this walks nodes.dmp to accept any Ensembl genome at
or BELOW the cohort's taxon. The CDS is then tagged with the COHORT's species taxid, so a host
call resolves at the rank the cohort records and the chapter reports. Same reasoning as using
1047171 for Z. tritici rather than the IPO323 strain node 336722.

**Ensembl Genomes numbers its plants releases separately from Ensembl core.** `current` is
plants release-63, which summary.txt calls "Release 116" because it loads into Ensembl v116
databases. `/pub/release-116/` does not exist.
"""

import argparse
import csv
import io
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parent.parent
RULE = HERE / "data/inputs/scope_rule.tsv"
COHORT = REPO / "release/01_cohort/data/outputs/supp_data_s2_cohort.tsv"
ENS = HERE / "data/raw/raw_ensembl_plants_species.txt"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
OUT = HERE / "data/curation/screen_hosts.tsv"

csv.field_size_limit(10 ** 7)


def read_rule():
    out = defaultdict(dict)
    with open(RULE) as fh:
        for row in csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                  delimiter="\t"):
            out[row["source"]][row["field"]] = row["value"]
    return out


def dmp_rows(zf, member):
    with zf.open(member) as raw:
        for line in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
            yield [f.strip() for f in line.split("|")]


def load_taxonomy(keys):
    idx, parent, rank, sci = {}, {}, {}, {}
    with zipfile.ZipFile(TAXDMP) as zf:
        for p in dmp_rows(zf, "names.dmp"):
            if len(p) < 4 or p[3] != "scientific name":
                continue
            tid = int(p[0])
            sci[tid] = p[1]
            if p[1].lower() in keys:
                idx.setdefault(p[1].lower(), tid)
        for p in dmp_rows(zf, "nodes.dmp"):
            if len(p) >= 3:
                parent[int(p[0])], rank[int(p[0])] = int(p[1]), p[2]
    return idx, parent, rank, sci


def genus_of(tid, parent, rank):
    seen = set()
    while tid and tid not in seen:
        if rank.get(tid) == "genus":
            return tid
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt is None or nxt == tid:
            return None
        tid = nxt
    return None


def is_at_or_below(tid, ancestor, parent):
    """Ensembl's taxid may be a SUBSPECIES of the cohort's (barley: 112509 under 4513)."""
    seen = set()
    while tid and tid not in seen:
        if tid == ancestor:
            return True
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt is None or nxt == tid:
            return False
        tid = nxt
    return False


def ensembl_rows():
    with open(ENS, encoding="utf-8", errors="replace") as fh:
        hdr = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [c.lstrip("#") for c in f]
                continue
            if line.strip():
                yield dict(zip(hdr, f))


# THE CDS FILENAME CANNOT BE CONSTRUCTED FROM THIS TABLE, so this screen does not try.
# The species table's `assembly` column and the FASTA's assembly token disagree: Landmark is
# `PGSBv2.0_Landmark` in the table and `PGSBv2.1` in the filename, Weebil is `EI_WeebilV1`
# against `WeebilV1`. Constructing URLs here 404s on 2 of 19 wheat cultivars. So the screen
# names the SPECIES and stays offline, and fetch_host_cds.py lists each cds/ directory to
# resolve the real filename. That split is the module's convention anyway: a screen decides
# from bytes already fetched, a fetch talks to the network.


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--min-run-frac", type=float, default=None,
                    help="override contrast_min_run_frac from scope_rule.tsv")
    args = ap.parse_args()

    rule = read_rule()
    host, ens = rule["host"], rule["ensembl_plants"]
    analysed_tid, outgroup_tid = int(host["analysed"]), int(host["outgroup"])
    frac = args.min_run_frac if args.min_run_frac is not None \
        else float(host["contrast_min_run_frac"])

    with open(COHORT, newline="", encoding="utf-8") as fh:
        cohort = list(csv.DictReader(fh, delimiter="\t"))
    counts = Counter((r["host"] or "").strip() for r in cohort)
    n_runs = len(cohort)
    floor = frac * n_runs

    idx, parent, rank, sci = load_taxonomy({h.lower() for h in counts if h})
    genomes = list(ensembl_rows())

    def available(tid):
        """Ensembl genomes at or below `tid`, reference first then cultivars by name.

        The REFERENCE is the entry whose `species` key carries no cultivar suffix, i.e. the
        shortest name that is a prefix of the others (`hordeum_vulgare` against
        `hordeum_vulgare_barke`). Ensembl has no column marking it.
        """
        hits = [g for g in genomes
                if g["taxonomy_id"].isdigit()
                and is_at_or_below(int(g["taxonomy_id"]), tid, parent)]
        if not hits:
            return []
        ref = min(hits, key=lambda g: len(g["species"]))
        rest = sorted((g for g in hits if g is not ref), key=lambda g: g["species"])
        return [ref] + rest

    analysed_genus = genus_of(analysed_tid, parent, rank)
    per_genus, resolved, best_in_genus = Counter(), {}, {}
    for h, n in counts.items():
        if not h:
            continue
        tid = idx.get(h.lower())
        resolved[h] = tid
        if tid:
            g = genus_of(tid, parent, rank)
            per_genus[g] += n
            cand = (bool(available(tid)), n)
            if g not in best_in_genus or cand > best_in_genus[g][1]:
                best_in_genus[g] = (h, cand)
    best_in_genus = {g: v[0] for g, v in best_in_genus.items()}

    rows = []

    def add(h, n, tid, decision, role, reason):
        got = available(tid) if tid else []
        # Breadth: every cultivar for the analysed host, the reference alone otherwise.
        if role != "analysed":
            got = got[:1]
        rows.append({
            "host": h, "cohort_runs": n, "ncbi_taxid": tid or "",
            "genus": sci.get(genus_of(tid, parent, rank), "") if tid else "",
            "decision": decision, "role": role,
            "n_genomes": len(got),
            "ensembl_species": ";".join(g["species"] for g in got),
            "ensembl_taxids": ";".join(sorted({g["taxonomy_id"] for g in got})),
            "assembly_accessions": ";".join(g["assembly_accession"] for g in got),
            "ensembl_assemblies": ";".join(g["assembly"] for g in got),
            "reason": reason,
        })

    for h, n in counts.most_common():
        if not h:
            add("(blank)", n, None, "drop", "",
                "no host recorded in the cohort; cannot be resolved to a taxid")
            continue
        tid = resolved[h]
        if tid is None:
            add(h, n, None, "drop", "",
                "not a scientific name in names.dmp; cannot be resolved to a taxid")
            continue
        g = genus_of(tid, parent, rank)
        if tid == analysed_tid:
            add(h, n, tid, "keep", "analysed",
                f"the analysed host, {n} cohort runs; EVERY Ensembl cultivar, because this "
                f"is where the reads are")
        elif g == analysed_genus:
            add(h, n, tid, "drop", "",
                "same genus as the analysed host: adding it migrates shared Triticum k-mers "
                "to the genus node and degrades the T. aestivum call (the db_v3 polyploid "
                "defect db_v4 consolidated away). Calls as T. aestivum instead")
        elif per_genus[g] >= floor and best_in_genus.get(g) == h:
            add(h, n, tid, "keep", "contrast",
                f"different genus, and {sci.get(g, g)} holds {per_genus[g]} of {n_runs} "
                f"cohort runs ({100 * per_genus[g] / n_runs:.1f}%), at or above the "
                f"{100 * frac:.0f}% floor; reference assembly only, not the pangenome")
        elif per_genus[g] >= floor:
            add(h, n, tid, "drop", "",
                f"{sci.get(g, g)} qualifies, but {best_in_genus[g]} represents it with more "
                f"runs. ONE SPECIES PER GENUS: two congeneric hosts split the genus node "
                f"between them, the same defect the same-genus exclusion prevents")
        else:
            add(h, n, tid, "drop", "",
                f"{sci.get(g, g)} holds {per_genus[g]} of {n_runs} cohort runs "
                f"({100 * per_genus[g] / n_runs:.1f}%), below the {100 * frac:.0f}% floor; "
                f"lands on a Triticeae neighbour, which is the correct neighbourhood")

    if not any(r["ncbi_taxid"] == outgroup_tid for r in rows):
        add(sci.get(outgroup_tid, str(outgroup_tid)), 0, outgroup_tid, "keep", "outgroup",
            "DECLARED dicot outgroup, not a cohort host. Host assignment ranks rather than "
            "thresholds, so without a distant alternative every library is called a cohort "
            "host by construction and a mislabelled library is undetectable. Reference only")

    for r in rows:
        if r["decision"] == "keep" and not r["ensembl_species"]:
            r["decision"] = "unbuildable"
            r["reason"] = ("selected by the rule but Ensembl Plants publishes no genome at "
                           "or below this taxon, so there is no CDS: " + r["reason"])

    cols = ["host", "cohort_runs", "ncbi_taxid", "genus", "decision", "role", "n_genomes",
            "ensembl_species", "ensembl_taxids", "ensembl_assemblies",
            "assembly_accessions", "reason"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    kept = [r for r in rows if r["decision"] == "keep"]
    print(f"{n_runs:,} cohort runs over {len(counts)} distinct host values; "
          f"contrast floor {100 * frac:.0f}% = {floor:.0f} runs")
    print(f"Ensembl Plants {ens['release']}, {len(genomes)} genomes\n")
    print(f"{'host':<24}{'runs':>6}  {'decision':<12}{'role':<10}{'genomes':>8}")
    for r in rows:
        print(f"{r['host'][:23]:<24}{r['cohort_runs']:>6}  {r['decision']:<12}"
              f"{r['role']:<10}{r['n_genomes']:>8}")
    print(f"\n  {len(kept)} hosts kept, {len(rows) - len(kept)} dropped")
    for r in kept:
        tids = r["ensembl_taxids"]
        note = f"  [Ensembl taxid {tids}, tagged as {r['ncbi_taxid']}]" \
            if tids != str(r["ncbi_taxid"]) else ""
        print(f"    {r['role']:<10}{r['host']:<22}{r['n_genomes']} genome(s){note}")
        for sp in r["ensembl_species"].split(";")[:3]:
            print(f"        {sp}")
        if r["n_genomes"] > 3:
            print(f"        ... and {r['n_genomes'] - 3} more")
    print(f"\n  -> {OUT.relative_to(REPO)}")
    print("Next: fetch_host_cds.py, then plan_build.py")


if __name__ == "__main__":
    main()
