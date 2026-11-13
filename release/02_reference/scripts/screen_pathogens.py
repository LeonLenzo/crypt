#!/usr/bin/env python3
"""Label which fungal species are known plant / wheat pathogens, from PHI-base.

    python release/02_reference/scripts/screen_pathogens.py

Reads  data/raw/raw_phibase.csv           PHI-base v4-19, curated pathogen-host interactions
       data/raw/raw_taxdmp.zip            nodes.dmp + names.dmp, for host kingdom and genus
       data/curation/screen_species.tsv   the USDA master, to flag which pathogens are in scope
       data/curation/screen_merges.tsv    provisional NCBI names folded onto their binomial
Writes data/curation/screen_pathogens.tsv one row per PHI-base FUNGAL pathogen species

A SEPARATE table on purpose. The pathogen label is not part of the reference set and never
enters the Kraken2 database (which is pure taxonomy); it is an annotation applied to a
detection afterwards. Keeping it here, keyed by species taxid, also gives one clean row per
pathogen instead of the label inheriting across every USDA name variant of a species in
screen_species.tsv.

PHI-base ships the NCBI taxid of BOTH pathogen (`Pathogen_NCBI_species_Taxonomy ID`) and host
(`Host_NCBI_Taxonomy_ID`), so every column here is a JOIN, never a name or substring match.
That is the whole reason to prefer it over FUNGuild/FungalTraits, which assign plant-pathogen
at genus by name substring and labelled chytrids, pollen parasites and the industrial yeast
*Komagataella phaffii* as plant pathogens in db_v4.

    is_plant_pathogen  >=1 PHI-base interaction on a Viridiplantae host THAT NAMES A DISEASE
    is_wheat_pathogen  >=1 PHI-base interaction on a Triticum host that names a disease
    in_usda_genus      the species' genus is in the USDA master, so it is reachable by scope.
                       The test the design must pass: every wheat pathogen should be yes.

The host condition is load-bearing: PHI-base scores pathogenicity against whatever host the
experiment used, and ~200 of its fungi are animal pathogens. The mutant-phenotype columns
describe a knockout, not whether the wild type is a pathogen, so presence of the interaction
suffices for THAT host.

THE DISEASE CONDITION IS THE SECOND HALF OF IT, and the host condition alone is not enough.
A PHI-base row records whatever host the experiment used, which is not always the host being
harmed: in a tripartite insect-fungus-plant assay, or an endophyte colonisation study, the
plant is present but is not the victim. Those rows carry `Disease_name = "no data found"`,
so requiring at least one NAMED disease separates them. Five fungi in this reference were
labelled plant pathogens on disease-less rows alone, and every one is a plant ASSOCIATE
rather than a pathogen: *Epichloe festucae* (the grass endophyte), *Aspergillus ochraceus*
(storage/ochratoxin), *Leptographium clavigerum* (beetle-vectored blue stain) and
*Metarhizium acridum* / *M. robertsii* (insect pathogens, and M. robertsii a root endophyte).
This is the same error the FUNGuild note above indicts, reached by a different route.

The condition is applied to the wheat label too, for consistency. It is a NO-OP there: all
22 PHI-base wheat species already carry at least one named disease, so the wheat count does
not move. Nothing is dropped from the table either way -- a species whose only plant rows
lack a disease stays, with `plant_pathogen_evidence = host_only` recording why it was not
labelled.

Labels are at SPECIES rank because PHI-base keys pathogens by species taxid. A detection
resolves only to species, so this is the right granularity; the one thing to read correctly is
that a species-complex label (e.g. every *Fusarium oxysporum* forma specialis) reflects the
species, not the individual f. sp.
"""

import csv
import io
import zipfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
PHIBASE = HERE / "data/raw/raw_phibase.csv"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
SPECIES = HERE / "data/curation/screen_species.tsv"
MERGES = HERE / "data/curation/screen_merges.tsv"
CURATED = HERE / "data/inputs/wheat_pathogens_curated.tsv"
OUT = HERE / "data/curation/screen_pathogens.tsv"
FUNGI, VIRIDIPLANTAE, TRITICUM = 4751, 33090, 4564  # 4564 = Triticum GENUS (was
# 4565 = T. aestivum, which missed PHI's 1,194 genus-level and 51 durum host records)

csv.field_size_limit(10 ** 7)
# PHI-base writes an absent disease as this literal; the rest are defensive
NO_DISEASE = {"no data found", "", "na", "n/a", "none", "unknown"}


def load_taxonomy():
    parent, rank, name = {}, {}, {}
    with zipfile.ZipFile(TAXDMP) as zf:
        with zf.open("nodes.dmp") as raw:
            for ln in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
                p = [f.strip() for f in ln.split("|")]
                parent[int(p[0])], rank[int(p[0])] = int(p[1]), p[2]
        with zf.open("names.dmp") as raw:
            for ln in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
                p = [f.strip() for f in ln.split("|")]
                if len(p) >= 4 and p[3] == "scientific name":
                    name[int(p[0])] = p[1]
    return parent, rank, name


def under(tid, root, parent):
    seen = set()
    while tid and tid not in seen:
        if tid == root:
            return True
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt == tid or nxt is None:
            return False
        tid = nxt
    return False


def ancestor_at(tid, want, parent, rank):
    seen = set()
    while tid and tid not in seen:
        if rank.get(tid) == want:
            return tid
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt == tid or nxt is None:
            return None
        tid = nxt
    return None


def species_of(tid, parent, rank):
    seen = set()
    while tid and tid not in seen:
        if rank.get(tid) == "species":
            return tid
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt == tid or nxt is None:
            return tid
        tid = nxt
    return tid


def load_merges():
    """-> {from_taxid: to_taxid}, so a label follows its genome onto the merged taxon."""
    merge = {}
    if not MERGES.exists():
        return merge
    for r in csv.DictReader((ln for ln in open(MERGES, encoding="utf-8")
                             if not ln.startswith("#")), delimiter="\t"):
        if r["action"] == "merge" and r["to_taxid"].isdigit():
            merge[int(r["from_taxid"])] = int(r["to_taxid"])
    return merge


def usda_genera(parent, rank):
    """The genus set the USDA master reaches, from screen_species.tsv's resolved taxids."""
    genera = set()
    if not SPECIES.exists():
        print(f"  note: {SPECIES.name} not found; in_usda_genus will be blank. "
              f"Run screen_species.py first.")
        return None
    with open(SPECIES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ncbi_taxid"]:
                g = ancestor_at(int(r["ncbi_taxid"]), "genus", parent, rank)
                if g:
                    genera.add(g)
    return genera


def main():
    parent, rank, name = load_taxonomy()
    genera = usda_genera(parent, rank)
    merge = load_merges()

    # One pass over PHI-base, aggregating to the pathogen species node.
    plant_hosts = defaultdict(set)    # species -> {plant host species}
    plant_named = defaultdict(int)    # plant rows that name a disease
    wheat_rows, wheat_named = defaultdict(int), defaultdict(int)
    n_rows = defaultdict(int)
    with open(PHIBASE, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh)
        pcol = next(c for c in rd.fieldnames if "Pathogen_NCBI_species" in c)
        hcol = next(c for c in rd.fieldnames if "Host_NCBI_Taxonomy_ID" in c)
        dcol = next(c for c in rd.fieldnames if c.strip().lower() == "disease_name")
        for r in rd:
            pt = (r[pcol] or "").strip()
            if not pt.isdigit():
                continue
            pt = int(pt)
            if pt not in rank or not under(pt, FUNGI, parent):
                continue
            sp = species_of(pt, parent, rank)
            sp = merge.get(sp, sp)
            n_rows[sp] += 1
            ht = (r[hcol] or "").strip()
            if not ht.isdigit():
                continue
            ht = int(ht)
            named = (r.get(dcol) or "").strip().lower() not in NO_DISEASE
            if ht in rank and under(ht, VIRIDIPLANTAE, parent):
                plant_hosts[sp].add(species_of(ht, parent, rank))
                plant_named[sp] += named
            if ht in rank and under(ht, TRITICUM, parent):
                wheat_rows[sp] += 1
                wheat_named[sp] += named

    # Curated wheat-pathogen supplement: canonical wheat diseases PHI-base does not curate.
    # Unioned into the labels, with wheat_pathogen_source tagging every call's provenance.
    curated = {}
    if CURATED.exists():
        for r in csv.DictReader((ln for ln in open(CURATED, encoding="utf-8")
                                 if not ln.startswith("#")), delimiter="\t"):
            if r.get("ncbi_taxid", "").strip().isdigit():
                curated[int(r["ncbi_taxid"])] = r.get("disease", "")

    rows = []
    for sp in sorted(set(n_rows) | set(curated), key=lambda t: name.get(t, str(t))):
        phi_plant = plant_named[sp] > 0          # a named disease on a plant host
        phi_wheat = wheat_named[sp] > 0
        cur_wheat = sp in curated
        is_plant = phi_plant or cur_wheat        # a curated wheat pathogen is a plant pathogen
        is_wheat = phi_wheat or cur_wheat
        evidence = ("named_disease" if phi_plant else
                    "curated" if cur_wheat else
                    "host_only" if sp in plant_hosts else "no_plant_host")
        src = "+".join(x for x in ("phibase" if phi_wheat else "",
                                   "curated" if cur_wheat else "") if x)
        g = ancestor_at(sp, "genus", parent, rank)
        in_usda = "" if genera is None else ("yes" if g in genera else "no")
        rows.append({
            "ncbi_taxid": sp, "species": name.get(sp, ""),
            "is_plant_pathogen": "yes" if is_plant else "no",
            "is_wheat_pathogen": "yes" if is_wheat else "no",
            "wheat_pathogen_source": src,
            "plant_pathogen_evidence": evidence,
            "n_plant_hosts": len(plant_hosts.get(sp, ())),
            "n_plant_named_disease": plant_named.get(sp, 0),
            "n_wheat_records": wheat_rows.get(sp, 0),
            "n_wheat_named_disease": wheat_named.get(sp, 0),
            "n_phibase_records": n_rows.get(sp, 0),
            "in_usda_genus": in_usda,
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["ncbi_taxid", "species", "is_plant_pathogen", "is_wheat_pathogen",
            "wheat_pathogen_source", "plant_pathogen_evidence", "n_plant_hosts",
            "n_plant_named_disease", "n_wheat_records", "n_wheat_named_disease",
            "n_phibase_records", "in_usda_genus"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    plant = [r for r in rows if r["is_plant_pathogen"] == "yes"]
    wheat = [r for r in rows if r["is_wheat_pathogen"] == "yes"]
    host_only = [r for r in rows if r["plant_pathogen_evidence"] == "host_only"]
    print(f"PHI-base fungal pathogen taxa: {len(rows):,}")
    print(f"  on a plant host, named disease   {len(plant):,}")
    print(f"  on wheat, named disease          {len(wheat):,}")
    print(f"  plant host but NO named disease  {len(host_only):,} (not labelled): "
          + ", ".join(r["species"] for r in host_only))
    if genera is not None:
        miss = [r for r in wheat if r["in_usda_genus"] == "no"]
        print(f"  wheat pathogens whose genus is NOT in the USDA master: {len(miss):,}"
              + (" — " + ", ".join(r["species"] for r in miss) if miss else ""))
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
