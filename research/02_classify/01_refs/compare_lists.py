#!/usr/bin/env python3
"""compare_lists.py — what each scope source uniquely contributes, measured against results.

    python 02_classify/01_refs/compare_lists.py

Produces the numbers behind `release/02_reference/README.md`. Run it before changing the
scope rule, because every estimate this project has made from availability counts has been
wrong and the only honest test of a membership list is what it contributed to the detections.

Three questions, in order, because the answers are not interchangeable:

  1 SET RELATIONS between the scope sources, at three projections. The sources sit at
    different ranks and a naive species comparison asks the genus-keyed ones for something
    they cannot give: FungalTraits is 1,069 genera and 0 species, residents is 30 genera and
    0 species, FUNGuild is a 1,500/1,062 genus/species mix. So: as given, at species, and at
    genus, plus a "buildable" projection that expands a genus member to the assembly-bearing
    species inside it, which is the only projection that says what a DATABASE would hold.

  2 COVERAGE OF THE RESULTS. For every fungal detection in gold/detections.tsv, is its
    species in USDA, is only its genus in USDA, or is it outside entirely. This is the test
    that decided the scope rule. Note what it does NOT prove: a detection can only exist for
    a taxon that was in db_v4, so this measures whether USDA covers what the UNION database
    found, not whether a USDA-only database would have found the same. The union is the
    generous case.

  3 WHAT WOULD BE LOST. The built taxa a genus-level USDA master drops, the reads they
    absorbed, the pathogen labels they carry and where those labels came from, and whether
    any DECLARED pathogen becomes unreachable. The declared test is the one that would break
    the chapter, so it is run against 01_literature's curated names, not against a guild.

Every taxid resolves through `_taxonomy`, so a forma specialis is promoted to its species
node and a synonym collapses onto the accepted one before any comparison.
"""

import argparse
import csv
import itertools
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _layout import (HOLOBASE, SILVER, GOLD, RESIDENTS, LIT_RUNS)   # noqa: E402
from _seeds import read_seed_rows                                   # noqa: E402
import _taxonomy                                                    # noqa: E402

SOURCES = ("funguild", "phibase", "fungaltraits", "usda", "residents")


def load_sources(tax, con):
    """{source name: set of internal tids}, fungi only, as each source gives them."""
    fungi = lambda t: t in tax.rank and tax.kingdom.get(t) == "Fungi"
    q = lambda sql: {t for (t,) in con.execute(sql) if fungi(t)}
    L = {}
    L["funguild"] = q("SELECT DISTINCT taxon_id FROM taxon_guild "
                      "WHERE guild LIKE '%Plant Pathogen%'")
    # The plant-host condition is load-bearing: PHI-base scores pathogenicity against
    # whatever host the experiment used, and 147 of its fungi are pathogens of Metazoa.
    L["phibase"] = q("SELECT DISTINCT a.fungus_id FROM host_association a "
                     "JOIN taxon h ON h.taxon_id = a.host_id "
                     "WHERE a.association_type='pathogenic' AND a.fungus_id IS NOT NULL "
                     "AND h.kingdom='Viridiplantae'")
    L["fungaltraits"] = q("SELECT DISTINCT taxon_id FROM taxon_trait "
                          "WHERE trait='primary_lifestyle' AND value='plant_pathogen'")
    L["usda"] = q("SELECT DISTINCT a.fungus_id FROM host_association a "
                  "JOIN taxon h ON h.taxon_id = a.host_id "
                  "WHERE a.source='usda_fungus_host' AND a.association_type='recorded_on' "
                  "AND a.fungus_id IS NOT NULL AND h.kingdom='Viridiplantae'")
    L["residents"] = {t for t in (tax.resolve(p[0].strip())
                                  for p in read_seed_rows(RESIDENTS, min_cols=3)
                                  if len(p) >= 3 and p[2].strip() in ("add", "miniprot"))
                      if t and fungi(t)}
    return L


def projections(tax, con, L):
    """The four views. `buildable` is the one that describes a database."""
    asm = {t for (t,) in con.execute("SELECT DISTINCT taxon_id FROM assembly")}
    asm_sp = {x for x in (tax.species_of(t) for t in asm) if x}

    def species(s):
        return {x for x in (tax.species_of(t) for t in s) if x}

    def genus(s):
        return {x for x in (tax.ancestor_at(t, "genus") for t in s) if x}

    def buildable(s):
        out = set()
        for t in s:
            if tax.rank.get(t) == "genus":
                out |= {d for d in tax.descendants(t)
                        if tax.rank.get(d) == "species" and d in asm_sp}
            else:
                x = tax.species_of(t)
                if x in asm_sp:
                    out.add(x)
        return out

    return (("AS GIVEN", lambda s: set(s)), ("SPECIES RANK", species),
            ("GENUS RANK", genus), ("BUILDABLE SPECIES", buildable))


def report_sets(tax, L):
    print("\nRANK CENSUS, fungi only")
    below = {"forma specialis", "varietas", "forma", "subspecies", "strain", "isolate"}
    print(f"  {'source':<16}{'total':>8}{'species':>9}{'genus':>8}{'below sp':>10}")
    for k in SOURCES:
        c = defaultdict(int)
        for t in L[k]:
            r = tax.rank.get(t) or "?"
            c["sp" if r == "species" else "gn" if r == "genus"
              else "bl" if r in below else "other"] += 1
        print(f"  {k:<16}{len(L[k]):>8,}{c['sp']:>9,}{c['gn']:>8,}{c['bl']:>10,}")


def report_projection(title, S):
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")
    print(f"  union of all five: {len(set().union(*S.values())):,}")
    for k in SOURCES:
        own = S[k] - set().union(*[S[j] for j in SOURCES if j != k])
        print(f"  {k:<14}{len(S[k]):>8,}   unique to it {len(own):>7,}")
    print(f"\n  |ROW \\ COL|" + "".join(f"{k[:12]:>14}" for k in SOURCES))
    for a in SOURCES:
        print(f"  {a:<12}" + "".join("             ." if a == b else f"{len(S[a] - S[b]):>14,}"
                                     for b in SOURCES))
    sup = [f"{a} > {b}" for a, b in itertools.permutations(SOURCES, 2)
           if S[b] and S[b] <= S[a]]
    print("  supersets: " + (", ".join(sup) if sup else "none"))


def usda_status(tax, U_SP, U_GN, ncbi):
    """Where a detected taxon sits relative to a GENUS-level USDA master."""
    t = tax.by_ncbi.get(ncbi)
    if not t:
        return "unresolved", None
    if tax.species_of(t) in U_SP:
        return "usda_species", t
    gn = tax.ancestor_at(t, "genus")
    if gn and gn in U_GN:
        return "usda_genus_only", t
    return "outside", t


def report_coverage(tax, U_SP, U_GN):
    n = defaultdict(int)
    lost = defaultdict(int)
    gonly_path = defaultdict(int)
    with open(GOLD / "detections.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kingdom"] != "Fungi":
                n["host"] += 1
                continue
            s, _ = usda_status(tax, U_SP, U_GN, int(r["taxid"]))
            n[s] += 1
            if s == "usda_genus_only" and r["pathogen"] == "yes":
                gonly_path[r["species"]] += 1
            if s in ("outside", "unresolved"):
                lost[f"{r['species']} (pathogen={r['pathogen']})"] += 1
    nf = sum(v for k, v in n.items() if k != "host")
    print(f"\n{'=' * 72}\nCOVERAGE OF THE DETECTIONS by a GENUS-level USDA master\n{'=' * 72}")
    print(f"  fungal detections                 {nf:,}  (plus {n['host']:,} host rows)")
    for k in ("usda_species", "usda_genus_only", "outside", "unresolved"):
        print(f"    {k:<22}{n[k]:>8,}  ({100 * n[k] / nf:.1f}%)")
    print(f"\n  outside entirely, by taxon:")
    for k, v in sorted(lost.items(), key=lambda kv: -kv[1]):
        print(f"    {k:<52}{v:>5}")
    print(f"\n  a SPECIES-rank master would additionally lose {n['usda_genus_only']:,} "
          f"detections, of which {sum(gonly_path.values()):,} are pathogen-labelled")


def report_cost(tax, con, U_SP, U_GN):
    path = {t for (t,) in con.execute(
        "SELECT DISTINCT taxon_id FROM taxon_guild WHERE guild LIKE '%Plant Pathogen%'")}
    path |= {t for (t,) in con.execute(
        "SELECT DISTINCT a.fungus_id FROM host_association a "
        "JOIN taxon h ON h.taxon_id = a.host_id "
        "WHERE a.association_type='pathogenic' AND a.fungus_id IS NOT NULL "
        "AND h.kingdom='Viridiplantae'")}

    built = {}
    with open(SILVER / "cds_manifest.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            built[int(r["taxid"])] = r["organism_name"]
    keep = defaultdict(list)
    for ncbi, name in built.items():
        s, t = usda_status(tax, U_SP, U_GN, ncbi)
        keep[s].append((ncbi, name, t))
    print(f"\n{'=' * 72}\nWHAT THE RULE DROPS\n{'=' * 72}")
    print(f"  db_v4 built taxa {len(built):,}")
    for k in ("usda_species", "usda_genus_only", "outside", "unresolved"):
        verb = "kept" if k.startswith("usda") else "DROPPED"
        print(f"    {k:<20}{len(keep[k]):>7,}  {verb}")

    # Where the labels on the dropped taxa come from. A genus-inherited label is not
    # evidence about the species: see `seeds/taxon_status.tsv`, the Tier-C curation seed.
    sp_assigned, gn_inherited = [], []
    for ncbi, name, t in keep["outside"] + keep["unresolved"]:
        if t is None:
            continue
        gn = tax.ancestor_at(t, "genus")
        if t in path or tax.species_of(t) in path:
            sp_assigned.append(name)
        elif gn in path:
            gn_inherited.append(name)
    print(f"\n  of the dropped, carrying a plant-pathogen label: "
          f"{len(sp_assigned) + len(gn_inherited)}")
    print(f"    assigned at SPECIES rank {len(sp_assigned):>4}  "
          f"{', '.join(sorted(sp_assigned)) or '-'}")
    print(f"    inherited from the GENUS {len(gn_inherited):>4}")

    reads = defaultdict(int)
    with open(SILVER / "_classified.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kingdom"] != "Fungi":
                continue
            try:
                n = int(float(r["reads"] or 0))
            except ValueError:
                n = 0
            reads[usda_status(tax, U_SP, U_GN, int(r["taxid"]))[0]] += n
    tot = sum(reads.values())
    drop = reads["outside"] + reads["unresolved"]
    print(f"\n  reads the dropped taxa absorbed: {drop:,} of {tot:,} fungal "
          f"({100 * drop / tot:.2f}%). They are REASSIGNED to the nearest surviving taxon,")
    print(f"  which this measurement cannot predict. Only a rebuild settles it.")


def report_declared(tax, U_SP, U_GN):
    """The test that would break the chapter. Curated names, never a guild."""
    declared = set()
    with open(LIT_RUNS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for nm in (r.get("declared_pathogens") or "").replace(";", ",").split(","):
                nm = nm.strip()
                if nm and nm.lower() not in ("", "-", "na", "none"):
                    declared.add(nm)
    print(f"\n{'=' * 72}\nDECLARED PATHOGENS: {len(declared)} distinct names\n{'=' * 72}")
    for nm in sorted(declared):
        t = tax.resolve(nm)
        if not t:
            verdict = "UNRESOLVED"
        elif tax.species_of(t) in U_SP:
            verdict = "USDA species"
        elif tax.ancestor_at(t, "genus") in U_GN:
            verdict = "USDA genus only"
        else:
            verdict = "NOT REACHABLE"
        print(f"  {nm:<40}{verdict}")


def main():
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args()
    tax = _taxonomy.load()
    con = sqlite3.connect(HOLOBASE)
    L = load_sources(tax, con)
    report_sets(tax, L)
    for title, proj in projections(tax, con, L):
        report_projection(title, {k: proj(L[k]) for k in SOURCES})
    U_SP = {x for x in (tax.species_of(t) for t in L["usda"]) if x}
    U_GN = {x for x in (tax.ancestor_at(t, "genus") for t in L["usda"]) if x}
    report_coverage(tax, U_SP, U_GN)
    report_cost(tax, con, U_SP, U_GN)
    report_declared(tax, U_SP, U_GN)


if __name__ == "__main__":
    main()
