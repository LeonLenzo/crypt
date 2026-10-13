#!/usr/bin/env python3
"""scope_residents.py — what the reference database is missing, measured three ways.

db_v3 was seeded from PHI-base and expanded to congeners, so it is a plant PATHOGEN database.
But its stated design rule is that breadth exists to give the classifier COMPETITORS, not to
cover more pathogens: Kraken2's LCA can only decline to resolve a read when two or more
database taxa share its k-mers. A read from an organism with no database representation at all
does not go unclassified. It goes to the nearest taxon that is present. That failure mode is
invisible to every diagnostic in 05_filter, because the reads are real, the organism is real,
and only the assignment is wrong: the accumulation curve rises, coverage grows, the score
passes.

This script measures how much of that exposure exists. It answers three questions and writes
one table each.

1. Would core_nt fix it?  (--compare)
   core_nt is the obvious alternative: 12.4x db_v3's total minimizer space, and the database
   Colvin & Chopra used on two of our own BioProjects. The answer is no, and the reason is not
   size. See `data/db_comparison.tsv` and `data/corent_vs_dbv3_species.tsv`.

2. What is actually on a cereal leaf?  (--residents)
   Taken from Colvin & Chopra's per-sample fungal abundance tables for field maize and
   sorghum, which is data rather than our guess at the community. Each taxon is checked
   against db_v3. See `data/residents_colvin.tsv`.

3. Which of the gaps can be filled?  (--annotation, needs network)
   db_v3 is CDS-only: kraken_build.py skips any assembly without annotation, because these are
   RNA-seq libraries. So an assembly is not enough, a gene model is required. This queries NCBI
   for assembly and annotation counts per genus. See `data/resident_annotation.tsv`.

   The genera with zero annotated assemblies are not dead ends. They are inputs to the miniprot
   homology projection already built for the 23 unannotated blind-spot congeners
   (utilities/, and `crypt-gene-prediction` in memory).

Caveat that must travel with question 2: Colvin's table is itself a Kraken2 + Bracken + core_nt
output, not ground truth, and the artefact this script is about is visible in it. In SORGHUM it
puts Parastagonospora nodorum at 2.86%, Alternaria dauci at 1.43%, Plenodomus lingam at 1.24%
and Ascochyta rabiei at 1.22%, all at 100% prevalence. A Nebraska sorghum leaf does not carry
wheat, carrot, canola and chickpea pathogens at those levels. Those are Pleosporales reads from
unrepresented residents scattering across whichever Pleosporales their database held. Read the
output as "these lineages are present", at genus and family level, never as species truth.

Run from the repository root:

    python 03_kraken/01_search/scope_residents.py --compare --residents
    python 03_kraken/01_search/scope_residents.py --annotation      # network, ~2 min
    python 03_kraken/01_search/scope_residents.py --all
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


def _root() -> Path:
    here = Path(__file__).resolve()
    for d in (here, *here.parents):
        if d.is_dir() and (d / "CLAUDE.md").exists():
            return d
    raise SystemExit("repo root not found")


sys.path.insert(0, str(_root()))
from _paths import CLASSIFIED, COLVIN_DATA, CORENT_INSPECT, DB_INSPECT  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

# 05_filter's adopted operating point. Both thresholds, not one: score is the artefact filter,
# coverage is the abundance floor. See 03_kraken/05_filter/README.md.
SCORE_FLOOR, COV_FLOOR = 0.70, 0.01

# The cohort's hosts. Sorghum carries no detections yet because the kraken pass ran on the
# 02_literature cohort; the 1,120 literature-first sorghum runs have never been classified.
CEREALS = {"Triticum aestivum", "Zea mays", "Oryza sativa", "Sorghum bicolor", "Hordeum vulgare"}

FUNGI, OOMYCOTA, VIRIDIPLANTAE = 4751, 4762, 33090
GENOME_SCALE = 1_000_000  # minimizers; below this a taxon is a barcode, not a genome

# Genera named as resident or core by the cereal phyllosphere literature. Colvin covers
# sorghum, maize and soybean only, so wheat and rice come from deep searches recorded in
# README.md: wheat core is Dioszegia, Cryptococcus, Cladosporium, Sporobolomyces and
# Cystofilobasidiaceae; rice core is Nigrospora and Papiliotrema.
RESIDENT_GENERA = [
    ("basidio_yeast", "Sporobolomyces"), ("basidio_yeast", "Sporidiobolus"),
    ("basidio_yeast", "Rhodotorula"), ("basidio_yeast", "Rhodosporidiobolus"),
    ("basidio_yeast", "Dioszegia"), ("basidio_yeast", "Filobasidium"),
    ("basidio_yeast", "Vishniacozyma"), ("basidio_yeast", "Papiliotrema"),
    ("basidio_yeast", "Holtermanniella"), ("basidio_yeast", "Cystofilobasidium"),
    ("basidio_yeast", "Bullera"), ("basidio_yeast", "Udeniomyces"),
    ("basidio_yeast", "Tilletiopsis"), ("basidio_yeast", "Cryptococcus"),
    ("asco_epiphyte", "Aureobasidium"), ("asco_epiphyte", "Cladosporium"),
    ("asco_epiphyte", "Epicoccum"), ("asco_epiphyte", "Boeremia"),
    ("asco_epiphyte", "Phoma"), ("asco_epiphyte", "Leptosphaerulina"),
    ("asco_epiphyte", "Phaeosphaeria"), ("asco_epiphyte", "Nigrospora"),
    ("asco_epiphyte", "Curvularia"), ("asco_epiphyte", "Chaetomium"),
    ("asco_epiphyte", "Periconia"), ("asco_epiphyte", "Acremonium"),
    ("asco_epiphyte", "Microdochium"), ("asco_epiphyte", "Ramularia"),
    ("asco_epiphyte", "Sarocladium"),
]

# Families Colvin names as consistently abundant in crop phyllosphere, across ITS and
# metagenomic surveys as well as their own data. The check is whether db_v3 holds each.
COLVIN_FAMILIES = ["Aspergillaceae", "Cladosporiaceae", "Nectriaceae", "Pleosporaceae",
                   "Didymellaceae", "Glomeraceae", "Debaryomycetaceae", "Ustilaginaceae",
                   "Hypoxylaceae"]

NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/"


def _passes(row):
    """Does this detection clear both 05_filter floors?

    Test POSITIVELY (`>=`), never by negating (`not <`). 3,382 of the 24,467 rows in
    `_classified.tsv` (13.8%) carry `score = nan`, all of them `cls = unclassified`: the
    classifier could not fit a curve for that host x species cloud at all. NaN compares False
    against everything, so `if score < FLOOR: continue` does NOT skip them, it ADMITS them.
    That inflated the cereal species count from 32 to 44 here before it was caught, pulling in
    Bipolaris oryzae, Claviceps purpurea and Aspergillus flavus on the strength of coverage
    alone with no score at all. 393 of the NaN rows clear the 1% coverage floor, so this is not
    a rounding-scale error.
    """
    try:
        score, coverage = float(row["score"]), float(row["coverage"])
    except (TypeError, ValueError):
        return False
    return score >= SCORE_FLOOR and coverage >= COV_FLOOR


# ── inspect tables ────────────────────────────────────────────────────────────
def walk_inspect(path):
    """Yield (taxid, clade_minimizers, rank, name, ancestor_taxids) per taxon line.

    Lineage comes from the two-spaces-per-level indentation of the name column, so no external
    taxonomy dump is needed. Comment lines (core_nt ships a `# Database options:` header, our
    own dumps do not) have no tabs and fall out on the column count.
    """
    stack = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            raw = p[5]
            depth = (len(raw) - len(raw.lstrip(" "))) // 2
            try:
                tid, clade = int(p[4].strip()), int(p[1])
            except ValueError:
                continue
            while stack and stack[-1][0] >= depth:
                stack.pop()
            anc = {t for _, t in stack}
            stack.append((depth, tid))
            yield tid, clade, p[3].strip(), raw.strip(), anc


def load_inspect(path):
    """(clade_minimizers_by_taxid, rank_by_taxid, name_by_taxid, taxid_by_species_name)."""
    clade, rank, name, by_name = {}, {}, {}, {}
    for tid, c, r, n, _ in walk_inspect(path):
        clade[tid], rank[tid], name[tid] = c, r, n
        if r in ("S", "S1"):
            by_name.setdefault(n, tid)
    return clade, rank, name, by_name


def kingdom_stats(path):
    """Per-kingdom clade minimizers and species-node counts, plus the fungal ceiling spread."""
    out = {k: [0, 0] for k in (FUNGI, OOMYCOTA, VIRIDIPLANTAE)}
    total, fungal_ceilings = 0, []
    for tid, clade, rank, _, anc in walk_inspect(path):
        if tid == 1:
            total = clade
        for k in out:
            if tid == k:
                out[k][0] = clade
            if k in anc and rank == "S":
                out[k][1] += 1
        if FUNGI in anc and rank == "S":
            fungal_ceilings.append(clade)
    return total, out, sorted(fungal_ceilings)


# ── question 1: core_nt versus db_v3 ──────────────────────────────────────────
def compare():
    v3c, v3r, _, v3_by_name = load_inspect(DB_INSPECT)
    ntc, _, _, _ = load_inspect(CORENT_INSPECT)
    v3_total, v3_king, v3_fung = kingdom_stats(DB_INSPECT)
    nt_total, nt_king, nt_fung = kingdom_stats(CORENT_INSPECT)

    rows = []

    def add(metric, v3, nt):
        rows.append({"metric": metric, "db_v3": v3, "core_nt": nt})

    add("total_minimizers", v3_total, nt_total)
    add("taxon_nodes", len(v3c), len(ntc))
    add("species_rank_nodes", sum(1 for t in v3r if v3r[t] == "S"), None)
    for k, label in ((FUNGI, "fungi"), (OOMYCOTA, "oomycota"), (VIRIDIPLANTAE, "viridiplantae")):
        add(f"{label}_minimizers", v3_king[k][0], nt_king[k][0])
        add(f"{label}_species_nodes", v3_king[k][1], nt_king[k][1])
    add("fungal_species_median_ceiling", int(st.median(v3_fung)), int(st.median(nt_fung)))
    add("fungal_species_genome_scale",
        sum(1 for c in v3_fung if c > GENOME_SCALE),
        sum(1 for c in nt_fung if c > GENOME_SCALE))
    _write(DATA / "db_comparison.tsv", ["metric", "db_v3", "core_nt"], rows)

    # Per-species ceilings for the pathogens 05_filter actually detects on cereal hosts.
    # Matching is by taxid resolved from db_v3's own tree, so both sides key on the same node.
    detections = {}
    with open(CLASSIFIED, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["host"] not in CEREALS:
                continue
            if not _passes(row):
                continue
            detections[row["species"]] = detections.get(row["species"], 0) + 1

    sp_rows = []
    for species, n in sorted(detections.items(), key=lambda kv: -kv[1]):
        tid = v3_by_name.get(species)
        a = v3c.get(tid, 0) if tid else 0
        b = ntc.get(tid, 0) if tid else 0
        sp_rows.append({
            "species": species, "taxid": tid or "", "cereal_detections": n,
            "db_v3_minimizers": a, "core_nt_minimizers": b,
            "ratio_nt_over_v3": f"{b / a:.4f}" if a else "",
            "core_nt_genome_scale": "yes" if b > GENOME_SCALE else "no",
        })
    _write(DATA / "corent_vs_dbv3_species.tsv",
           ["species", "taxid", "cereal_detections", "db_v3_minimizers",
            "core_nt_minimizers", "ratio_nt_over_v3", "core_nt_genome_scale"], sp_rows)

    # What db_v3 would GAIN: genome-scale fungal species core_nt holds and we do not.
    v3_gs = {t for t, c in v3c.items() if c > GENOME_SCALE and v3r.get(t) == "S"}
    gain = []
    for tid, clade, rank, nm, anc in walk_inspect(CORENT_INSPECT):
        if rank == "S" and FUNGI in anc and clade > GENOME_SCALE and tid not in v3_gs:
            gain.append({"taxid": tid, "species": nm, "core_nt_minimizers": clade})
    gain.sort(key=lambda r: -r["core_nt_minimizers"])
    _write(DATA / "corent_only_genome_scale.tsv",
           ["taxid", "species", "core_nt_minimizers"], gain)

    print(f"  db_comparison.tsv              {len(rows)} metrics")
    print(f"  corent_vs_dbv3_species.tsv     {len(sp_rows)} detected species")
    print(f"  corent_only_genome_scale.tsv   {len(gain)} fungal species db_v3 lacks")
    print(f"  core_nt is {nt_total / v3_total:.1f}x db_v3 overall but holds "
          f"{nt_king[FUNGI][0] / v3_king[FUNGI][0]:.2f}x its fungal minimizers")


# ── question 2: what is actually on a cereal leaf ─────────────────────────────
def residents():
    v3_names, v3_taxids = set(), set()
    for tid, _, _, nm, _ in walk_inspect(DB_INSPECT):
        v3_names.add(nm)
        v3_taxids.add(str(tid))

    rows = []
    for fname, crop in (("dataset_S1_maize.csv", "maize"), ("dataset_S2_sorghum.csv", "sorghum")):
        path = COLVIN_DATA / fname
        if not path.exists():
            raise SystemExit(f"missing {path}; see README.md for the supplement URLs")
        with open(path, newline="") as fh:
            r = csv.reader(fh)
            next(r)
            for row in r:
                level, taxon, taxid = row[1], row[2], row[3]
                vals = [float(x) for x in row[4:] if x not in ("", "NA")]
                if not vals:
                    continue
                rows.append({
                    "crop": crop, "level": level, "taxon": taxon, "taxid": taxid,
                    "mean_rel_abundance": f"{st.mean(vals):.6f}",
                    "prevalence": f"{sum(1 for v in vals if v > 0) / len(vals):.4f}",
                    "in_db_v3": "yes" if (taxon in v3_names or taxid in v3_taxids) else "no",
                })
    rows.sort(key=lambda r: (r["crop"], r["level"], -float(r["mean_rel_abundance"])))
    _write(DATA / "residents_colvin.tsv",
           ["crop", "level", "taxon", "taxid", "mean_rel_abundance", "prevalence", "in_db_v3"],
           rows)

    # The family check needs ancestor NAMES rather than taxids, so walk the tree again keeping
    # the indentation stack as names. Cheap: db_v3's inspect table is 1,889 lines.
    fam_species, fam_present = {f: 0 for f in COLVIN_FAMILIES}, {f: False for f in COLVIN_FAMILIES}
    stack = []
    with open(DB_INSPECT) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            raw = p[5]
            depth = (len(raw) - len(raw.lstrip(" "))) // 2
            while stack and stack[-1][0] >= depth:
                stack.pop()
            anc_names = {n for _, n in stack}
            stack.append((depth, raw.strip()))
            if raw.strip() in fam_present:
                fam_present[raw.strip()] = True
            if p[3].strip() == "S":
                for f in COLVIN_FAMILIES:
                    if f in anc_names:
                        fam_species[f] += 1
    fam = [{"family": f, "in_db_v3": "yes" if fam_present[f] else "no",
            "db_v3_species": fam_species[f]} for f in COLVIN_FAMILIES]
    _write(DATA / "colvin_families.tsv", ["family", "in_db_v3", "db_v3_species"], fam)

    for crop in ("maize", "sorghum"):
        sp = [r for r in rows if r["crop"] == crop and r["level"] == "Species"]
        tot = sum(float(r["mean_rel_abundance"]) for r in sp)
        miss = sum(float(r["mean_rel_abundance"]) for r in sp if r["in_db_v3"] == "no")
        print(f"  {crop:8s} {len(sp):3d} species, "
              f"{sum(1 for r in sp if r['in_db_v3'] == 'no'):3d} absent from db_v3, "
              f"{miss / tot * 100:4.1f}% of abundance on absent taxa")
    print(f"  residents_colvin.tsv           {len(rows)} rows")
    print(f"  colvin_families.tsv            {sum(1 for f in fam if f['in_db_v3'] == 'no')}"
          f" of {len(fam)} dominant families absent")


# ── question 3: can the gap be filled ────────────────────────────────────────
def _ncbi_total(taxon, annotated):
    q = {"page_size": 1, "filters.assembly_version": "current"}
    if annotated:
        q["filters.has_annotation"] = "true"
    url = NCBI + urllib.parse.quote(taxon) + "/dataset_report?" + urllib.parse.urlencode(q)
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.load(fh).get("total_count", 0)
        except Exception:
            time.sleep(4)
    return None


def annotation():
    """db_v3 is CDS-only, so an assembly without a gene model cannot be added as-is."""
    rows = []
    for group, genus in RESIDENT_GENERA:
        total, ann = _ncbi_total(genus, False), _ncbi_total(genus, True)
        rows.append({
            "group": group, "genus": genus, "assemblies": total, "annotated": ann,
            "pct_annotated": f"{ann / total * 100:.0f}" if total and ann is not None else "",
            "route": "add" if (ann or 0) > 0 else "miniprot" if (total or 0) > 0 else "none",
        })
        time.sleep(0.4)
    _write(DATA / "resident_annotation.tsv",
           ["group", "genus", "assemblies", "annotated", "pct_annotated", "route"], rows)
    for route in ("add", "miniprot", "none"):
        n = sum(1 for r in rows if r["route"] == route)
        print(f"  route={route:9s} {n:2d} genera")
    print(f"  resident_annotation.tsv        {len(rows)} genera")


def _write(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compare", action="store_true", help="core_nt versus db_v3 content")
    ap.add_argument("--residents", action="store_true", help="Colvin resident abundances")
    ap.add_argument("--annotation", action="store_true", help="NCBI annotation counts (network)")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        a.compare = a.residents = a.annotation = True
    if not (a.compare or a.residents or a.annotation):
        ap.error("pick at least one of --compare / --residents / --annotation / --all")
    if a.compare:
        print("core_nt versus db_v3:")
        compare()
    if a.residents:
        print("resident community (Colvin & Chopra 2026):")
        residents()
    if a.annotation:
        print("annotation availability:")
        annotation()


if __name__ == "__main__":
    main()
