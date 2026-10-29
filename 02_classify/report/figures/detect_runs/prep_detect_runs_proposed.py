#!/usr/bin/env python3
"""prep_detect_runs_proposed.py -- PROPOSAL ONLY, 2026-10-28. Re-tiers detect_runs.tsv.

Not a replacement for prep_detect_runs.py. It reads that script's output and rewrites two
columns so the proposed colour scheme can be rendered on the SAME RUNS and compared side by
side. If the scheme is adopted, this logic moves into prep_detect_runs.py and this file goes.

WHAT CHANGES. The current `tier` answers "which database listed this taxon": PHI-base, then
FUNGuild, then the membership routes. The proposed one answers "what is this organism to THIS
plant", and evidence grade moves off colour into its own table column.

    colour   pathogen of this crop / pathogen of another host / recorded on this crop,
             not a pathogen / recorded on other plants only / no plant record / host
    letter   A PHI-base, B FUNGuild species-rank, C curated (taxon_status.tsv)

The host axis comes from USDA `recorded_on` in holobase, matched at GENUS. It resolves 93% of
called pathogen rows, and it contradicts the current labels: 147 of the 265 rows coloured
"PHI-base: pathogen of another plant" are recorded by USDA on the crop that was sampled. That
tier reports which host a PHI-base experiment happened to use, not where the fungus lives.

`taxon_status.tsv` is applied as the agreed override: status=='pathogen' makes a taxon a
pathogen, every other status makes it not one, regardless of what PHI-base or FUNGuild say.
That is the same wiring 05_detect is due to get; doing it here keeps the figure honest about
what the curation decided without pre-empting the detect.py change.

    python 02_classify/report/figures/detect_runs/prep_detect_runs_proposed.py
"""

import csv
import glob
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import CLASSIFIED, HOLOBASE, REFDATA, SEEDS  # noqa: E402
from _util import atomic_write  # noqa: E402

import argparse
_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--all", action="store_true", help="re-tier the full detect_runs_all.tsv")
_ARGS = _ap.parse_args()
SRC = HERE / ("detect_runs_all.tsv" if _ARGS.all else "detect_runs.tsv")
TSV = HERE / ("detect_runs_all_proposed.tsv" if _ARGS.all else "detect_runs_proposed.tsv")
STATUS = SEEDS / "taxon_status.tsv"


def usda_hosts() -> dict:
    """fungal ncbi_taxid -> {plant host GENERA it has been recorded on}.

    `recorded_on`, not `pathogenic`: USDA is the breadth source and carries no pathogenicity,
    which is exactly right here. The question this answers is "has this fungus ever been found
    on this crop", and the pathogen/not decision is made separately and before it.
    """
    con = sqlite3.connect(f"file:{HOLOBASE}?mode=ro", uri=True)
    out = defaultdict(set)
    for tid, host in con.execute(
            "SELECT f.ncbi_taxid, h.scientific_name FROM host_association a "
            "JOIN taxon f ON f.taxon_id = a.fungus_id "
            "JOIN taxon h ON h.taxon_id = a.host_id "
            "WHERE a.source = 'usda_fungus_host' AND h.kingdom = 'Viridiplantae' "
            "AND f.ncbi_taxid IS NOT NULL"):
        if host:
            out[int(tid)].add(host.split()[0])
    con.close()
    return out


# ORDERS WHOSE MEMBERS HAVE NO SAPROTROPHIC PHASE. Rusts, powdery mildews and smuts cannot
# grow on dead tissue and cannot be cultured, so transcripts from one at depth in an RNA-seq
# library require living, infected host tissue. This is the only class of detection to which
# "DNA presence is not proof of active infection" does not apply, which is why it is worth
# separating. Textbook mycology, and unlike pathogenicity an ORDER is safely inheritable to
# every species in it.
#   Smuts are the soft edge: Ustilaginales are biotrophic in planta but have a culturable
# haploid yeast phase. They are 8 of 1,617 calls, so the grouping does not turn on them.
BIOTROPH_ORDERS = {"Pucciniales", "Erysiphales", "Ustilaginales", "Urocystidales",
                   "Taphrinales", "Entylomatales"}


def fungaltraits() -> dict:
    """fungal GENUS -> (order, primary_lifestyle), from Polme et al. 2020.

    Genus-keyed, and that is fine HERE in a way it is not for pathogen status. The trap the
    species-rank rule exists to stop is inheriting "plant pathogen" from a genus that contains
    both pathogens and biocontrol yeasts. Trophic habit does not work like that: every member
    of Pucciniales is an obligate biotroph, because that is what the order IS. The lifestyle
    fallback for the non-biotroph split is weaker and only decides opportunist-versus-not.
    """
    hits = sorted(glob.glob(str(REFDATA / "databases/fungaltraits/*/"
                                "polme2020-s1-fungal-traits-genera.tsv")))
    if not hits:
        print("  no FungalTraits table; trophic class will be blank", file=sys.stderr)
        return {}
    out = {}
    with open(hits[-1]) as fh:                       # last by date, the newest pull
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["GENUS"]] = (r["Order"], r["primary_lifestyle"] or "")
    return out


def trophic_of(species: str, ft: dict) -> str:
    """biotroph | necrotroph | opportunist, for a taxon already judged a pathogen."""
    hit = ft.get(species.split()[0])
    if not hit:
        return ""
    order, life = hit
    if order in BIOTROPH_ORDERS:
        return "biotroph"
    # NOT a plant_pathogen genus, yet something called it a pathogen at species rank. These
    # are the ubiquitous colonists a reviewer singles out first (Cladosporium, Aspergillus,
    # Penicillium), and saying so is better than letting them sit beside Exserohilum.
    return "necrotroph" if life == "plant_pathogen" else "opportunist"


def curated() -> dict:
    """species -> status, from the curated Tier C override."""
    out = {}
    with open(STATUS) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            if r.get("species") and r.get("status"):
                out[r["species"]] = r["status"]
    return out


def main() -> None:
    hosts = usda_hosts()
    cur = curated()
    ft = fungaltraits()
    print(f"  FungalTraits: {len(ft):,} genera")
    print(f"  USDA recorded_on: {len(hosts):,} fungal taxids with a plant host")
    print(f"  curated override: {len(cur)} species, "
          f"{sum(1 for v in cur.values() if v == 'pathogen')} pathogen")

    rows = list(csv.DictReader(open(SRC), delimiter="\t"))
    runs = {r["run"] for r in rows}
    # detect_runs.tsv does not carry taxid, and the USDA join needs one. Pull it from the
    # classified table for just these runs rather than widening the upstream script.
    taxid, routes = {}, {}
    for r in csv.DictReader(open(CLASSIFIED), delimiter="\t"):
        if r["run"] in runs:
            taxid[(r["run"], r["species"])] = r["taxid"]
            routes[(r["run"], r["species"])] = {x for x in (r["routes"] or "").split(";") if x}

    moved, counts = 0, defaultdict(int)
    for r in rows:
        old = r["tier"]
        sp, src = r["species"], (r["pathogen_source"] or "").strip()
        tid = taxid.get((r["run"], sp))
        host_genus = (r["host"] or "").split()[0] if r.get("host") else ""
        recorded = hosts.get(int(tid), set()) if tid else set()
        # USDA FIRST, membership route as the FALLBACK where USDA is silent. USDA is the
        # better evidence (an explicit host record) but it is not the only way a taxon is
        # known to sit on these crops: db_v4 seeded the phyllosphere residents from measured
        # abundance in field maize and sorghum, and those carry a route and no USDA row.
        # Without the fallback, Dioszegia hungarica and Sporobolomyces ruberrimus -- the
        # wheat and maize resident core -- land in "no plant host record", which is false.
        rts = routes.get((r["run"], sp), set())
        on_this = bool(host_genus and host_genus in recorded) or "host_associate" in rts
        has_record = bool(recorded) or bool(rts & {"host_associate", "plant_associate"})

        status = cur.get(sp)
        # The ladder, in order. PHI-base and FUNGuild still decide pathogen status for
        # everything the curation did not reach; the curation overrides both where it did.
        if status == "pathogen":
            is_path, grade = True, "C"
        elif status:                       # curated as mycoparasite/biocontrol/resident/...
            is_path, grade = False, "C"
        elif old in ("phibase_crop", "phibase_plant", "funguild"):
            is_path = True
            grade = "A" if src in ("phibase", "both") else "B"
        else:
            is_path, grade = False, ""

        if (r["kingdom"] or "").strip() == "Viridiplantae":
            tier, grade = "host", ""
        elif is_path and on_this:
            tier = "pathogen_crop"
        elif is_path:
            tier = "pathogen_other"
        elif on_this:
            tier = "resident_crop"
        elif has_record:
            tier = "resident_other"
        else:
            tier = "no_record"

        # TROPHIC CLASS, on pathogens only. It answers "what kind of organism is this",
        # which is a different question from the tier's "what is it to this plant", so the
        # figure gives it a different channel (shape) rather than another hue.
        r["tier"], r["evidence"] = tier, grade
        r["trophic"] = trophic_of(sp, ft) if tier in ("pathogen_crop", "pathogen_other") else ""
        if r["call"] == "yes":
            counts[(old, tier)] += 1
            moved += old != tier

    with atomic_write(TSV) as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    called = sum(1 for r in rows if r["call"] == "yes")
    print(f"-> {TSV.name}: {len(rows):,} candidates, {called:,} called, "
          f"{moved:,} called rows change tier")
    for (old, new), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        flag = "" if old == new else "  <-- moved"
        print(f"   {old:20s} -> {new:16s} {n:5,}{flag}")
    tro = defaultdict(int)
    for r in rows:
        if r["call"] == "yes" and r["trophic"]:
            tro[r["trophic"]] += 1
    print("   trophic class of called pathogens: "
          + ", ".join(f"{k} {v}" for k, v in sorted(tro.items(), key=lambda kv: -kv[1])))


if __name__ == "__main__":
    main()
