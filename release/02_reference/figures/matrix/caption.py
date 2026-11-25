#!/usr/bin/env python3
"""caption.py — the caption for panel A, written from the curation tables.

    python release/02_reference/figures/matrix/caption.py

Reads  matrix_numbers.tsv                              (beside this script)
       ../../data/outputs/supp_data_s3_reference.tsv    the published reference table
       ../../data/curation/screen_pathogens.tsv, screen_merges.tsv   for provenance only
       ../../data/curation/screen_species.tsv           to split the genus-only row
Writes matrix_caption.txt

The counts come from S3, the same table the panel is drawn from, so the caption cannot
narrate a different figure. The two curation reads are for things S3 cannot carry: the
`aff.` nodes that were deliberately NOT merged, and a labelled taxon with no in-scope
assembly, neither of which is a row in the reference.

Nothing here is prose-with-numbers-typed-in: the narrative is a template and every figure in
it is read at write time. Rerun after any pipeline change.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]          # release/02_reference
FIG = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "scripts"))
from _caption import rd, para, write                # noqa: E402

CUR = HERE / "data/curation"
S3 = HERE / "data/outputs/supp_data_s3_reference.tsv"
SP_NODE = re.compile(r"(?<![A-Za-z])sp\.")          # an unidentified NCBI node


def provisional_counts(ref):
    """How much of each tier is an unidentified `sp.` node, read off the taxon names."""
    gen = [r for r in ref if r["route"] == "usda_genus"]
    obs = [r for r in ref if r["route"] != "usda_genus"]
    return (sum(1 for r in gen if SP_NODE.search(r["taxon"])), len(gen),
            sum(1 for r in obs if SP_NODE.search(r["taxon"])), len(obs))


def genus_tier_split(ref, species):
    """-> (sibling, unbuildable): what the genus ROUTE actually holds.

    One label over two different situations, and the distinction is what the old word
    "congener" hid. Either USDA names a SIBLING species of this taxon and that sibling is
    itself in the reference, which is a congener in the strict sense; or USDA names species
    in the genus and none of them could be built, so the genus qualified on an observation
    that never became a row of its own. A third case, a genus USDA never reaches at all, is
    EMPTY by construction of the scope rule and is asserted here, so a scope leak fails
    loudly rather than being absorbed into one of the other two.

    USDA names are read through the pipeline's own resolution, never by splitting the raw
    string: a raw-name match reports *Dioszegia* as absent from USDA when it is there as
    *Bullera*, and synonymy is load-bearing throughout this module.
    """
    in_ref = {r["genus"] for r in ref if r["route"] == "usda_species"}
    usda_sp = {r["matched_name"].strip("[]").split()[0] for r in species
               if r["matched_name"] and r["rank"] not in ("", "genus")}
    sib = unb = leak = 0
    for r in ref:
        if r["route"] != "usda_genus":
            continue
        if r["genus"] in in_ref:
            sib += 1
        elif r["genus"] in usda_sp:
            unb += 1
        else:
            leak += 1
    if leak:
        raise SystemExit(f"{leak} genus-route taxa sit in genera USDA never names. The scope "
                         f"rule cannot produce that; fix the scope, not the caption.")
    return sib, unb


def main():
    M = {(r["row"], r["col"]): (int(r["taxa"]), int(r["assemblies"]), int(r["cds"]),
                                int(r["predicted"]))
         for r in rd(FIG / "matrix_numbers.tsv")}
    # Fungi only: S3 also publishes the host sink (route=host, kingdom=Viridiplantae). The
    # three host rows would not change any number this caption quotes -- no pathogen taxid is
    # a plant and no host name parses as an `sp.` node -- but the caption describes the fungal
    # scope, so it reads the fungal scope.
    ref = [r for r in rd(S3) if r["kingdom"] == "Fungi"]
    pth = rd(CUR / "screen_pathogens.tsv")
    mrg = rd(CUR / "screen_merges.tsv")
    spc = rd(CUR / "screen_species.tsv")
    scope = {int(r["ncbi_taxid"]) for r in ref}

    coh_t, coh_a, coh_c, coh_p = M[("all", "cohort")]
    gen_sp, gen_n, obs_sp, obs_n = provisional_counts(ref)
    sib_n, unb_n = genus_tier_split(ref, spc)
    pred = [r for r in ref if int(r["n_predicted"]) > 0]
    pred_wheat = sum(1 for r in pred if r["usda_row"] == "wheat_resident")
    # Every predicted taxon had no NCBI annotation; that is what put it on the worklist.
    # Asserted rather than asserted-in-prose, since the caption makes the claim below.
    assert all(int(r["n_annotated"]) == 0 for r in pred), \
        "a predicted taxon already carried an NCBI annotation; the Pred. wording is wrong"
    n_merge = sum(1 for r in mrg if r["action"] == "merge")
    n_reroute = sum(1 for r in mrg if r["action"] == "reroute")
    n_keep = sum(1 for r in mrg if r["action"] == "keep")
    asm_moved = sum(int(r["assemblies_moved"]) for r in mrg if r["action"] != "keep")

    wheat = {int(r["ncbi_taxid"]): r for r in pth if r["is_wheat_pathogen"] == "yes"}
    n_cur = sum(1 for r in wheat.values() if r["wheat_pathogen_source"] == "curated")
    n_phi = len(wheat) - n_cur
    absent = [r["species"] for t, r in wheat.items() if t not in scope]
    host_only = [r["species"] for r in pth if r["plant_pathogen_evidence"] == "host_only"]

    rows = []
    for rn, lab in (("wheat_resident", "Wheat Resident"),
                    ("plant_resident", "Plant Resident")):
        for cn, cl in (("none", "Non pathogen"), ("vagrant", "Plant pathogen"),
                       ("wheat_path", "Wheat pathogen")):
            t, a, c, _ = M[(rn, cn)]
            if t:
                rows.append(f"  {lab:<16}{cl:<17}{t:>8,}{a:>13,}{c:>9,}{M[(rn, cn)][3]:>8,}")
    table = "\n".join(rows)
    absent_note = (
        f" The remaining {len(absent)} ({', '.join(absent)}) "
        f"{'has' if len(absent) == 1 else 'have'} no in-scope assembly, so "
        f"{'it carries' if len(absent) == 1 else 'they carry'} the label but "
        f"{'is' if len(absent) == 1 else 'are'} not in the reference and "
        f"{'does' if len(absent) == 1 else 'do'} not appear in the panel."
        if absent else "")

    write(FIG / "matrix_caption.txt", para(
        "The Reference Database, by Record and by Pathogenicity",
        "Caption for matrix.png (panel A).",
        "",
        "GENERATED by caption.py — do not hand-edit, the numbers will drift.",
        "",
        "",
        "FULL CAPTION",
        "",
        f"Composition of the reference database: {coh_t:,} fungal taxa carrying {coh_a:,} "
        f"genome assemblies, of which {coh_c:,} have a gene model and can contribute CDS. "
        f"Each taxon is placed by two independent facts and the matrix crosses them. Rows "
        f"are whether the USDA Plant-Fungal Database names that species on wheat or does "
        f"not. Columns are what PHI-base scores: no pathogenicity record, a named disease "
        f"on some plant, or a named disease on wheat. Cells are disjoint and sum to the "
        f"cohort: every taxon is in exactly one.",
        "",
        "Each cell reports four counts. Taxa is distinct NCBI taxa; Asm. is the genome "
        "assemblies available for them; CDS is the subset of those assemblies carrying a "
        "gene model DOWNLOADED FROM NCBI, which is what the database could be built from "
        "before this work. The gap between Asm. and CDS is the annotation bottleneck that "
        "motivates gene prediction.",
        "",
        f"Pred. is the fourth count and is the part of that gap this work closed: gene "
        f"models PREDICTED here, ab initio, for taxa NCBI has never annotated "
        f"(research/predict; validated against BUSCO and against an annotated congener's "
        f"k-mer space). All {len(pred)} are taxa whose NCBI annotation count is zero, so "
        f"Pred. never overlaps CDS and a predicted taxon was uncallable before db2. "
        f"{pred_wheat} of the {len(pred)} are in the wheat row, which is where the worklist "
        f"aimed them: a fungus USDA records on wheat, in the reference, with no gene model "
        f"is the most expensive kind of silence this database can hold.",
        "",
        "  Row             Column               Taxa   Assemblies      CDS   Pred.",
        table,
        "",
        "Colour carries both axes, so the panel needs no key. Hue is the column: green for "
        "no pathogenicity record, amber for a pathogen of other plants, orange-red for a "
        "wheat pathogen. Saturation is the row: full strength where USDA names the taxon "
        "on wheat, pale where it does not. Panel B uses the same five colours for the same "
        "five cells.",
        "",
        "",
        "WHAT THE EMPTY CELL MEANS",
        "",
        f"One of the six combinations does not occur and is NOT DRAWN. No wheat pathogen in "
        f"scope is off-wheat in USDA: all {M[('wheat_resident', 'wheat_path')][0]} sit in "
        f"the top row. USDA host records and PHI-base pathogenicity scores are independent "
        f"compilations made for different purposes, so their agreeing here is an empirical "
        f"check on the scope rule rather than a definition, and it is the test the design "
        f"has to pass. The gap is a result and not a layout choice.",
        "",
        f"A second structural absence is real but is no longer visible as a gap, because "
        f"the rows no longer separate it. Of the {gen_n:,} taxa that entered the reference "
        f"on their GENUS rather than under their own name, not one carries a pathogen "
        f"label of either kind. A fungus PHI-base scores with a named disease is a studied "
        f"fungus and USDA records studied fungi, so the label cannot reach a taxon nobody "
        f"has recorded. No detection can be called a pathogen on genus membership alone.",
        "",
        "",
        "HOW THE TWO AXES ARE DEFINED",
        "",
        f"The unit is a TAXON, not a described species, and the panel says so. {gen_sp:,} "
        f"of the {gen_n:,} genus-only taxa are NCBI `sp.` nodes: unidentified isolates "
        f"carrying "
        f"a strain code and no binomial. They are legitimate reference entries, because an "
        f"unidentified Alternaria gives the classifier as real a competitor as a named one, "
        f"but they are not species in any taxonomic sense. The observation rows do not have "
        f"the problem ({obs_sp} of {obs_n:,}), because USDA indexes only formal names.",
        "",
        f"THE PLANT ROW POOLS TWO ROUTES INTO THE REFERENCE AND DOES NOT SEPARATE THEM, "
        f"deliberately. {gen_n:,} of its taxa are there because USDA recorded their GENUS "
        f"and not the species itself. How a taxon reached the reference is a fact about the "
        f"database build, recorded per taxon in the `route` column of "
        f"supp_data_s3_reference.tsv, and drawing it as a row made the figure assert "
        f"something about the ORGANISM that the route cannot support: USDA Fungus-Host is a "
        f"herbarium records database, so an unrecorded name is as likely to be a gap in its "
        f"coverage as a fact about the fungus. Two cases arrive by that route and nothing "
        f"separates them: {sib_n:,} taxa whose genus is represented here by a USDA-recorded "
        f"sibling species, and {unb_n:,} whose genus qualified on USDA records that could "
        f"not be built into a row of their own, mostly lichens, mushrooms and mycorrhizas "
        f"with no annotated assembly. Treating the route as weak evidence would discard "
        f"real detections, and a draft of a downstream figure did exactly that: it dropped "
        f"*Dioszegia hungarica*, a common wheat phyllosphere yeast USDA holds only under "
        f"its former name *Bullera*.",
        "",
        f"A PHI-base row counts only where it NAMES A DISEASE. PHI-base scores an "
        f"interaction against whatever host the experiment used, and the plant is not "
        f"always the host being harmed: in a tripartite insect-fungus-plant assay, or an "
        f"endophyte colonisation study, the plant is present but is not the victim, and "
        f"those rows carry Disease_name = \"no data found\". {len(host_only)} fungi in this "
        f"reference were labelled plant pathogens on such rows alone and are not labelled "
        f"here (" + ", ".join(host_only) + f"). Each is a plant associate rather than a "
        f"pathogen. The condition applies to the wheat column too, where it changes nothing.",
        "",
        f"The wheat pathogen label covers {len(wheat)} taxa: {n_phi} from PHI-base and "
        f"{n_cur} from a curated supplement (data/inputs/wheat_pathogens_curated.tsv) "
        f"carrying the canonical wheat diseases PHI-base does not curate. "
        f"{M[('wheat_resident', 'wheat_path')][0]} of them are in the reference."
        + absent_note +
        " Provenance is recorded per taxon in wheat_pathogen_source; no call is anonymous.",
        "",
        f"Before counting, {n_merge} provisional NCBI names were folded onto the USDA "
        f"binomial they belong to and {n_reroute} was promoted in place, moving "
        f"{asm_moved} assemblies. USDA and NCBI can hold one fungus under two names, which "
        f"otherwise splits it into an observation with no genome and a genus-only row with "
        f"no host record. {n_keep} `aff.` nodes were deliberately NOT merged, because `aff.` "
        f"denotes a related undescribed species and those belong in the genus-only row. "
        f"Every decision is a row in data/curation/screen_merges.tsv.",
        "",
        "",
        "WHAT IS NOT SHOWN",
        "",
        "Scope is ANY plant, never wheat. The database deliberately holds fungi with no "
        "wheat record and no pathogenicity record, because a reference absent from the "
        "database is invisible: its reads go silently to the nearest taxon present. The "
        "columns are annotations applied to a detection, never filters on the database.",
        "",
        "The panel counts what the reference CONTAINS, not what was detected. No cell is a "
        "claim about abundance, prevalence or co-infection.",
        "",
        "",
        "SHORT CAPTION",
        "",
        f"Composition of the reference database: {coh_t:,} fungal taxa carrying {coh_a:,} "
        f"genome assemblies, {coh_c:,} of them with a gene model from NCBI and {len(pred)} "
        f"with one predicted here. Each taxon is placed by whether the USDA Plant-Fungal "
        f"Database names it on wheat (rows) crossed with what PHI-base scores (columns: no "
        f"pathogenicity record, a named disease on another plant, a named disease on "
        f"wheat). Cells are disjoint and sum to the cohort; hue is the column and "
        f"saturation the row, and panel B reuses the same five colours. One combination "
        f"does not occur and is not drawn: no wheat pathogen in scope is off-wheat in "
        f"USDA.",
    ))


if __name__ == "__main__":
    main()
