#!/usr/bin/env python3
"""_exclusions.py — BioProjects judged unusable, named once, with the reason.

Some BioProjects are in the cohort because the pipeline put them there and are wrong on
inspection. Until this file existed, each such finding lived only as a sentence in a README,
which meant it had to be re-derived every time someone asked "can we use this one?" and it
never actually changed a number. This is the registry: a decision made once, applied at the
analysis stage, and visible in the output.

**Exclude here, never upstream.** Per the clean-run rule, every sample is downloaded,
classified and scored identically; nothing is skipped because of what we expect to find.
Exclusion is an analysis-stage judgment about what a sample can support, so it is applied when
a rate is computed, not when reads are fetched. That keeps the excluded samples measurable:
you can always report the cohort with and without them, which is the only way to show the
exclusion is not doing the work of the result.

Every entry needs `reason` (what the study actually is), `evidence` (how we know, specifically
enough that someone else can check it), and `decided` (when). An entry without evidence is an
opinion and does not belong here.

Usage:

    import sys; sys.path.insert(0, str(_find_root()))
    from _exclusions import EXCLUDED, is_excluded, split_excluded

    kept, dropped = split_excluded(rows, key="BioProject")

Related: `_paths.py` for cross-module paths, `_util.py` for domain-free helpers.
"""
from pathlib import Path

_ROOT_MARKERS = ("CLAUDE.md", ".gitignore")


def _find_root(start=None):
    here = Path(start or __file__).resolve()
    for d in (here, *here.parents):
        if d.is_dir() and any((d / m).exists() for m in _ROOT_MARKERS):
            return d
    raise RuntimeError(f"repo root not found above {here}; expected one of {_ROOT_MARKERS}")


# ── The registry ─────────────────────────────────────────────────────────────────────────────
#
# `kind` is why the study cannot support a field co-infection claim, not how bad it is:
#
#   inoculated   deliberate pathogen introduction, so any secondary signal is the design
#   in-vitro     no plant host in the sequenced material, or an artificial surface
#   archival     culture-collection or herbarium material, not a contemporaneous collection
#
# `n_biosamples` is the count in the classified cohort at the time of the decision. It is
# recorded so a later change in cohort size is visible rather than silent; it is not used to
# match anything.

EXCLUDED = {
    "PRJNA1217477": dict(
        kind="inoculated",
        n_biosamples=211,
        reason=(
            "Multiplant transcriptomic atlas of Botrytis cinerea infection at UC Davis. Every "
            "plant was deliberately inoculated, so secondary signal is the experiment rather "
            "than an incidental co-infection, and the cohort's largest apparent field site is "
            "a laboratory."
        ),
        evidence=(
            "Classified llm_study_setting == 'field' and geo_loc_name 'USA: California, Davis', "
            "which is the institution, not a collection site. Its PNAS supplement "
            "(pnas.2601719123.sd07.xlsx) tabulates the origins of the inoculum ISOLATES, "
            "Californian vineyards, not the glasshouse plants that were sequenced; merging it "
            "would have stamped vineyard coordinates onto Arabidopsis. 8% of the analysed "
            "field cohort on its own."
        ),
        decided="2026-10-02",
    ),
    "PRJNA526829": dict(
        kind="in-vitro",
        n_biosamples=60,
        reason=(
            "Magnaporthe oryzae MoAa91 chitin-binding study, conducted in vitro on artificial "
            "surfaces. There is no field-collected plant tissue in it."
        ),
        evidence="Classified llm_study_setting == 'field'; the title and methods describe in vitro work.",
        decided="2026-10-02",
    ),
    "PRJNA1314945": dict(kind="wrong-assay", n_biosamples=306,
        reason="Not RNA-seq at all. SLAF-seq is Specific-Locus Amplified Fragment sequencing, a "
               "reduced-representation GENOMIC DNA method for SNP genotyping. There is no "
               "transcriptome here, so there are no fungal transcripts to find.",
        evidence="Paper is Yang et al. 2026, BMC Genomics, 'SLAF-seq efficiently identifies SNP "
                 "markers for wheat (Triticum aestivum L.) improvement' (10.1186/s12864-026-12921-6, "
                 "PMC13471324), which names PRJNA1314945 as its deposit. It describes 'a "
                 "reduced-representation genome sequencing technology' in which 'two optimal "
                 "restriction enzymes are chosen for double digestion' and 'Qualified genomic DNA "
                 "samples were individually digested'. The word RNA does not appear anywhere in the "
                 "paper. The ENA study title is literally 'SLAF-seq Date of 306 wheat varieties'. "
                 "DESPITE THIS the archive records library_strategy=RNA-Seq, library_selection=cDNA "
                 "and library_source=TRANSCRIPTOMIC on all 306 runs, so offtarget_screen.py passed it "
                 "as 'plant organism, all RNA-Seq'. The submitter's assay metadata is simply wrong.",
        decided="2026-10-06"),
    "PRJNA1428298": dict(kind="inoculated", n_biosamples=340,
        reason="Pilot experiment for PRJNA1217477, the same UC Davis Botrytis atlas, and excluded "
               "for the same reason: every plant was deliberately inoculated, so a secondary "
               "detection is the experiment rather than a cryptic co-infection.",
        evidence="NCBI project_description names its parent outright: 'This dataset includes pilot "
                 "experiment data generated to support large-scale comparative analyses across "
                 "multiple eudicot species infected with 72 Botrytis isolates, i.e., PRJNA1217477: "
                 "Co-transcriptomic networks to identify conserved and lineage-specific mechanisms "
                 "for core eudicot defenses against a generalist pathogen'. Same submitter "
                 "(University of California, Davis). Title: 'Pilot co transcriptomic experiment for "
                 "large scale analysis of eudicot defenses against 72 Botrytis cinerea isolates'. "
                 "Found 2026-10-06 by offtarget_screen.py, which flagged it as organism 'mixed "
                 "sample' with 340 RNA-Seq runs; it shares a paper (10.1073/pnas.2601719123) with "
                 "PRJNA1217477 but was never linked to the exclusion.",
        decided="2026-10-06"),
    "PRJNA328045": dict(
        kind="inoculated",
        n_biosamples=39,
        reason=(
            "Comparative transcriptomics of Leptosphaeria maculans virulence factors, using "
            "controlled compatible and incompatible inoculations."
        ),
        evidence="Classified llm_study_setting == 'field'; the design is deliberate inoculation.",
        decided="2026-10-02",
    ),
}


# ── Tissue: FOLIAR, which is narrower than the cohort was curated ───────────────────────────
#
# The cohort was screened for field-grown AERIAL tissue, and the Undermind prompt says so
# explicitly: "leaf, stem, spike, head, panicle, flower, fruit". An ear and a kernel are
# aerial, so they are in scope by that definition and are NOT a curation error.
#
# The ANALYSIS assumes foliar, and that assumption is load-bearing in two places: the
# FungalTraits pathogenic-capacity axis is leaf/fruit/seed, and the resident community the
# db_v4 scope was built against is the PHYLLOSPHERE (Dioszegia, Cryptococcus, Cladosporium,
# Sporobolomyces). A developing kernel is a different organ with a different mycobiome.
#
# Leon's call 2026-10-29: narrow to foliar. This NARROWS A DOCUMENTED SCOPE DECISION, so it is
# recorded here rather than applied silently, and it is applied at analysis stage like every
# other exclusion so the set stays measurable with and without.
#
# Measured cost: 19 of the 3,624 retained runs. 17 `immature kernel` from PRJNA383416 and 2
# `ear` from PRJNA291919, both maize, both projects otherwise valid and retained in full for
# their leaf samples. The other non-foliar tissue in the wider cohort (254 kernel in
# PRJNA383416, 112 in PRJNA395908, 45 seeds, 36 kernel in PRJNA721299) never reached the
# retained set: it was already removed by the depth floor and the four-crop restriction.
NON_FOLIAR_TISSUE = {
    "immature kernel", "kernel", "Kernel", "kernel (developing, in situ)",
    "seed", "seeds", "ear", "tassel", "spike",
}


def is_foliar(tissue: str) -> bool:
    """True unless the curated tissue names a non-leaf organ.

    Matches the curated string EXACTLY, deliberately. A substring test reads "ear leaf",
    "leaf (ear leaf 16 blade, midrib removed)" and "stem/leaf (rust lesion)" as non-foliar,
    and all three are leaf samples -- that mistake cost a 727-run overcount the first time.
    """
    return (tissue or "").strip() not in NON_FOLIAR_TISSUE


def is_excluded(bioproject: str) -> bool:
    """True if this BioProject is in the registry."""
    return bioproject in EXCLUDED


def reason_for(bioproject: str) -> str | None:
    """The recorded reason, or None if the BioProject is not excluded."""
    e = EXCLUDED.get(bioproject)
    return e["reason"] if e else None


def split_excluded(rows, key="BioProject"):
    """Partition an iterable of dict rows into (kept, dropped).

    Returns both halves rather than filtering in place, because a rate reported without the
    excluded samples should be reportable with them too. A caller that throws the second half
    away has made that choice explicitly.
    """
    kept, dropped = [], []
    for r in rows:
        (dropped if is_excluded(r.get(key, "")) else kept).append(r)
    return kept, dropped


def summary() -> str:
    """One line per entry, for a script to print alongside whatever it computed."""
    if not EXCLUDED:
        return "exclusions: none"
    lines = [f"exclusions: {len(EXCLUDED)} BioProjects "
             f"({sum(e['n_biosamples'] for e in EXCLUDED.values())} BioSamples at time of decision)"]
    for bp, e in sorted(EXCLUDED.items(), key=lambda kv: -kv[1]["n_biosamples"]):
        lines.append(f"  {bp}  {e['n_biosamples']:>4}  {e['kind']:<10}  {e['reason'].split('.')[0]}.")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
