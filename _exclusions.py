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
