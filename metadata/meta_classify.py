#!/usr/bin/env python3
"""
metadata/meta_classify.py — LLM BioProject classification on FULL manuscript text.

Successor to classify_metadata.py. Two structural changes:

1. Entry gate: only BioProjects with full_text (from meta_search.py + meta_text.py)
   get classified at all. classify_metadata.py classified every BioProject
   regardless of whether there was any real text to read, quietly falling back to
   abstract-only (or nothing) for most of them. samples.tsv here only ever
   contains BioSamples whose BioProject has a manuscript — matches sample_funnel
   v3's entry node.

2. Full text as PRIMARY evidence, not a 3000-char methods_text snippet (that field
   doesn't even exist anymore — meta_search.py stopped writing it 2026-08-25).
   BioSample XML metadata (tissue, geo_loc_name, host) demotes to secondary
   corroboration/fallback — the manuscript is now the trustworthy source.

Six LLM calls per BioProject — five judgment dimensions (own prompt, own rules,
own confidence+rationale each) plus one fact-extraction call for the remaining facts:

  stress            biotic|abiotic|none|unclear
  setting           field|greenhouse|growth_chamber|detached_leaf_assay|in_vitro|unclear
  tissue            aerial|non-aerial|unclear
  coinfection_intent intentional_multi_pathogen|single_pathogen_focus|not_disease_focused
  hostpath          named_pathogens, named_hosts — own call (split out of extract
                    2026-08-28), each with its own confidence+rationale. Evidence
                    matters here specifically because this field is BioProject-level
                    (see host_resolved below) — a confident, well-quoted multi-host
                    list usually means a genuine multi-host study (e.g. a rust
                    fungus's alternate-host cycle), not extraction noise.
  extract           symptom_status, exposure_type, geographic_location, library_prep,
                    host_cultivar, host_resistance (all "not_stated" if the manuscript
                    doesn't say — fact-finding, not judgment, so no per-field confidence)

named_hosts is extracted ONCE PER BIOPROJECT (from the full manuscript, which may
describe several host species across different BioSamples) and, like every other
LLM field here, gets copied onto every biosample_representative row in that
BioProject — this is correct for genuinely BioProject-constant judgments (stress,
setting, tissue) but NOT for named_hosts/named_pathogens when a BioProject spans
multiple host species across different samples (found 2026-08-28: 147 BioProjects,
1,119/6,467 samples, list llm_named_hosts as multi-valued — confirmed identical
across every BioSample in all 147, i.e. genuinely BioProject-level, not a per-sample
attribution the original call could ever have gotten right).

host_resolved fixes this with a second, targeted pass (--disambiguate-hosts):
for BioProjects where named_hosts has >1 value, one extra per-BioSample call feeds
that sample's bs_host/bs_description/tissue/isolation_source (structured BioSample
metadata, genuinely per-sample) alongside the candidate host list and asks which one
this specific sample matches, with its own confidence+rationale. Single-host
BioProjects need no extra call — host_resolved is just named_hosts[0], computed in
Python. Cache: host_disambig_cache.jsonl, keyed by BioSample (not BioProject).

Downstream computed fields (plain Python, AFTER the LLM calls — not LLM judgments,
deliberately, so the comparison logic can be revised without re-running expensive
calls as the "cryptic" definition gets refined):

  pathogen_match_status   compares llm_named_pathogens against STAT stat_pathogens
  tissue_agreement        llm_tissue vs a keyword-bucketed BioSample XML tissue value
  geo_agreement           llm_geographic_location vs BioSample XML geo_loc_name
  host_agreement          llm_host_resolved vs BioSample XML host/bs_host
  metadata_disagreement_flag   true if any of the above == "disagree" — a real
                          disagreement between manuscript and database usually means
                          something's wrong (wrong paper matched, mislabelled
                          BioSample, etc.), not that either source is "more right".
                          Caveat: a BioProject can legitimately span multiple
                          tissues/sites across its BioSamples, so a flag here isn't
                          automatically a data error — just worth a look.

Output: metadata/output/meta_classify/data/samples.tsv
        one row per biosample_representative BioSample, ONLY for full-text BioProjects

Cache: metadata/output/meta_classify/data/classify_cache.jsonl        (per-BioProject)
       metadata/output/meta_classify/data/host_disambig_cache.jsonl   (per-BioSample)
       both append-only JSONL

Requires: OPENAI_API_KEY env var

Run from crypt/:
  python metadata/meta_classify.py
  python metadata/meta_classify.py --workers 8
  python metadata/meta_classify.py --rerun-all
  python metadata/meta_classify.py --focus setting     # rerun one dimension only
  python metadata/meta_classify.py --disambiguate-hosts  # per-BioSample host pass
"""

import argparse
import csv
import json
import os
import re
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import (_Tee, link_latest, load_json, make_log_dir,
                   resolve_taxon_name, HOST_NAME_ALIASES)

try:
    from openai import OpenAI
except ImportError:
    sys.exit("ERROR: openai package not installed — pip install openai")

# ── Settings ──────────────────────────────────────────────────────────────────

BIOPROJECTS = Path("metadata/output/meta_search/data/bioprojects.json")
BIOSAMPLES  = Path("metadata/output/meta_search/data/biosamples.json")
RUNS_TSV    = Path("stat/output/stat_filter/data/runs.tsv")
PHIBASE_DB  = Path("stat/output/stat_build/data/phibase_db.json")
OUT_DIR     = Path("metadata/output/meta_classify")
CACHE_PATH  = OUT_DIR / "data" / "classify_cache.jsonl"
DISAMBIG_CACHE_PATH = OUT_DIR / "data" / "host_disambig_cache.jsonl"

MODEL       = "gpt-4o-mini"
# The real bottleneck now is tokens/min, not requests/min: each call carries up to
# ~15-20k tokens of full manuscript text (classify_metadata.py's old prompts were
# ~1-2k tokens, where request-count throttling was the right lever). At the
# observed 200,000 TPM org limit, that's only ~10-13 calls/minute sustainable
# SYSTEM-WIDE regardless of worker count — more workers just means more threads
# blocked on the same token budget, not more throughput. Keep workers modest.
MAX_WORKERS = 4
MAX_RETRIES = 6
_TPM_LIMIT  = 170_000   # safety margin under the observed 200,000 hard cap
_CHARS_PER_TOKEN = 4    # rough estimate for prompt sizing, not exact

_TPM_LOCK   = threading.Lock()
_tpm_window: list[tuple] = []   # [(timestamp, estimated_tokens), ...] within the last 60s


def _tpm_wait(estimated_tokens: int) -> None:
    """Block until sending `estimated_tokens` more would stay under _TPM_LIMIT for
    the trailing 60s window. Global across all worker threads — this is what
    per-request rate limiting (the old _RATE_MIN_GAP approach) can't do, since the
    actual OpenAI limit here is measured in tokens, not requests."""
    with _TPM_LOCK:
        while True:
            now = time.time()
            _tpm_window[:] = [(t, n) for t, n in _tpm_window if now - t < 60]
            used = sum(n for _, n in _tpm_window)
            if used + estimated_tokens <= _TPM_LIMIT:
                _tpm_window.append((now, estimated_tokens))
                return
            oldest_t = _tpm_window[0][0] if _tpm_window else now
            wait_for = max(0.5, 60 - (now - oldest_t) + 0.25)
            time.sleep(wait_for)


def _estimate_tokens(prompt: str, max_tokens: int) -> int:
    return len(prompt) // _CHARS_PER_TOKEN + max_tokens

FULL_TEXT_CAP = 60_000   # meta_text.py's own cap; re-guarded here defensively

_FOCUS_VALUES = {"all", "stress", "setting", "tissue", "coinfection", "hostpath", "extract"}

_STRESS_VALUES  = {"biotic", "abiotic", "none", "unclear"}
_SETTING_VALUES = {"field", "greenhouse", "growth_chamber", "detached_leaf_assay",
                   "in_vitro", "unclear"}
_TISSUE_VALUES  = {"aerial", "non-aerial", "unclear"}
_COINF_VALUES   = {"intentional_multi_pathogen", "single_pathogen_focus",
                   "not_disease_focused", "unclear"}
_CONF_VALUES    = {"high", "medium", "low"}

_SYMPTOM_VALUES     = {"symptomatic", "asymptomatic", "mixed", "not_stated"}
_EXPOSURE_VALUES    = {"natural", "deliberate_inoculation", "not_stated"}
_LIBPREP_VALUES     = {"polyA_selection", "rRNA_depletion", "total_RNA", "not_stated"}
_RESISTANCE_VALUES  = {"resistant", "susceptible", "tolerant", "not_stated"}

_SAMPLE_FIELDS = [
    # BioSample identity
    "BioSample", "BioProject", "mode", "biosample_n_runs",
    # STAT
    "stat_host", "host_pct", "stat_pathogens", "interaction_status",
    "co_infection_flag", "same_genus_secondary", "n_pathogens",
    "fungi_pct", "oomycete_pct", "nematode_pct",
    # BioSample XML
    "tissue", "geo_loc_name", "collection_date", "dev_stage",
    "isolation_source", "lat_lon", "bs_host",
    # LLM — stress
    "llm_stress", "llm_stress_confidence", "llm_stress_rationale",
    # LLM — setting
    "llm_study_setting", "llm_setting_confidence", "llm_setting_rationale",
    # LLM — tissue
    "llm_tissue", "llm_tissue_confidence", "llm_tissue_rationale",
    # LLM — coinfection intent
    "llm_coinfection_intent", "llm_coinfection_confidence", "llm_coinfection_rationale",
    # LLM — named hosts/pathogens (own call, own confidence — see module docstring)
    "llm_named_pathogens", "llm_pathogens_confidence", "llm_pathogens_rationale",
    "llm_named_hosts", "llm_hosts_confidence", "llm_hosts_rationale",
    # host_resolved: single-value disambiguation of llm_named_hosts per BioSample —
    # equals named_hosts[0] directly (no LLM call) when named_hosts has one value;
    # for multi-value BioProjects, only populated after --disambiguate-hosts
    "llm_host_resolved", "llm_host_resolved_confidence", "llm_host_resolved_rationale",
    # NCBI taxid resolution (plain Python, deterministic lookup against
    # phibase_db.json + a common-name alias table — NOT asked of the LLM, which
    # would invite confident-looking wrong numbers; see _util.resolve_taxon_name)
    "llm_named_pathogens_taxids", "llm_named_hosts_taxids", "llm_host_resolved_taxid",
    # LLM — extraction (remaining facts, no per-field confidence — see module docstring)
    "llm_symptom_status", "llm_exposure_type",
    "llm_geographic_location", "llm_library_prep", "llm_host_cultivar", "llm_host_resistance",
    # Downstream computed (Python, not LLM — see module docstring)
    "pathogen_match_status", "tissue_agreement", "geo_agreement", "host_agreement",
    "metadata_disagreement_flag",
    # Literature
    "title", "pmid", "doi", "pmid_source", "submission_date", "pub_date",
]

_MODE_DESC = {
    "mal":   "MAL (microbe-as-library): library source organism is a PLANT PATHOGEN — "
             "reads are primarily from the pathogen.",
    "hal":   "HAL (host-as-library): library source organism is a PLANT HOST — "
             "reads are primarily from the host transcriptome.",
    "mixed": "Mixed MAL+HAL: BioProject contains both pathogen-library and host-library runs.",
}

# ── Prompt templates ───────────────────────────────────────────────────────────
#
# Shared manuscript context goes FIRST in every prompt (fixed prefix across all 5
# calls for the same BP), dimension-specific instructions LAST — this lets
# OpenAI's automatic prompt caching discount the repeated resend of the same long
# text across calls issued close together in time.

_CTX = """\
--- BioProject: {bp} ---
Library mode: {mode_desc}
Title: {title}
Description: {description}
Abstract: {abstract}

Full manuscript text (primary source of evidence — trust this over BioSample
metadata below, which is secondary/fallback only):
{full_text}
"""

_PROMPT_STRESS = """\
{ctx}
---
You are classifying the STRESS TYPE this RNA-seq BioProject investigated, using the
full manuscript text above. Pick ONE:

biotic   — exposure to any pathogen, pest, or parasite; deliberate inoculation or
           naturally occurring infection; disease resistance or susceptibility trials;
           any study where a named pathogen is the biological focus; co-infection designs
abiotic  — physical or chemical stress ONLY: drought, heat, cold, salt, flooding, UV,
           wounding, heavy metals, nutrient deficiency; NO pathogen involved
none     — no stress: pure host biology, developmental atlas, genotype comparison,
           QTL/GWAS transcriptomics, phenology; neither pathogen nor abiotic treatment
unclear  — the manuscript genuinely never states a stress/treatment anywhere — check
           the full text before defaulting here, don't use this just because it
           wasn't in the abstract

RULES:
1. MAL mode: library organism IS the pathogen → almost certainly biotic
2. Named pathogen anywhere in the manuscript → biotic
3. Hormone treatments (jasmonate, salicylate, ethylene) WITHOUT a pathogen → abiotic
4. Trust an explicit statement in Methods over the title/abstract alone

Return ONLY valid JSON, no surrounding text:
{{"stress": "<biotic|abiotic|none|unclear>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences, quote a short phrase from the manuscript>"}}"""


_PROMPT_SETTING = """\
{ctx}
--- BioSample collection signals ({n_biosamples} BioSamples, secondary/fallback only) ---
XML geo_loc_name (countries):        {geolocnames}
XML collection_date present:         {n_with_collection_date} of {n_biosamples} BioSamples
BioSample descriptions:              {bs_descriptions}

---
You are classifying the STUDY SETTING of this RNA-seq BioProject, using the full
manuscript text above as your PRIMARY source of evidence. BioSample metadata is
secondary — use it only to corroborate or as a fallback when the manuscript is silent.

Where were the sequenced plant samples grown/collected, and was wild pathogen exposure
possible? Search the full Materials and Methods before deciding. Pick ONE:

field                — natural or agricultural conditions with unrestricted pathogen
                       exposure: farm, orchard, commercial crop, natural population,
                       field survey, disease surveillance. Includes field trials that
                       were deliberately inoculated, AS LONG AS the growing environment
                       itself was open (see exposure_type in the separate extraction
                       call for the inoculation distinction — don't fold it in here).
greenhouse           — greenhouse-grown, polytunnel, screenhouse — semi-controlled but
                       NOT sealed against unplanned/wild pathogen entry.
growth_chamber       — sealed growth chamber, controlled-environment room, phytotron,
                       climate cabinet — pathogen exposure is architecturally excluded.
detached_leaf_assay  — excised leaf/tissue inoculated ex planta (petri dish, leaf disk).
in_vitro             — axenic culture, cell suspension, callus, protoplast, or the
                       pathogen/microbe grown alone with NO host tissue present at all.
unclear              — the manuscript genuinely never states growing conditions
                       anywhere in Methods/Materials. Don't use this just because the
                       signal isn't in the abstract — check the full text first.

RULES:
1. Trust an explicit statement in Methods over any BioSample XML signal or species prior.
2. "Field isolate" describing only where the PATHOGEN STRAIN originated (later used in
   a lab/greenhouse/growth-chamber inoculation) does NOT make the experiment field.
3. Arabidopsis/Nicotiana benthamiana/Brachypodium are NOT automatically growth_chamber —
   check the text; only fall back on this prior if the manuscript is truly silent.
4. If Methods distinguishes field vs greenhouse for different sample subsets, pick
   whichever the SEQUENCED samples (not just any samples in the paper) came from.

Return ONLY valid JSON, no surrounding text:
{{"study_setting": "<field|greenhouse|growth_chamber|detached_leaf_assay|in_vitro|unclear>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences, must quote a short phrase from the manuscript if not unclear>"}}"""


_PROMPT_TISSUE = """\
{ctx}
--- BioSample tissue signals ({n_biosamples} BioSamples, secondary/fallback only) ---
XML tissue values:         {tissues}
XML dev_stage values:      {dev_stages}
BioSample descriptions:    {bs_descriptions}

---
You are classifying the TISSUE TYPE sequenced in this RNA-seq BioProject, using the
full manuscript text above as your PRIMARY source. Pick ONE:

aerial      — above-ground tissue: leaf, leaflet, blade, stem, shoot, hypocotyl,
              flower, sepal, petal, bud, inflorescence, spike, panicle, rachis,
              seed, grain, endosperm, embryo, fruit, pericarp, pollen, anther,
              petiole, canopy, pod, flag leaf, urediniospore (rust spores)
non-aerial  — below-ground or storage tissue: root, lateral root, root tip, tuber,
              bulb, corm, rhizome, crown, stolon
unclear     — whole plant / seedling without tissue specification; mixed shoot+root;
              callus, protoplast, suspension culture; the manuscript genuinely never
              specifies — check the full Methods before defaulting here

RULES:
1. Trust an explicit tissue statement in Methods over XML tissue values or title.
2. "Leaf rust", "stripe rust on leaf", "flag leaf", "infected leaf" → aerial
3. "Root rot", "crown rot", "Phytophthora root infection" → non-aerial
4. MAL mode: tissue is where the pathogen resides — rust/powdery mildew → aerial;
   Phytophthora/Pythium root/crown rot → non-aerial
5. Seeds and grain are aerial (harvested above ground)

Return ONLY valid JSON, no surrounding text:
{{"tissue": "<aerial|non-aerial|unclear>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences, quote a short phrase from the manuscript if not unclear>"}}"""


_PROMPT_COINFECTION = """\
{ctx}
STAT-detected pathogens in sequencing data: {stat_pathogen_summary}

---
You are determining whether this RNA-seq study was DESIGNED to investigate multiple
pathogens at once, using the full manuscript text above. Pick ONE:

intentional_multi_pathogen — the study explicitly investigates two or more pathogens/
    pests together as its design (co-infection experiment, pathogen interaction study,
    disease complex, sequential or simultaneous inoculation with multiple agents).
single_pathogen_focus      — the study names and investigates ONE pathogen as its
    focus, even if other organisms are mentioned in passing (background microbiome,
    prior literature, biocontrol agents not applied to these samples).
not_disease_focused        — the study is not about pathogen infection at all (abiotic
    stress, developmental biology, genotype comparison, etc.) — no pathogen is the
    intended subject.

RULES:
1. A pathogen appearing only in Introduction/Discussion as background/citation does
   NOT count — only pathogens actually applied to or investigated in THESE samples.
2. Biocontrol/beneficial microbe studies count as intentional_multi_pathogen only if
   disease-causing organisms are ALSO part of the design, not beneficial-microbe-alone.
3. This is independent of what STAT actually detected — STAT pathogens are shown only
   as context; classify the AUTHORS' STATED design intent, not the sequencing result.

Return ONLY valid JSON, no surrounding text:
{{"coinfection_intent": "<intentional_multi_pathogen|single_pathogen_focus|not_disease_focused>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences quoting the manuscript>"}}"""


_PROMPT_HOSTPATHOGEN = """\
{ctx}
---
You are extracting the named PATHOGENS and HOSTS stated in this manuscript. Only
report what the text actually says — empty list if the manuscript doesn't name any.
Do not infer or guess from the species prior of the mode/BioSample metadata.

named_pathogens — list of pathogen species/strains the authors explicitly name as
    present in or relevant to THESE samples (not just cited from other studies,
    not background/biocontrol organisms). Empty list if none named.
named_hosts     — list of host plant SPECIES explicitly named for these samples.
    A manuscript studying a pathogen's full host range or life cycle (e.g. a rust
    fungus alternating between a crop and a wild reservoir host) may genuinely name
    several — list all of them; do not collapse to one. This field is extracted
    once for the whole BioProject and may cover more than one BioSample's actual
    host — a downstream pass disambiguates per BioSample, so completeness here
    matters more than picking a single "best" one.

    Formatting rules for BOTH lists:
    - Prefer the scientific (binomial) name over a common name, e.g. "Triticum
      aestivum" not "wheat" — use the common name only if the manuscript never
      gives a scientific name anywhere for that organism.
    - One entry per distinct species. If the manuscript uses both a common name
      and the scientific name for the SAME organism, that's still ONE entry
      (the scientific name) — do not list it twice.
    - Species only — do not include cultivar/genotype/accession names as
      separate list entries (e.g. "wheat cultivar Fielder" is the species
      Triticum aestivum; the cultivar name belongs in host_cultivar, a
      separate field extracted elsewhere, not here).

For EACH of the two fields, also give ONE confidence + rationale for the field as
a whole (not per list item): how well-supported is the list you extracted, and
why — quote a short phrase.

Return ONLY valid JSON, no surrounding text:
{{"named_pathogens": ["..."],
  "pathogens_confidence": "<high|medium|low>",
  "pathogens_rationale": "<1-2 sentences, quote a short phrase from the manuscript>",
  "named_hosts": ["..."],
  "hosts_confidence": "<high|medium|low>",
  "hosts_rationale": "<1-2 sentences, quote a short phrase from the manuscript>"}}"""


_PROMPT_EXTRACT = """\
{ctx}
---
You are extracting FACTS stated in this manuscript. Only report what the text
actually says — use "not_stated" for anything the manuscript doesn't mention. Do
not infer or guess.

symptom_status       — "symptomatic" (visible disease signs described), "asymptomatic"
                       (explicitly healthy/no symptoms), "mixed" (both sampled), or
                       "not_stated".
exposure_type        — "natural" (spontaneous/wild infection, disease survey),
                       "deliberate_inoculation" (authors inoculated with a specific
                       pathogen, regardless of field/greenhouse/lab setting), or
                       "not_stated".
geographic_location  — country/region/site as described in Methods (e.g. "Punjab,
                       India"), or "not_stated".
library_prep         — "polyA_selection", "rRNA_depletion", "total_RNA", or "not_stated".
host_cultivar        — specific cultivar/genotype/accession name if given, else
                       "not_stated".
host_resistance      — "resistant", "susceptible", "tolerant", or "not_stated" — as
                       DESCRIBED BY THE AUTHORS for the cultivar used, not inferred.

Return ONLY valid JSON, no surrounding text:
{{"symptom_status": "<symptomatic|asymptomatic|mixed|not_stated>",
  "exposure_type": "<natural|deliberate_inoculation|not_stated>",
  "geographic_location": "<text or not_stated>",
  "library_prep": "<polyA_selection|rRNA_depletion|total_RNA|not_stated>",
  "host_cultivar": "<text or not_stated>",
  "host_resistance": "<resistant|susceptible|tolerant|not_stated>"}}"""


_PROMPT_HOST_DISAMBIGUATE = """\
--- BioProject: {bp} ---
This BioProject's manuscript names multiple candidate host species (extracted
separately, from the full text): {candidates}

You are deciding which ONE of these candidates the specific BioSample below is —
using its structured BioSample metadata (from the SRA/ENA record, NOT the
manuscript). This is per-sample metadata, not manuscript text, so treat it as
direct evidence of what was actually sequenced for this sample:

BioSample:          {biosample}
XML host field:      {bs_host}
BioSample description: {bs_description}
Tissue:              {tissue}
Isolation source:     {isolation_source}
Geographic location:  {geo_loc_name}

Pick the ONE candidate from the list above that this BioSample's metadata best
supports. If the metadata doesn't clearly point to any single candidate (e.g. it's
empty, or it's ambiguous/generic), say so — do not guess.

Return ONLY valid JSON, no surrounding text:
{{"host_resolved": "<one of the candidates, exactly as given, or 'unresolved'>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences citing which BioSample field supports the pick>"}}"""


# ── BioSample summary ──────────────────────────────────────────────────────────

def _bp_summary(rows: list[dict], bs_meta: dict) -> dict:
    tiss    = sorted({(bs_meta.get(r["BioSample"], {}).get("tissue", "") or "")
                      for r in rows} - {""})
    geo     = sorted({(bs_meta.get(r["BioSample"], {}).get("geo_loc_name", "") or "")
                      .split(":")[0].strip()
                      for r in rows} - {""})
    descs   = sorted({(bs_meta.get(r["BioSample"], {}).get("bs_description", "") or "")
                      for r in rows} - {""})
    devs    = sorted({(bs_meta.get(r["BioSample"], {}).get("dev_stage", "") or "")
                      for r in rows} - {""})
    n_cdates = sum(1 for r in rows
                   if (bs_meta.get(r["BioSample"], {}).get("collection_date", "") or "").strip())
    stat_pats: Counter = Counter()
    for r in rows:
        for p in (r.get("stat_pathogens", "") or "").split(";"):
            p = p.strip()
            if p:
                stat_pats[p] += 1
    modes = {r.get("mode", "") for r in rows}
    mode  = "hal" if modes == {"hal"} else "mal" if modes == {"mal"} else "mixed"
    return {
        "n_biosamples":           len(rows),
        "tissues":                ", ".join(tiss[:10])  or "none in XML",
        "geolocnames":            ", ".join(geo[:8])    or "none in XML",
        "bs_descriptions":        ", ".join(descs[:5])  or "none in XML",
        "dev_stages":             ", ".join(devs[:5])   or "none in XML",
        "n_with_collection_date": n_cdates,
        "stat_pathogen_summary":  ", ".join(p for p, _ in stat_pats.most_common(5))
                                  or "none detected by STAT",
        "mode": mode,
    }


# ── Prompt builders ────────────────────────────────────────────────────────────

def _ctx(bp: str, summary: dict, meta: dict) -> str:
    return _CTX.format(
        bp          = bp,
        mode_desc   = _MODE_DESC.get(summary["mode"], ""),
        title       = (meta.get("title",       "") or "")[:300],
        description = (meta.get("description", "") or "")[:500],
        abstract    = (meta.get("abstract",    "") or "")[:3000],
        full_text   = (meta.get("full_text",   "") or "")[:FULL_TEXT_CAP],
    )


def _stress_prompt(bp, summary, meta):
    return _PROMPT_STRESS.format(ctx=_ctx(bp, summary, meta))


def _setting_prompt(bp, summary, meta):
    return _PROMPT_SETTING.format(
        ctx=_ctx(bp, summary, meta),
        n_biosamples=summary["n_biosamples"], geolocnames=summary["geolocnames"],
        n_with_collection_date=summary["n_with_collection_date"],
        bs_descriptions=summary["bs_descriptions"],
    )


def _tissue_prompt(bp, summary, meta):
    return _PROMPT_TISSUE.format(
        ctx=_ctx(bp, summary, meta),
        n_biosamples=summary["n_biosamples"], tissues=summary["tissues"],
        dev_stages=summary["dev_stages"], bs_descriptions=summary["bs_descriptions"],
    )


def _coinfection_prompt(bp, summary, meta):
    return _PROMPT_COINFECTION.format(
        ctx=_ctx(bp, summary, meta), stat_pathogen_summary=summary["stat_pathogen_summary"],
    )


def _hostpath_prompt(bp, summary, meta):
    return _PROMPT_HOSTPATHOGEN.format(ctx=_ctx(bp, summary, meta))


def _extract_prompt(bp, summary, meta):
    return _PROMPT_EXTRACT.format(ctx=_ctx(bp, summary, meta))


def _host_disambiguate_prompt(bp, candidates, bsid, bs_xml):
    return _PROMPT_HOST_DISAMBIGUATE.format(
        bp=bp, candidates="; ".join(candidates), biosample=bsid,
        bs_host=bs_xml.get("host", "") or "not stated",
        bs_description=bs_xml.get("bs_description", "") or "not stated",
        tissue=bs_xml.get("tissue", "") or "not stated",
        isolation_source=bs_xml.get("isolation_source", "") or "not stated",
        geo_loc_name=bs_xml.get("geo_loc_name", "") or "not stated",
    )


# ── API call helpers ───────────────────────────────────────────────────────────

def _parse_json(text: str) -> dict:
    try:
        return json.loads(text.strip())
    except json.JSONDecodeError:
        pass
    m = re.search(r'\{.*\}', text, re.DOTALL)
    if m:
        try:
            return json.loads(m.group())
        except json.JSONDecodeError:
            pass
    return {}


def _call_dim(prompt: str, client: "OpenAI", bp: str, dim: str, max_tokens: int = 350) -> dict | None:
    """Returns the parsed dict, or None if the call/parse never produced anything
    usable after all retries. None is a distinct failure signal — callers must NOT
    fold it into "unclear"/"not_stated", or a real API failure becomes
    indistinguishable from a genuine "the paper doesn't say" classification."""
    last_exc = None
    est_tokens = _estimate_tokens(prompt, max_tokens)
    for attempt in range(MAX_RETRIES):
        _tpm_wait(est_tokens)
        try:
            resp = client.chat.completions.create(
                model      = MODEL,
                max_tokens = max_tokens,
                messages   = [{"role": "user", "content": prompt}],
            )
            parsed = _parse_json(resp.choices[0].message.content)
            if parsed:
                return parsed
            last_exc = "empty/unparseable JSON response"
        except Exception as exc:
            last_exc = exc
        if attempt < MAX_RETRIES - 1:
            time.sleep(min(2 ** attempt, 60))
    print(f"  FAILED {bp} [{dim}] after {MAX_RETRIES} attempts: {last_exc}", flush=True)
    return None


def _classify_bp(bp: str, summary: dict, meta: dict, client: "OpenAI", focus: str,
                 cached: dict | None = None) -> dict:
    """cached: the existing cache entry (if any) — dimensions already _is_ok in it
    are skipped rather than re-run, so a retry pass only pays for what failed."""
    result: dict = {"bp": bp}
    conf = _CONF_VALUES
    cached = cached or {}

    if focus in ("all", "stress") and not _is_ok(cached, "llm_stress"):
        raw = _call_dim(_stress_prompt(bp, summary, meta), client, bp, "stress")
        if raw is None:
            result["llm_stress"] = "api_error"
            result["llm_stress_confidence"] = ""
            result["llm_stress_rationale"] = ""
        else:
            stress = raw.get("stress", "")
            result["llm_stress"]            = stress if stress in _STRESS_VALUES else "unclear"
            result["llm_stress_confidence"] = raw.get("confidence", "") if raw.get("confidence") in conf else "low"
            result["llm_stress_rationale"]  = str(raw.get("rationale", "") or "")

    if focus in ("all", "setting") and not _is_ok(cached, "llm_study_setting"):
        raw = _call_dim(_setting_prompt(bp, summary, meta), client, bp, "setting")
        if raw is None:
            result["llm_study_setting"] = "api_error"
            result["llm_setting_confidence"] = ""
            result["llm_setting_rationale"] = ""
        else:
            setting = raw.get("study_setting", "")
            result["llm_study_setting"]      = setting if setting in _SETTING_VALUES else "unclear"
            result["llm_setting_confidence"] = raw.get("confidence", "") if raw.get("confidence") in conf else "low"
            result["llm_setting_rationale"]  = str(raw.get("rationale", "") or "")

    if focus in ("all", "tissue") and not _is_ok(cached, "llm_tissue"):
        raw = _call_dim(_tissue_prompt(bp, summary, meta), client, bp, "tissue")
        if raw is None:
            result["llm_tissue"] = "api_error"
            result["llm_tissue_confidence"] = ""
            result["llm_tissue_rationale"] = ""
        else:
            tissue = raw.get("tissue", "")
            result["llm_tissue"]            = tissue if tissue in _TISSUE_VALUES else "unclear"
            result["llm_tissue_confidence"] = raw.get("confidence", "") if raw.get("confidence") in conf else "low"
            result["llm_tissue_rationale"]  = str(raw.get("rationale", "") or "")

    if focus in ("all", "coinfection") and not _is_ok(cached, "llm_coinfection_intent"):
        raw = _call_dim(_coinfection_prompt(bp, summary, meta), client, bp, "coinfection")
        if raw is None:
            result["llm_coinfection_intent"] = "api_error"
            result["llm_coinfection_confidence"] = ""
            result["llm_coinfection_rationale"] = ""
        else:
            ci = raw.get("coinfection_intent", "")
            result["llm_coinfection_intent"]     = ci if ci in _COINF_VALUES else "unclear"
            result["llm_coinfection_confidence"] = raw.get("confidence", "") if raw.get("confidence") in conf else "low"
            result["llm_coinfection_rationale"]  = str(raw.get("rationale", "") or "")

    if focus in ("all", "hostpath") and not _complete(cached, "hostpath"):
        raw = _call_dim(_hostpath_prompt(bp, summary, meta), client, bp, "hostpath", max_tokens=500)
        if raw is None:
            result["llm_named_pathogens"] = ["api_error"]
            result["llm_pathogens_confidence"] = ""
            result["llm_pathogens_rationale"] = ""
            result["llm_named_hosts"] = ["api_error"]
            result["llm_hosts_confidence"] = ""
            result["llm_hosts_rationale"] = ""
        else:
            pathogens = raw.get("named_pathogens", [])
            hosts     = raw.get("named_hosts", [])
            result["llm_named_pathogens"] = [str(p) for p in pathogens] if isinstance(pathogens, list) else []
            result["llm_pathogens_confidence"] = raw.get("pathogens_confidence", "") if raw.get("pathogens_confidence") in conf else "low"
            result["llm_pathogens_rationale"]  = str(raw.get("pathogens_rationale", "") or "")
            result["llm_named_hosts"]     = [str(h) for h in hosts]     if isinstance(hosts, list)     else []
            result["llm_hosts_confidence"] = raw.get("hosts_confidence", "") if raw.get("hosts_confidence") in conf else "low"
            result["llm_hosts_rationale"]  = str(raw.get("hosts_rationale", "") or "")

    if focus in ("all", "extract") and not _is_ok(cached, "llm_symptom_status"):
        raw = _call_dim(_extract_prompt(bp, summary, meta), client, bp, "extract", max_tokens=350)
        if raw is None:
            result["llm_symptom_status"] = "api_error"
            result["llm_exposure_type"] = "api_error"
            result["llm_geographic_location"] = "api_error"
            result["llm_library_prep"] = "api_error"
            result["llm_host_cultivar"] = "api_error"
            result["llm_host_resistance"] = "api_error"
        else:
            sym = raw.get("symptom_status", "")
            result["llm_symptom_status"] = sym if sym in _SYMPTOM_VALUES else "not_stated"
            exp = raw.get("exposure_type", "")
            result["llm_exposure_type"] = exp if exp in _EXPOSURE_VALUES else "not_stated"
            result["llm_geographic_location"] = str(raw.get("geographic_location", "") or "not_stated")
            lp = raw.get("library_prep", "")
            result["llm_library_prep"] = lp if lp in _LIBPREP_VALUES else "not_stated"
            result["llm_host_cultivar"] = str(raw.get("host_cultivar", "") or "not_stated")
            res = raw.get("host_resistance", "")
            result["llm_host_resistance"] = res if res in _RESISTANCE_VALUES else "not_stated"

    return result


def _disambiguate_host(bp: str, candidates: list, bsid: str, bs_xml: dict,
                       client: "OpenAI") -> dict:
    """Per-BioSample host disambiguation — only called for BioProjects where
    named_hosts has >1 value (see module docstring). Returns
    {host_resolved, confidence, rationale}; host_resolved is one of `candidates`
    verbatim, or 'unresolved' if the BioSample metadata doesn't clearly point to
    one. None fields on API failure (caller tags api_error, same convention as
    _call_dim)."""
    raw = _call_dim(_host_disambiguate_prompt(bp, candidates, bsid, bs_xml),
                    client, bp, f"host_disambig:{bsid}", max_tokens=200)
    if raw is None:
        return {"host_resolved": "api_error", "confidence": "", "rationale": ""}
    picked = str(raw.get("host_resolved", "") or "")
    if picked not in candidates:
        picked = "unresolved"
    conf = raw.get("confidence", "")
    return {
        "host_resolved": picked,
        "confidence": conf if conf in _CONF_VALUES else "low",
        "rationale": str(raw.get("rationale", "") or ""),
    }


# ── Downstream computed fields (plain Python, not LLM) ─────────────────────────

_AERIAL_KW = ("leaf", "leaflet", "blade", "stem", "shoot", "hypocotyl", "flower",
              "sepal", "petal", "bud", "inflorescence", "spike", "panicle", "rachis",
              "seed", "grain", "endosperm", "embryo", "fruit", "pericarp", "pollen",
              "anther", "petiole", "canopy", "pod", "flag leaf", "urediniospore")
_NONAERIAL_KW = ("root", "tuber", "bulb", "corm", "rhizome", "crown", "stolon")


def _bucket_xml_tissue(xml_tissue: str) -> str:
    t = (xml_tissue or "").lower()
    if not t:
        return ""
    if any(k in t for k in _NONAERIAL_KW):
        return "non-aerial"
    if any(k in t for k in _AERIAL_KW):
        return "aerial"
    return ""


def _tissue_agreement(llm_tissue: str, xml_tissue: str) -> str:
    bucket = _bucket_xml_tissue(xml_tissue)
    if not bucket or llm_tissue not in ("aerial", "non-aerial"):
        return "not_comparable"
    return "agree" if bucket == llm_tissue else "disagree"


_STOPWORDS = {"the", "of", "and", "in", "near", "at", "a", "an"}


def _geo_agreement(llm_geo: str, xml_geo: str) -> str:
    if not xml_geo or not xml_geo.strip():
        return "not_comparable"
    if not llm_geo or llm_geo == "not_stated":
        return "not_comparable"
    llm_tokens = {w for w in re.findall(r"[a-z]+", llm_geo.lower()) if w not in _STOPWORDS}
    xml_tokens = {w for w in re.findall(r"[a-z]+", xml_geo.lower()) if w not in _STOPWORDS}
    if not llm_tokens or not xml_tokens:
        return "not_comparable"
    return "agree" if llm_tokens & xml_tokens else "disagree"


def _host_agreement(llm_hosts: list, xml_host: str) -> str:
    if not xml_host or not xml_host.strip():
        return "not_comparable"
    if not llm_hosts:
        return "not_comparable"
    xml_l = xml_host.lower().strip()
    xml_genus = xml_l.split()[0] if xml_l.split() else xml_l
    for h in llm_hosts:
        h_l = (h or "").lower().strip()
        if not h_l:
            continue
        h_genus = h_l.split()[0] if h_l.split() else h_l
        if h_l in xml_l or xml_l in h_l or h_genus == xml_genus:
            return "agree"
    return "disagree"


def _pathogen_match_status(llm_pathogens: list, stat_pathogens_str: str) -> str:
    stat_list = [p.strip() for p in (stat_pathogens_str or "").split(";") if p.strip()]
    named = [p for p in (llm_pathogens or []) if p and p != "not_stated"]
    if not named and not stat_list:
        return "no_pathogens_either_source"
    if not stat_list:
        return "named_not_detected_by_stat"
    if not named:
        return "stat_only_no_named"          # STAT found pathogen(s), authors named none
    llm_genera  = {p.split()[0].lower() for p in named if p.split()}
    stat_genera = {p.split()[0].lower() for p in stat_list if p.split()}
    overlap    = llm_genera & stat_genera
    extra_stat = stat_genera - llm_genera
    if overlap and not extra_stat:
        return "match"
    if overlap and extra_stat:
        return "partial_match_plus_undeclared"
    if not overlap and extra_stat:
        return "no_match_stat_found_different"
    return "no_match"


# ── Completeness check ─────────────────────────────────────────────────────────

_ALL_REQUIRED = {
    "llm_stress", "llm_stress_confidence",
    "llm_study_setting", "llm_setting_confidence",
    "llm_tissue", "llm_tissue_confidence",
    "llm_coinfection_intent", "llm_coinfection_confidence",
    "llm_named_pathogens", "llm_pathogens_confidence",
    "llm_named_hosts", "llm_hosts_confidence",
    "llm_symptom_status",
}
_DIM_REQUIRED = {
    "stress":      {"llm_stress",             "llm_stress_confidence"},
    "setting":     {"llm_study_setting",       "llm_setting_confidence"},
    "tissue":      {"llm_tissue",              "llm_tissue_confidence"},
    "coinfection": {"llm_coinfection_intent",  "llm_coinfection_confidence"},
    # llm_pathogens_confidence/llm_hosts_confidence are new fields (added when
    # named_pathogens/named_hosts were split into their own call) — requiring
    # them here means old cache entries (pre-split, no confidence fields) are
    # correctly seen as incomplete and get the hostpath call re-run, rather
    # than being silently treated as already-done.
    "hostpath":    {"llm_named_pathogens", "llm_pathogens_confidence",
                     "llm_named_hosts", "llm_hosts_confidence"},
    "extract":     {"llm_symptom_status"},
}
_VAL_KEY = {
    "stress": "llm_stress", "setting": "llm_study_setting",
    "tissue": "llm_tissue", "coinfection": "llm_coinfection_intent",
    "hostpath": "llm_named_pathogens", "extract": "llm_symptom_status",
}


def _is_ok(entry: dict, val_key: str) -> bool:
    val = entry.get(val_key)
    if val is None or val == "":
        return False
    if val == "api_error" or val == ["api_error"]:
        return False   # a failed call must be retried, not treated as cached-complete
    return True


def _complete(entry: dict, focus: str = "all") -> bool:
    fields = _ALL_REQUIRED if focus == "all" else _DIM_REQUIRED[focus]
    if not fields.issubset(entry.keys()):
        return False
    if focus == "all":
        return all(_is_ok(entry, k) for k in _VAL_KEY.values())
    return _is_ok(entry, _VAL_KEY[focus])


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers",   type=int, default=MAX_WORKERS,
                    help=f"parallel BioProject workers (default {MAX_WORKERS})")
    ap.add_argument("--rerun-all", action="store_true",
                    help="re-classify all BioProjects, ignoring cache")
    ap.add_argument("--focus",     default="all", choices=sorted(_FOCUS_VALUES),
                    help="classify only one dimension (updates that dimension in cache)")
    ap.add_argument("--disambiguate-hosts", action="store_true",
                    help="after BioProject-level classification, run the per-BioSample "
                         "host disambiguation pass for multi-value named_hosts "
                         "BioProjects (see module docstring)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    sys.stdout = _Tee(log_dir / "classify.log")
    link_latest(logs_base, log_dir / "classify.log")

    api_key = os.environ.get("OPENAI_API_KEY", "")
    if not api_key:
        sys.exit("ERROR: OPENAI_API_KEY not set")
    client = OpenAI(api_key=api_key)
    print(f"Model: {MODEL}  workers: {args.workers}  focus: {args.focus}", flush=True)

    for path, hint in [
        (BIOPROJECTS, "run metadata/meta_search.py first"),
        (BIOSAMPLES,  "run metadata/meta_search.py first"),
        (RUNS_TSV,    "run stat/stat_filter.py first"),
        (PHIBASE_DB,  "run stat/stat_build.py first"),
    ]:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found — {hint}")

    # ── Load inputs ───────────────────────────────────────────────────────────
    bp_meta: dict[str, dict] = load_json(BIOPROJECTS)
    bs_meta: dict[str, dict] = load_json(BIOSAMPLES)
    name_to_taxid: dict = load_json(PHIBASE_DB).get("name_to_taxid", {})
    host_aliases = {k.lower(): v for k, v in HOST_NAME_ALIASES.items()}

    with open(RUNS_TSV, newline="") as f:
        all_runs = list(csv.DictReader(f, delimiter="\t"))

    rep_rows = [r for r in all_runs if r.get("biosample_representative") == "True"]
    bp_to_rep: dict[str, list[dict]] = defaultdict(list)
    for r in rep_rows:
        bp_to_rep[r["BioProject"]].append(r)

    # Entry gate: only BioProjects with full manuscript text — see module docstring
    n_before = len(bp_to_rep)
    bp_to_rep = {bp: rows for bp, rows in bp_to_rep.items()
                 if (bp_meta.get(bp, {}) or {}).get("full_text")}
    print(f"BioProjects: {n_before:,} total  →  {len(bp_to_rep):,} with full text "
          f"(entry gate)  representative BioSamples: "
          f"{sum(len(v) for v in bp_to_rep.values()):,}", flush=True)

    # ── Load cache ────────────────────────────────────────────────────────────
    cache: dict[str, dict] = {}
    if CACHE_PATH.exists() and not args.rerun_all:
        with open(CACHE_PATH) as f:
            for line in f:
                line = line.strip()
                if line:
                    rec = json.loads(line)
                    cache[rec["bp"]] = rec
    print(f"Cache: {len(cache):,} entries", flush=True)

    # ── Classify ──────────────────────────────────────────────────────────────
    todo   = [bp for bp in bp_to_rep
              if args.rerun_all or not _complete(cache.get(bp, {}), args.focus)]
    n_skip = len(bp_to_rep) - len(todo)
    print(f"To classify: {len(todo):,}  skipping cached: {n_skip:,}", flush=True)

    done = 0
    dim_fail_counts: Counter = Counter()   # per-dimension api_error tally, for visibility
    if todo:
        with open(CACHE_PATH, "a") as cache_fh, \
             ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(
                    _classify_bp, bp,
                    _bp_summary(bp_to_rep[bp], bs_meta),
                    bp_meta.get(bp, {}),
                    client,
                    args.focus,
                    cache.get(bp, {}),
                ): bp
                for bp in todo
            }
            for fut in as_completed(futures):
                bp     = futures[fut]
                result = fut.result()
                merged = dict(cache.get(bp, {"bp": bp}))
                merged.update(result)
                cache[bp] = merged
                cache_fh.write(json.dumps(merged) + "\n")
                cache_fh.flush()
                done += 1

                bp_errors = [dim for dim, key in _VAL_KEY.items()
                            if result.get(key) in ("api_error", ["api_error"])]
                for dim in bp_errors:
                    dim_fail_counts[dim] += 1
                if bp_errors:
                    print(f"  [{done:,}/{len(todo)}] {bp}: FAILED dims = {bp_errors}", flush=True)

                if done % 20 == 0 or done == len(todo):
                    total_fail = sum(dim_fail_counts.values())
                    print(f"  {done:,} / {len(todo):,}  "
                          f"total dim failures so far: {total_fail}  {dict(dim_fail_counts)}",
                          flush=True)

    if dim_fail_counts:
        print(f"\nFailure breakdown by dimension: {dict(dim_fail_counts)}", flush=True)
        print("Re-run the same command (no --rerun-all) to retry only the failed "
              "dimensions — successful ones are skipped via the cache.", flush=True)

    # ── Host disambiguation (per-BioSample, multi-value named_hosts only) ─────
    disambig_cache: dict[str, dict] = {}
    if DISAMBIG_CACHE_PATH.exists():
        with open(DISAMBIG_CACHE_PATH) as f:
            for line in f:
                line = line.strip()
                if line:
                    rec = json.loads(line)
                    disambig_cache[rec["bsid"]] = rec

    if args.disambiguate_hosts:
        # Target: BioSamples in a BioProject whose named_hosts has >1 distinct
        # value — single-value BioProjects need no LLM call (host_resolved is
        # computed directly in the row-building loop below).
        disambig_todo = []   # (bp, candidates, bsid)
        for bp, run_list in bp_to_rep.items():
            named_hosts = cache.get(bp, {}).get("llm_named_hosts", []) or []
            candidates = sorted(set(named_hosts))
            if len(candidates) <= 1:
                continue
            for run in run_list:
                bsid = run["BioSample"]
                cached_entry = disambig_cache.get(bsid)
                # Invalidate if named_hosts changed since this BioSample was last
                # resolved (e.g. a hostpath prompt fix re-shapes candidates) — a
                # cached host_resolved value is only meaningful against the exact
                # candidate list it was picked from. Old cache entries (pre-dating
                # this fix) have no "candidates" key at all, so they correctly miss
                # too and get re-resolved once.
                if (cached_entry
                        and cached_entry.get("host_resolved") not in (None, "", "api_error")
                        and cached_entry.get("candidates") == candidates):
                    continue
                disambig_todo.append((bp, candidates, bsid))

        print(f"\nHost disambiguation: {len(disambig_todo):,} BioSamples to resolve "
              f"(multi-value named_hosts BioProjects only)", flush=True)

        if disambig_todo:
            n_disambig_done = 0
            bsid_candidates = {bsid: candidates for bp, candidates, bsid in disambig_todo}
            with open(DISAMBIG_CACHE_PATH, "a") as dfh, \
                 ThreadPoolExecutor(max_workers=args.workers) as pool:
                futures = {
                    pool.submit(_disambiguate_host, bp, candidates, bsid,
                               bs_meta.get(bsid, {}), client): bsid
                    for bp, candidates, bsid in disambig_todo
                }
                for fut in as_completed(futures):
                    bsid = futures[fut]
                    result = fut.result()
                    rec = {"bsid": bsid, "candidates": bsid_candidates[bsid], **result}
                    disambig_cache[bsid] = rec
                    dfh.write(json.dumps(rec) + "\n")
                    dfh.flush()
                    n_disambig_done += 1
                    if n_disambig_done % 50 == 0 or n_disambig_done == len(disambig_todo):
                        print(f"  {n_disambig_done:,} / {len(disambig_todo):,}", flush=True)

    # ── Write samples.tsv ─────────────────────────────────────────────────────
    out_tsv = OUT_DIR / "data" / "samples.tsv"
    rows = []
    for bp, run_list in bp_to_rep.items():
        llm     = cache.get(bp, {})
        bp_info = bp_meta.get(bp, {})

        named_pathogens = llm.get("llm_named_pathogens", []) or []
        named_hosts     = llm.get("llm_named_hosts", []) or []

        # NCBI taxid resolution — one per named item, deterministic lookup (see
        # _util.resolve_taxon_name); "" for anything that doesn't resolve.
        pathogen_taxids = []
        for p in named_pathogens:
            tid, _, _ = resolve_taxon_name(p, name_to_taxid, host_aliases)
            pathogen_taxids.append(str(tid) if tid else "")
        host_taxids = []
        for h in named_hosts:
            tid, _, _ = resolve_taxon_name(h, name_to_taxid, host_aliases)
            host_taxids.append(str(tid) if tid else "")

        llm_cols = {
            "llm_stress":                llm.get("llm_stress",                ""),
            "llm_stress_confidence":     llm.get("llm_stress_confidence",     ""),
            "llm_stress_rationale":      llm.get("llm_stress_rationale",      ""),
            "llm_study_setting":         llm.get("llm_study_setting",         ""),
            "llm_setting_confidence":    llm.get("llm_setting_confidence",    ""),
            "llm_setting_rationale":     llm.get("llm_setting_rationale",     ""),
            "llm_tissue":                llm.get("llm_tissue",                ""),
            "llm_tissue_confidence":     llm.get("llm_tissue_confidence",     ""),
            "llm_tissue_rationale":      llm.get("llm_tissue_rationale",      ""),
            "llm_coinfection_intent":    llm.get("llm_coinfection_intent",    ""),
            "llm_coinfection_confidence": llm.get("llm_coinfection_confidence", ""),
            "llm_coinfection_rationale": llm.get("llm_coinfection_rationale", ""),
            "llm_named_pathogens":       "; ".join(named_pathogens),
            "llm_pathogens_confidence":  llm.get("llm_pathogens_confidence",  ""),
            "llm_pathogens_rationale":   llm.get("llm_pathogens_rationale",   ""),
            "llm_named_hosts":           "; ".join(named_hosts),
            "llm_hosts_confidence":      llm.get("llm_hosts_confidence",      ""),
            "llm_hosts_rationale":       llm.get("llm_hosts_rationale",       ""),
            "llm_symptom_status":        llm.get("llm_symptom_status",        ""),
            "llm_exposure_type":         llm.get("llm_exposure_type",         ""),
            "llm_geographic_location":   llm.get("llm_geographic_location",   ""),
            "llm_library_prep":          llm.get("llm_library_prep",          ""),
            "llm_host_cultivar":         llm.get("llm_host_cultivar",         ""),
            "llm_host_resistance":       llm.get("llm_host_resistance",       ""),
            "llm_named_pathogens_taxids": "; ".join(pathogen_taxids),
            "llm_named_hosts_taxids":     "; ".join(host_taxids),
        }
        lit_cols = {
            "title":           bp_info.get("title",           ""),
            "pmid":            bp_info.get("pmid",            ""),
            "doi":             bp_info.get("doi",             ""),
            "pmid_source":     bp_info.get("pmid_source",     ""),
            "submission_date": bp_info.get("submission_date", ""),
            "pub_date":        bp_info.get("pub_date",        ""),
        }

        distinct_hosts = sorted(set(named_hosts))

        for run in run_list:
            bsid   = run["BioSample"]
            bs_xml = bs_meta.get(bsid, {})

            # host_resolved: deterministic for 0/1 named hosts; looked up from the
            # per-BioSample disambiguation cache for multi-value BioProjects (only
            # populated after --disambiguate-hosts has been run for this BP).
            if len(distinct_hosts) <= 1:
                host_resolved = distinct_hosts[0] if distinct_hosts else ""
                host_resolved_conf = "high" if distinct_hosts else ""
                host_resolved_rationale = "single named host, no disambiguation needed" if distinct_hosts else ""
            else:
                d = disambig_cache.get(bsid, {})
                host_resolved = d.get("host_resolved", "")
                if host_resolved in ("api_error", "unresolved"):
                    host_resolved = "" if host_resolved == "api_error" else "unresolved"
                host_resolved_conf = d.get("confidence", "")
                host_resolved_rationale = d.get("rationale", "")

            host_resolved_taxid = ""
            if host_resolved and host_resolved != "unresolved":
                tid, _, _ = resolve_taxon_name(host_resolved, name_to_taxid, host_aliases)
                host_resolved_taxid = str(tid) if tid else ""

            tiss_agree = _tissue_agreement(llm_cols["llm_tissue"], bs_xml.get("tissue", ""))
            geo_agree  = _geo_agreement(llm_cols["llm_geographic_location"], bs_xml.get("geo_loc_name", ""))
            host_agree = _host_agreement([host_resolved] if host_resolved else named_hosts,
                                        bs_xml.get("host", ""))
            disagree_flag = "disagree" in (tiss_agree, geo_agree, host_agree)
            pathogen_match = _pathogen_match_status(named_pathogens, run.get("stat_pathogens", ""))

            rows.append({
                "BioSample":            bsid,
                "BioProject":           bp,
                "mode":                 run.get("mode",                  ""),
                "biosample_n_runs":     run.get("biosample_n_runs",      ""),
                "stat_host":            run.get("host",                  ""),
                "host_pct":             run.get("host_pct",              ""),
                "stat_pathogens":       run.get("stat_pathogens",        ""),
                "interaction_status":   run.get("interaction_status",    ""),
                "co_infection_flag":    run.get("co_infection_flag",     ""),
                "same_genus_secondary": run.get("same_genus_secondary",  ""),
                "n_pathogens":          run.get("n_pathogens",           ""),
                "fungi_pct":            run.get("fungi_pct",             ""),
                "oomycete_pct":         run.get("oomycete_pct",          ""),
                "nematode_pct":         run.get("nematode_pct",          ""),
                "tissue":               bs_xml.get("tissue",             ""),
                "geo_loc_name":         bs_xml.get("geo_loc_name",       ""),
                "collection_date":      bs_xml.get("collection_date",    ""),
                "dev_stage":            bs_xml.get("dev_stage",          ""),
                "isolation_source":     bs_xml.get("isolation_source",   ""),
                "lat_lon":              bs_xml.get("lat_lon",            ""),
                "bs_host":              bs_xml.get("host",               ""),
                **llm_cols,
                "llm_host_resolved":            host_resolved,
                "llm_host_resolved_confidence": host_resolved_conf,
                "llm_host_resolved_rationale":  host_resolved_rationale,
                "llm_host_resolved_taxid":      host_resolved_taxid,
                "pathogen_match_status":      pathogen_match,
                "tissue_agreement":           tiss_agree,
                "geo_agreement":              geo_agree,
                "host_agreement":             host_agree,
                "metadata_disagreement_flag": disagree_flag,
                **lit_cols,
            })

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_SAMPLE_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # ── Summary ───────────────────────────────────────────────────────────────
    def dist(field: str) -> str:
        c = Counter(r[field] for r in rows if r.get(field) not in (None, ""))
        return "  ".join(f"{k}={v}" for k, v in c.most_common())

    n_disagree = sum(1 for r in rows if r["metadata_disagreement_flag"])
    n_multi_host_rows = sum(1 for r in rows if len(set((r["llm_named_hosts"] or "").split("; ")) - {""}) > 1)
    n_multi_resolved  = sum(1 for r in rows
                            if len(set((r["llm_named_hosts"] or "").split("; ")) - {""}) > 1
                            and r["llm_host_resolved"] not in ("", "unresolved"))
    summary = (
        f"meta_classify summary\n"
        f"BioProjects classified: {len(cache):,} / {len(bp_to_rep):,} (full-text entry gate)\n"
        f"BioSample rows:         {len(rows):,}\n"
        f"Stress:       {dist('llm_stress')}\n"
        f"Setting:      {dist('llm_study_setting')}\n"
        f"Tissue:       {dist('llm_tissue')}\n"
        f"Coinfection intent: {dist('llm_coinfection_intent')}\n"
        f"Symptom:      {dist('llm_symptom_status')}\n"
        f"Exposure:     {dist('llm_exposure_type')}\n"
        f"Library prep: {dist('llm_library_prep')}\n"
        f"Pathogen match status: {dist('pathogen_match_status')}\n"
        f"Metadata disagreement: {n_disagree:,} / {len(rows):,} rows "
        f"({100*n_disagree/max(len(rows),1):.1f}%)\n"
        f"Multi-value named_hosts rows: {n_multi_host_rows:,}  "
        f"resolved (via host_resolved): {n_multi_resolved:,} "
        f"({100*n_multi_resolved/max(n_multi_host_rows,1):.1f}%)"
        + ("" if args.disambiguate_hosts else "  [run --disambiguate-hosts to resolve]") + "\n"
    )
    summary_path = log_dir / "classify_summary.txt"
    summary_path.write_text(summary)
    link_latest(logs_base, summary_path)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}  ({len(rows):,} rows)")
    print(f"Cache:   {CACHE_PATH}  ({len(cache):,} entries)")


if __name__ == "__main__":
    main()
