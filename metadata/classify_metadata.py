#!/usr/bin/env python3
"""
metadata/classify_metadata.py — LLM BioProject classification + BioSample output table.

Reads bioprojects.json + biosamples.json + runs.tsv directly.
Makes THREE separate LLM calls per BioProject — one per classification dimension — so
each gets its own focused prompt, rules, confidence score, and rationale:

  stress   (biotic|abiotic|none|unclear)   + llm_stress_confidence/rationale
  setting  (field|lab|unclear)             + llm_setting_confidence/rationale
  tissue   (aerial|non-aerial|unclear)     + llm_tissue_confidence/rationale

Output: metadata/output/classify_metadata/data/samples.tsv
        one row per biosample_representative BioSample

Cache: metadata/output/classify_metadata/data/classify_cache.jsonl
       append-only JSONL, one JSON line per BioProject

Requires: OPENAI_API_KEY env var

Run from crypt/:
  python metadata/classify_metadata.py
  python metadata/classify_metadata.py --workers 8
  python metadata/classify_metadata.py --rerun-all
  python metadata/classify_metadata.py --focus setting   # rerun one dimension only
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
from _util import _Tee, link_latest, load_json, make_log_dir

try:
    from openai import OpenAI
except ImportError:
    sys.exit("ERROR: openai package not installed — pip install openai")

# ── Settings ──────────────────────────────────────────────────────────────────

BIOPROJECTS = Path("metadata/output/meta_search/data/bioprojects.json")
BIOSAMPLES  = Path("metadata/output/meta_search/data/biosamples.json")
RUNS_TSV    = Path("stat/output/stat_filter/data/runs.tsv")
OUT_DIR     = Path("metadata/output/classify_metadata")
CACHE_PATH  = OUT_DIR / "data" / "classify_cache.jsonl"

MODEL       = "gpt-4o-mini"
MAX_WORKERS = 8
MAX_RETRIES = 6
# 3 calls per BP. gpt-4o-mini allows ~2k RPM on tier 1; cap at 150 req/min to be safe.
_RATE_LOCK    = threading.Lock()
_RATE_MIN_GAP = 60.0 / 150       # ~0.4 s between calls
_last_call    = 0.0

_FOCUS_VALUES = {"all", "stress", "setting", "tissue"}

_STRESS_VALUES  = {"biotic", "abiotic", "none", "unclear"}
_SETTING_VALUES = {"field", "lab", "unclear"}
_TISSUE_VALUES  = {"aerial", "non-aerial", "unclear"}
_CONF_VALUES    = {"high", "medium", "low"}

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
    "llm_named_pathogen", "llm_named_host",
    # LLM — setting
    "llm_study_setting", "llm_setting_confidence", "llm_setting_rationale",
    # LLM — tissue
    "llm_tissue", "llm_tissue_confidence", "llm_tissue_rationale",
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

_CTX = """\
--- BioProject: {bp} ---
Library mode: {mode_desc}
Title: {title}
Description: {description}
Abstract: {abstract}
Methods excerpt: {methods}
"""

_PROMPT_STRESS = """\
You are classifying the STRESS TYPE of an RNA-seq BioProject from plant pathology research.

{ctx}
STAT-detected pathogens in sequencing data: {stat_pathogen_summary}

---
What stress treatment were the plants subjected to? Pick ONE:

biotic   — exposure to any pathogen, pest, or parasite; deliberate inoculation or
           naturally occurring infection; disease resistance or susceptibility trials;
           any study where a named pathogen is the biological focus; co-infection designs
abiotic  — physical or chemical stress ONLY: drought, heat, cold, salt, flooding, UV,
           wounding, heavy metals, nutrient deficiency; NO pathogen involved
none     — no stress: pure host biology, developmental atlas, genotype comparison,
           QTL/GWAS transcriptomics, phenology; neither pathogen nor abiotic treatment
unclear  — insufficient information to determine

RULES:
1. MAL mode: library organism IS the pathogen → almost certainly biotic
2. Named pathogen in title/abstract → biotic
3. Hormone treatments (jasmonate, salicylate, ethylene) WITHOUT a pathogen → abiotic
4. "Disease resistance", "biotic stress", "infection" → biotic
5. Return unclear ONLY if no pathogen, disease, or stress language appears at all

Return ONLY valid JSON, no surrounding text:
{{"stress": "<biotic|abiotic|none|unclear>",
  "named_pathogen": "<primary pathogen species named by researcher, or empty string>",
  "named_host": "<host plant species named by researcher, or empty string>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences citing specific text above>"}}"""


_PROMPT_SETTING = """\
You are classifying the STUDY SETTING (field vs laboratory) of an RNA-seq BioProject
from plant pathology research.

{ctx}
--- BioSample collection signals ({n_biosamples} BioSamples) ---
XML geo_loc_name (countries):        {geolocnames}
XML collection_date present:         {n_with_collection_date} of {n_biosamples} BioSamples
BioSample descriptions:              {bs_descriptions}

---
Where were the sequenced plant samples collected or grown? Pick ONE:

field   — plants in NATURAL or AGRICULTURAL conditions: farm, orchard, commercial crop,
          natural population, field survey; tissue obtained from field-growing plants.
          Includes field pathogenomics / disease surveillance / portable diagnostics where
          infected plant tissue is collected directly in the field without a lab inoculation
          step (e.g. samples stored in RNAlater in the field, MARPLE diagnostics,
          global disease monitoring collecting leaves across farms or countries).
lab     — controlled conditions: greenhouse, growth chamber, pot experiment, axenic
          culture, detached leaf assay, in vitro.
          Classify lab ONLY when text EXPLICITLY states greenhouse, growth chamber,
          controlled inoculation, or similar. Do NOT default to lab when absent.
unclear — cannot determine — USE THIS as the DEFAULT when no explicit field or lab
          language is present in the provided text.

RULES:
1. BioSample descriptions containing "field sample", "field-collected", "collected from
   [farm/location]" → field
2. "Field isolate" or "field-derived isolate" describing only the PATHOGEN STRAIN that
   was then used in lab inoculation → lab (not field)
3. geo_loc_name country alone does NOT indicate field; lab studies record institution
   country. Multiple distinct countries WITH collection dates → stronger field signal.
4. "Collected from", "collected by [person]", "stored in RNAlater", "portable sequencer",
   "point-of-care", "disease surveillance", "MARPLE" → field
5. Arabidopsis thaliana, Nicotiana benthamiana, Brachypodium distachyon → almost always
   lab. Classify field ONLY if text explicitly describes field collection of plant material.
6. Do NOT write "inoculation" or "controlled conditions" in your rationale unless those
   exact words appear in the text above. When evidence is sparse, return unclear.

Return ONLY valid JSON, no surrounding text:
{{"study_setting": "<field|lab|unclear>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences citing ONLY evidence present in the text or signals above>"}}"""


_PROMPT_TISSUE = """\
You are classifying the TISSUE TYPE sequenced in an RNA-seq BioProject from plant
pathology research.

{ctx}
--- BioSample tissue signals ({n_biosamples} BioSamples) ---
XML tissue values:         {tissues}
XML dev_stage values:      {dev_stages}
BioSample descriptions:    {bs_descriptions}

---
What plant tissue or organ was the primary source of RNA? Pick ONE:

aerial      — above-ground tissue: leaf, leaflet, blade, stem, shoot, hypocotyl,
              flower, sepal, petal, bud, inflorescence, spike, panicle, rachis,
              seed, grain, endosperm, embryo, fruit, pericarp, pollen, anther,
              petiole, canopy, pod, flag leaf, urediniospore (rust spores)
non-aerial  — below-ground or storage tissue: root, lateral root, root tip, tuber,
              bulb, corm, rhizome, crown, stolon
unclear     — whole plant / seedling without tissue specification; mixed shoot+root;
              callus, protoplast, suspension culture; not stated; cannot determine

RULES:
1. XML tissue values are the most reliable signal — if all/most agree, use that
2. "Leaf rust", "stripe rust on leaf", "flag leaf", "infected leaf" → aerial
3. "Root rot", "crown rot", "Phytophthora root infection" → non-aerial
4. MAL mode: tissue is where the pathogen resides — rust/powdery mildew → aerial;
   Phytophthora/Pythium root/crown rot → non-aerial; check title if unclear
5. Seeds and grain are aerial (harvested above ground)
6. If methods describe a specific tissue dissection, trust that over XML tissue values

Return ONLY valid JSON, no surrounding text:
{{"tissue": "<aerial|non-aerial|unclear>",
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences>"}}"""


# ── BioSample summary ──────────────────────────────────────────────────────────

def _bp_summary(rows: list[dict], bs_meta: dict) -> dict:
    co      = Counter(r.get("co_infection_flag", "") for r in rows)
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
        "n_single":               co.get("single", 0),
        "n_multi":                co.get("multi_species", 0),
        "n_mk":                   co.get("multi_kingdom", 0),
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
        bp        = bp,
        mode_desc = _MODE_DESC.get(summary["mode"], ""),
        title       = (meta.get("title",        "") or "")[:300],
        description = (meta.get("description",  "") or "")[:500],
        abstract    = (meta.get("abstract",     "") or "")[:2000],
        methods     = (meta.get("methods_text", "") or "")[:3000],
    )


def _stress_prompt(bp: str, summary: dict, meta: dict) -> str:
    return _PROMPT_STRESS.format(
        ctx                  = _ctx(bp, summary, meta),
        stat_pathogen_summary = summary["stat_pathogen_summary"],
    )


def _setting_prompt(bp: str, summary: dict, meta: dict) -> str:
    return _PROMPT_SETTING.format(
        ctx                   = _ctx(bp, summary, meta),
        n_biosamples          = summary["n_biosamples"],
        geolocnames           = summary["geolocnames"],
        n_with_collection_date = summary["n_with_collection_date"],
        bs_descriptions       = summary["bs_descriptions"],
    )


def _tissue_prompt(bp: str, summary: dict, meta: dict) -> str:
    return _PROMPT_TISSUE.format(
        ctx             = _ctx(bp, summary, meta),
        n_biosamples    = summary["n_biosamples"],
        tissues         = summary["tissues"],
        dev_stages      = summary["dev_stages"],
        bs_descriptions = summary["bs_descriptions"],
    )


# ── API call helpers ───────────────────────────────────────────────────────────

def _rate_wait() -> None:
    global _last_call
    with _RATE_LOCK:
        gap = time.time() - _last_call
        if gap < _RATE_MIN_GAP:
            time.sleep(_RATE_MIN_GAP - gap)
        _last_call = time.time()


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


def _call_dim(prompt: str, client: "OpenAI", bp: str, dim: str) -> dict:
    for attempt in range(MAX_RETRIES):
        _rate_wait()
        try:
            resp = client.chat.completions.create(
                model      = MODEL,
                max_tokens = 300,
                messages   = [{"role": "user", "content": prompt}],
            )
            return _parse_json(resp.choices[0].message.content)
        except Exception as exc:
            if attempt == MAX_RETRIES - 1:
                print(f"  FAILED {bp} [{dim}]: {exc}", flush=True)
                return {}
            time.sleep(min(2 ** attempt, 60))
    return {}


def _classify_bp(bp: str, summary: dict, meta: dict,
                 client: "OpenAI", focus: str) -> dict:
    result: dict = {"bp": bp}
    conf  = _CONF_VALUES

    if focus in ("all", "stress"):
        raw = _call_dim(_stress_prompt(bp, summary, meta), client, bp, "stress")
        stress = raw.get("stress", "")
        result["llm_stress"]           = stress     if stress     in _STRESS_VALUES else "unclear"
        result["llm_named_pathogen"]   = str(raw.get("named_pathogen", "") or "")
        result["llm_named_host"]       = str(raw.get("named_host",     "") or "")
        result["llm_stress_confidence"] = raw.get("confidence", "") if raw.get("confidence") in conf else "low"
        result["llm_stress_rationale"]  = str(raw.get("rationale", "") or "")

    if focus in ("all", "setting"):
        raw = _call_dim(_setting_prompt(bp, summary, meta), client, bp, "setting")
        setting = raw.get("study_setting", "")
        result["llm_study_setting"]      = setting if setting in _SETTING_VALUES else "unclear"
        result["llm_setting_confidence"] = raw.get("confidence", "") if raw.get("confidence") in conf else "low"
        result["llm_setting_rationale"]  = str(raw.get("rationale", "") or "")

    if focus in ("all", "tissue"):
        raw = _call_dim(_tissue_prompt(bp, summary, meta), client, bp, "tissue")
        tissue = raw.get("tissue", "")
        result["llm_tissue"]             = tissue if tissue in _TISSUE_VALUES else "unclear"
        result["llm_tissue_confidence"]  = raw.get("confidence", "") if raw.get("confidence") in conf else "low"
        result["llm_tissue_rationale"]   = str(raw.get("rationale", "") or "")

    return result


# ── Completeness check ─────────────────────────────────────────────────────────

_ALL_REQUIRED = {
    "llm_stress", "llm_stress_confidence",
    "llm_study_setting", "llm_setting_confidence",
    "llm_tissue", "llm_tissue_confidence",
    "llm_named_pathogen", "llm_named_host",
}
_DIM_REQUIRED = {
    "stress":  {"llm_stress",        "llm_stress_confidence"},
    "setting": {"llm_study_setting", "llm_setting_confidence"},
    "tissue":  {"llm_tissue",        "llm_tissue_confidence"},
}


def _complete(entry: dict, focus: str = "all") -> bool:
    fields = _ALL_REQUIRED if focus == "all" else _DIM_REQUIRED[focus]
    if not fields.issubset(entry.keys()):
        return False
    # value must not be empty for the primary classification field
    val_key = {"all": "llm_stress", "stress": "llm_stress",
                "setting": "llm_study_setting", "tissue": "llm_tissue"}.get(focus)
    return bool(entry.get(val_key))


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
    ]:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found — {hint}")

    # ── Load inputs ───────────────────────────────────────────────────────────
    bp_meta: dict[str, dict] = load_json(BIOPROJECTS)
    bs_meta: dict[str, dict] = load_json(BIOSAMPLES)

    with open(RUNS_TSV, newline="") as f:
        all_runs = list(csv.DictReader(f, delimiter="\t"))

    rep_rows = [r for r in all_runs if r.get("biosample_representative") == "True"]
    bp_to_rep: dict[str, list[dict]] = defaultdict(list)
    for r in rep_rows:
        bp_to_rep[r["BioProject"]].append(r)

    print(f"BioProjects: {len(bp_to_rep):,}  "
          f"representative BioSamples: {len(rep_rows):,}", flush=True)

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

    done = failed = 0
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
                ): bp
                for bp in todo
            }
            for fut in as_completed(futures):
                bp     = futures[fut]
                result = fut.result()
                # merge focused result into existing cache entry
                merged = dict(cache.get(bp, {"bp": bp}))
                merged.update(result)
                cache[bp] = merged
                cache_fh.write(json.dumps(merged) + "\n")
                cache_fh.flush()
                done += 1
                if not result.get("llm_stress") and args.focus in ("all", "stress"):
                    failed += 1
                if done % 100 == 0 or done == len(todo):
                    print(f"  {done:,} / {len(todo):,}  failed: {failed}", flush=True)

    # ── Write samples.tsv ─────────────────────────────────────────────────────
    out_tsv = OUT_DIR / "data" / "samples.tsv"
    rows = []
    for bp, run_list in bp_to_rep.items():
        llm     = cache.get(bp, {})
        bp_info = bp_meta.get(bp, {})
        llm_cols = {
            "llm_stress":              llm.get("llm_stress",              ""),
            "llm_stress_confidence":   llm.get("llm_stress_confidence",   ""),
            "llm_stress_rationale":    llm.get("llm_stress_rationale",    ""),
            "llm_named_pathogen":      llm.get("llm_named_pathogen",      ""),
            "llm_named_host":          llm.get("llm_named_host",          ""),
            "llm_study_setting":       llm.get("llm_study_setting",       ""),
            "llm_setting_confidence":  llm.get("llm_setting_confidence",  ""),
            "llm_setting_rationale":   llm.get("llm_setting_rationale",   ""),
            "llm_tissue":              llm.get("llm_tissue",              ""),
            "llm_tissue_confidence":   llm.get("llm_tissue_confidence",   ""),
            "llm_tissue_rationale":    llm.get("llm_tissue_rationale",    ""),
        }
        lit_cols = {
            "title":           bp_info.get("title",           ""),
            "pmid":            bp_info.get("pmid",            ""),
            "doi":             bp_info.get("doi",             ""),
            "pmid_source":     bp_info.get("pmid_source",     ""),
            "submission_date": bp_info.get("submission_date", ""),
            "pub_date":        bp_info.get("pub_date",        ""),
        }
        for run in run_list:
            bsid   = run["BioSample"]
            bs_xml = bs_meta.get(bsid, {})
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
                **lit_cols,
            })

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_SAMPLE_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # ── Summary ───────────────────────────────────────────────────────────────
    stress_counts  = Counter(r["llm_stress"]        for r in rows if r["llm_stress"])
    set_counts     = Counter(r["llm_study_setting"] for r in rows if r["llm_study_setting"])
    tiss_counts    = Counter(r["llm_tissue"]        for r in rows if r["llm_tissue"])

    # confidence breakdown per dimension
    def conf_dist(field: str) -> str:
        c = Counter(r[field] for r in rows if r.get(field))
        return "  ".join(f"{k}={v}" for k, v in c.most_common())

    summary = (
        f"classify_metadata summary\n"
        f"BioProjects classified: {len(cache):,} / {len(bp_to_rep):,}\n"
        f"BioSample rows:         {len(rows):,}\n"
        f"Stress:\n"
        + "".join(f"  {k:<22} {v:,}\n" for k, v in stress_counts.most_common())
        + f"  confidence: {conf_dist('llm_stress_confidence')}\n"
        + "Study setting:\n"
        + "".join(f"  {k:<22} {v:,}\n" for k, v in set_counts.most_common())
        + f"  confidence: {conf_dist('llm_setting_confidence')}\n"
        + "Tissue:\n"
        + "".join(f"  {k:<22} {v:,}\n" for k, v in tiss_counts.most_common())
        + f"  confidence: {conf_dist('llm_tissue_confidence')}\n"
    )
    summary_path = log_dir / "classify_summary.txt"
    summary_path.write_text(summary)
    link_latest(logs_base, summary_path)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}  ({len(rows):,} rows)")
    print(f"Cache:   {CACHE_PATH}  ({len(cache):,} entries)")


if __name__ == "__main__":
    main()
