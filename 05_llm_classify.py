#!/usr/bin/env python3
"""
05_llm_classify.py — LLM-based BioProject classification via Claude API.

Reads biosample_kw.tsv + bioprojects.json + literature.json.
For each BioProject, calls Claude with:
  - BioProject text (title, description, abstract, methods)
  - Keyword classifier output from 04_filter_kw.py (treatment, study_setting,
    named_host) — LLM reviews and flags agreement/disagreement
  - BioSample summary: n, co_infection distribution, XML tissue/geo_loc coverage

Output per BioProject:
  llm_treatment           single|host_study|abiotic_stress|coinf_experiment|
                          surveillance|unclear
  llm_study_setting       field|lab|mixed|unclear
  llm_named_pathogen      primary declared pathogen (from text)
  llm_named_host          primary declared host (from text)
  llm_tissue              inferred tissue type(s) — fallback where XML tissue blank
  llm_kw_treatment_agree  true/false — agrees with keyword treatment
  llm_kw_setting_agree    true/false — agrees with keyword study_setting
  llm_confidence          high|medium|low
  llm_rationale           1-2 sentence explanation

Cache: output/05_llm_classify/data/classify_cache.jsonl (append-only, resumable)
Output: output/05_llm_classify/data/bioproject_llm.tsv

Requires: ANTHROPIC_API_KEY env var

Run from crypt/:
  python 05_llm_classify.py
  python 05_llm_classify.py --workers 16
  python 05_llm_classify.py --rerun-all   # ignore cache, re-classify everything
"""

import argparse
import csv
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from _util import _Tee, link_latest, load_json, make_log_dir

try:
    from openai import OpenAI
except ImportError:
    sys.exit("ERROR: openai package not installed — pip install openai")

# ── Settings ──────────────────────────────────────────────────────────────────

BIOSAMPLE_KW = Path("output/04_filter_kw/data/biosample_kw.tsv")
BP_XML       = Path("output/03a_fetch_xml/data/bioprojects.json")
LITERATURE   = Path("output/03b_fetch_literature/data/literature.json")
OUT_DIR      = Path("output/05_llm_classify")
CACHE_PATH   = OUT_DIR / "data" / "classify_cache.jsonl"

MODEL       = "gpt-4o-mini"
MAX_WORKERS = 8
MAX_RETRIES = 3

# Cache entries missing any of these fields are re-classified
_REQUIRED_FIELDS = {
    "llm_treatment", "llm_study_setting", "llm_named_pathogen", "llm_named_host",
    "llm_tissue", "llm_kw_treatment_agree", "llm_kw_setting_agree",
    "llm_confidence", "llm_rationale",
}

OUTPUT_FIELDS = [
    "BioProject",
    "llm_treatment", "llm_study_setting",
    "llm_named_pathogen", "llm_named_host", "llm_tissue",
    "llm_kw_treatment_agree", "llm_kw_setting_agree",
    "llm_confidence", "llm_rationale",
]

_TREATMENT_VALUES = {
    "single", "host_study", "abiotic_stress",
    "coinf_experiment", "surveillance", "unclear",
}
_SETTING_VALUES   = {"field", "lab", "mixed", "unclear"}
_CONFIDENCE_VALUES = {"high", "medium", "low"}

# ── Prompt ────────────────────────────────────────────────────────────────────

_PROMPT = """\
You are classifying an RNA-seq BioProject from NCBI SRA for a plant cryptic co-infection study.
The study mines public RNA-seq data for evidence of unreported co-infecting organisms detected by
NCBI STAT k-mer taxonomy. Your classification helps distinguish intentional single-pathogen
experiments from field surveys, host studies, and abiotic stress experiments — contexts that
change the interpretation of any detected co-infections.

--- BioProject: {bp} ---
Title: {title}
Description: {description}
Abstract: {abstract}
Methods excerpt: {methods}

--- Keyword classifier output (04_filter_kw.py) ---
treatment:     {kw_treatment}
study_setting: {kw_setting}
named_host:    {kw_named_host}

--- BioSample summary ({n_biosamples} BioSamples) ---
co_infection_flag: single={n_single}, multi_species={n_multi}, multi_kingdom={n_mk}
XML tissue values present:  {tissues}
XML geo_loc_name values:    {geolocnames}

---
TREATMENT (pick one):
  single           host response to a single declared pathogen, or single-pathogen inoculation
  host_study       pure host biology — development, transcriptome assembly, physiology; no pathogen
  abiotic_stress   drought, heat, salt, cold, flooding, UV, nutrient deficiency; no pathogen
  coinf_experiment intentional co-infection or multi-pathogen inoculation experiment
  surveillance     field survey or disease monitoring; no inoculation
  unclear          insufficient information to determine

STUDY_SETTING (pick one):
  field   field-collected samples, natural infection, farm/orchard/commercial crop
  lab     controlled: greenhouse, growth chamber, inoculation, transgenic/mutant lines
  mixed   both field collection and lab manipulation
  unclear cannot determine

Return ONLY valid JSON, no surrounding text:
{{
  "treatment": "<value>",
  "study_setting": "<value>",
  "named_pathogen": "<primary pathogen species declared by the researcher, or empty string>",
  "named_host": "<host plant species declared by the researcher, or empty string>",
  "tissue": "<primary tissue type(s) sampled inferred from text, comma-separated if multiple, or empty string>",
  "kw_treatment_agree": <true|false>,
  "kw_setting_agree": <true|false>,
  "confidence": "<high|medium|low>",
  "rationale": "<1-2 sentences on key decisions>"
}}"""


# ── Helpers ───────────────────────────────────────────────────────────────────

def _bp_summary(rows: list[dict]) -> dict:
    """Aggregate per-BP stats from biosample_kw.tsv rows for the prompt."""
    co   = Counter(r.get("co_infection_flag", "") for r in rows)
    tiss = sorted({r["tissue"] for r in rows if r.get("tissue")})
    geo  = sorted({r["geo_loc_name"].split(":")[0].strip()
                   for r in rows if r.get("geo_loc_name")})
    ref  = rows[0]
    return {
        "n_biosamples":  len(rows),
        "n_single":      co.get("single", 0),
        "n_multi":       co.get("multi_species", 0),
        "n_mk":          co.get("multi_kingdom", 0),
        "tissues":       ", ".join(tiss[:10]) or "none in XML",
        "geolocnames":   ", ".join(geo[:8])   or "none in XML",
        "kw_treatment":  ref.get("treatment", ""),
        "kw_setting":    ref.get("study_setting", ""),
        "kw_named_host": ref.get("named_host", ""),
    }


def _build_prompt(bp: str, summary: dict, meta: dict) -> str:
    return _PROMPT.format(
        bp           = bp,
        title        = (meta.get("title",        "") or "")[:300],
        description  = (meta.get("description",  "") or "")[:500],
        abstract     = (meta.get("abstract",     "") or "")[:2000],
        methods      = (meta.get("methods_text", "") or "")[:3000],
        n_biosamples = summary["n_biosamples"],
        n_single     = summary["n_single"],
        n_multi      = summary["n_multi"],
        n_mk         = summary["n_mk"],
        tissues      = summary["tissues"],
        geolocnames  = summary["geolocnames"],
        kw_treatment = summary["kw_treatment"],
        kw_setting   = summary["kw_setting"],
        kw_named_host= summary["kw_named_host"],
    )


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


def _normalise(raw: dict) -> dict:
    """Coerce raw model output to valid vocabulary; fill missing keys."""
    def _bool_str(v) -> str:
        if isinstance(v, bool):
            return str(v).lower()
        if isinstance(v, str) and v.lower() in ("true", "false"):
            return v.lower()
        return ""

    treatment  = raw.get("treatment", "")
    setting    = raw.get("study_setting", "")
    confidence = raw.get("confidence", "")
    return {
        "llm_treatment":          treatment  if treatment  in _TREATMENT_VALUES  else "unclear",
        "llm_study_setting":      setting    if setting    in _SETTING_VALUES    else "unclear",
        "llm_named_pathogen":     str(raw.get("named_pathogen", "") or ""),
        "llm_named_host":         str(raw.get("named_host",     "") or ""),
        "llm_tissue":             str(raw.get("tissue",         "") or ""),
        "llm_kw_treatment_agree": _bool_str(raw.get("kw_treatment_agree")),
        "llm_kw_setting_agree":   _bool_str(raw.get("kw_setting_agree")),
        "llm_confidence":         confidence if confidence in _CONFIDENCE_VALUES else "low",
        "llm_rationale":          str(raw.get("rationale", "") or ""),
    }


def _classify_one(bp: str, summary: dict, meta: dict,
                  client: "OpenAI") -> dict:
    prompt = _build_prompt(bp, summary, meta)
    for attempt in range(MAX_RETRIES):
        try:
            resp = client.chat.completions.create(
                model      = MODEL,
                max_tokens = 512,
                messages   = [{"role": "user", "content": prompt}],
            )
            raw    = _parse_json(resp.choices[0].message.content)
            result = _normalise(raw)
            result["bp"] = bp
            return result
        except Exception as exc:
            if attempt == MAX_RETRIES - 1:
                print(f"  FAILED {bp}: {exc}", flush=True)
                return {"bp": bp, **{k: "" for k in _REQUIRED_FIELDS}}
            time.sleep(2 ** attempt)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers",   type=int, default=MAX_WORKERS,
                    help=f"parallel API calls (default {MAX_WORKERS})")
    ap.add_argument("--rerun-all", action="store_true",
                    help="re-classify all BioProjects, ignoring cache")
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
    print(f"Model: {MODEL}  workers: {args.workers}", flush=True)

    for path, hint in [
        (BIOSAMPLE_KW, "run 04_filter_kw.py first"),
        (BP_XML,       "run 03a_fetch_xml.py first"),
        (LITERATURE,   "run 03b_fetch_literature.py first"),
    ]:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found — {hint}")

    # ── Load biosample_kw.tsv, group by BioProject ────────────────────────────
    with open(BIOSAMPLE_KW, newline="") as f:
        kw_rows = list(csv.DictReader(f, delimiter="\t"))
    bp_to_rows: dict[str, list] = defaultdict(list)
    for r in kw_rows:
        bp_to_rows[r["BioProject"]].append(r)
    print(f"Loaded {len(kw_rows):,} BioSample rows, "
          f"{len(bp_to_rows):,} unique BioProjects", flush=True)

    # ── Load BP text metadata ─────────────────────────────────────────────────
    bp_xml = load_json(BP_XML)
    lit    = load_json(LITERATURE)
    bp_meta: dict[str, dict] = {}
    for bp in set(bp_xml) | set(lit):
        bp_meta[bp] = {
            "title":        bp_xml.get(bp, {}).get("title",        ""),
            "description":  bp_xml.get(bp, {}).get("description",  ""),
            "abstract":     lit.get(bp,    {}).get("abstract",     ""),
            "methods_text": lit.get(bp,    {}).get("methods_text", ""),
        }

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

    # ── Determine work queue ──────────────────────────────────────────────────
    todo = [bp for bp in bp_to_rows
            if args.rerun_all or not _REQUIRED_FIELDS.issubset(cache.get(bp, {}).keys())]
    n_skip = len(bp_to_rows) - len(todo)
    print(f"To classify: {len(todo):,}  skipping cached: {n_skip:,}", flush=True)

    # ── Classify ──────────────────────────────────────────────────────────────
    done = failed = 0
    if todo:
        with open(CACHE_PATH, "a") as cache_fh, \
             ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(
                    _classify_one,
                    bp,
                    _bp_summary(bp_to_rows[bp]),
                    bp_meta.get(bp, {}),
                    client,
                ): bp
                for bp in todo
            }
            for fut in as_completed(futures):
                bp     = futures[fut]
                result = fut.result()
                cache[bp] = result
                cache_fh.write(json.dumps(result) + "\n")
                cache_fh.flush()
                done += 1
                if not result.get("llm_treatment"):
                    failed += 1
                if done % 100 == 0 or done == len(todo):
                    print(f"  {done:,} / {len(todo):,}  failed: {failed}",
                          flush=True)

    # ── Write output TSV ──────────────────────────────────────────────────────
    out_tsv = OUT_DIR / "data" / "bioproject_llm.tsv"
    rows = []
    for bp in bp_to_rows:
        e = cache.get(bp, {})
        rows.append({
            "BioProject":             bp,
            "llm_treatment":          e.get("llm_treatment",          ""),
            "llm_study_setting":      e.get("llm_study_setting",      ""),
            "llm_named_pathogen":     e.get("llm_named_pathogen",     ""),
            "llm_named_host":         e.get("llm_named_host",         ""),
            "llm_tissue":             e.get("llm_tissue",             ""),
            "llm_kw_treatment_agree": e.get("llm_kw_treatment_agree", ""),
            "llm_kw_setting_agree":   e.get("llm_kw_setting_agree",   ""),
            "llm_confidence":         e.get("llm_confidence",         ""),
            "llm_rationale":          e.get("llm_rationale",          ""),
        })

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_classified = sum(1 for r in rows if r["llm_treatment"])
    treat_counts = Counter(r["llm_treatment"]     for r in rows if r["llm_treatment"])
    set_counts   = Counter(r["llm_study_setting"] for r in rows if r["llm_study_setting"])
    n_ta = sum(1 for r in rows if r["llm_kw_treatment_agree"] == "true")
    n_sa = sum(1 for r in rows if r["llm_kw_setting_agree"]   == "true")
    n    = max(n_classified, 1)

    summary = (
        f"05_llm_classify summary\n"
        f"BioProjects classified:    {n_classified:,} / {len(rows):,}\n"
        f"Treatment:\n"
        + "".join(f"  {k:<22} {v:,}\n" for k, v in treat_counts.most_common())
        + f"Study setting:\n"
        + "".join(f"  {k:<22} {v:,}\n" for k, v in set_counts.most_common())
        + f"Keyword agreement:\n"
        f"  treatment:             {n_ta:,} / {n_classified:,}  ({100*n_ta/n:.1f}%)\n"
        f"  study_setting:         {n_sa:,} / {n_classified:,}  ({100*n_sa/n:.1f}%)\n"
    )
    summary_path = log_dir / "classify_summary.txt"
    summary_path.write_text(summary)
    link_latest(logs_base, summary_path)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}  ({len(rows):,} rows)")
    print(f"Cache:   {CACHE_PATH}  ({len(cache):,} entries)")


if __name__ == "__main__":
    main()
