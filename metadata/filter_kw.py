#!/usr/bin/env python3
"""
04_filter_kw.py — join BioSample metadata with keyword study design classification.

Reads runs.tsv (filtered to biosample_representative rows), biosamples.json,
bioprojects.json, and literature.json to produce biosample_meta.tsv:
one row per BioSample (biosample_representative run).

BioSample is the unit of organisation. BioProject metadata (title, PMID,
study design) is joined as repeated columns.

Two-axis classification (BioProject-level, applied per BioSample):

  treatment (priority order):
    coinf_experiment  intentional co-infection experiment
    abiotic_stress    drought/heat/salt/cold experiment
    host_study        pure host biology (development, assembly, physiology)
    single            single-pathogen or host-response study
    unclear           no text at all

  study_setting (priority order):
    lab               controlled conditions keyword found
    field             field/survey/natural infection keyword found
                      OR geo_loc_name present in BioSample XML
    unclear           text present but no discriminating keywords
    no_data           no text available

named_host: resolved from BioSample XML host attribute first, then
            BioProject title / abstract / description / methods scanning.

LLM classification is a separate step — see llm_classifications.jsonl.

Output:
  output/04_filter_kw/data/biosample_kw.tsv

Usage:
  python 04_filter_kw.py
  python 04_filter_kw.py --hc    # restrict to same_genus_secondary=False
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, link_latest, load_json, make_log_dir
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_vocab import normalise_tissue

# ── Settings ──────────────────────────────────────────────────────────────────

RUNS_TSV         = Path("stat/output/filter_runs/data/runs.tsv")
BIOPROJECTS_PATH = Path("metadata/output/fetch_xml/data/bioprojects.json")
BIOSAMPLES_PATH  = Path("metadata/output/fetch_xml/data/biosamples.json")
LITERATURE_PATH  = Path("metadata/output/fetch_lit/data/literature.json")
LLM_PATH         = Path("metadata/output/llm_classify/data/bioproject_llm.tsv")
PHIBASE_DB       = Path("stat/output/build/data/phibase_db.json")
OUT_DIR          = Path("metadata/output/filter_kw")

OUTPUT_FIELDS = [
    # BioSample identity
    "BioSample", "BioProject", "SRAStudy", "Run", "mode",
    # BioSample XML attributes
    "tissue", "tissue_category", "geo_loc_name", "collection_date", "lat_lon",
    "bs_host", "isolation_source", "dev_stage",
    # Infection summary (from biosample_representative run)
    "n_runs", "co_infection_flag", "same_genus_secondary",
    "library_organism", "library_detected",
    "stat_pathogens", "stat_hosts", "n_pathogens",
    "interaction_status", "host_pct", "analyzed",
    # Named host (resolved canonical name)
    "named_host", "named_host_source",
    # BioProject study design (keyword-based)
    "treatment", "treatment_keywords", "study_setting", "setting_keywords",
    # BioProject literature metadata
    "title", "description", "abstract",
    "primary_pmid", "primary_pub_date", "primary_publication",
    "n_papers_found", "pmid_source", "bp_submission_date",
]

# ── Keyword host extraction ────────────────────────────────────────────────────

_HOST_VERNACULAR: dict[str, str] = {
    "wheat":             "Triticum aestivum",
    "bread wheat":       "Triticum aestivum",
    "spring wheat":      "Triticum aestivum",
    "winter wheat":      "Triticum aestivum",
    "durum wheat":       "Triticum durum",
    "durum":             "Triticum durum",
    "emmer":             "Triticum dicoccoides",
    "einkorn":           "Triticum monococcum",
    "barley":            "Hordeum vulgare",
    "maize":             "Zea mays",
    "corn":              "Zea mays",
    "rice":              "Oryza sativa",
    "sorghum":           "Sorghum bicolor",
    "oat":               "Avena sativa",
    "oats":              "Avena sativa",
    "rye":               "Secale cereale",
    "triticale":         "Triticum aestivum",
    "tomato":            "Solanum lycopersicum",
    "potato":            "Solanum tuberosum",
    "soybean":           "Glycine max",
    "soya":              "Glycine max",
    "arabidopsis":       "Arabidopsis thaliana",
    "tobacco":           "Nicotiana tabacum",
    "rapeseed":          "Brassica napus",
    "canola":            "Brassica napus",
    "oilseed rape":      "Brassica napus",
    "cotton":            "Gossypium hirsutum",
    "grape":             "Vitis vinifera",
    "grapevine":         "Vitis vinifera",
    "sunflower":         "Helianthus annuus",
    "banana":            "Musa acuminata",
    "cassava":           "Manihot esculenta",
    "sugarcane":         "Saccharum officinarum",
    "sugar cane":        "Saccharum officinarum",
    "strawberry":        "Fragaria ananassa",
    "apple":             "Malus domestica",
    "citrus":            "Citrus sinensis",
    "orange":            "Citrus sinensis",
    "peach":             "Prunus persica",
    "lettuce":           "Lactuca sativa",
    "cucumber":          "Cucumis sativus",
    "pepper":            "Capsicum annuum",
    "pea":               "Pisum sativum",
    "bean":              "Phaseolus vulgaris",
    "chickpea":          "Cicer arietinum",
    "lupin":             "Lupinus angustifolius",
    "alfalfa":           "Medicago sativa",
    "poplar":            "Populus trichocarpa",
    "flax":              "Linum usitatissimum",
    "linseed":           "Linum usitatissimum",
    "coffee":            "Coffea arabica",
    "cacao":             "Theobroma cacao",
    "cocoa":             "Theobroma cacao",
    "brachypodium":      "Brachypodium distachyon",
    "setaria":           "Setaria italica",
}


def _build_host_index(db: dict) -> list[tuple[str, str]]:
    """Return (name_lower, canonical_name) sorted longest-first for fast search."""
    host_taxids: set[int] = {int(t) for t in db.get("host_to_seed", {})}
    t2n: dict[int, str]   = {int(k): v for k, v in db.get("taxid_to_name", {}).items()}
    entries: dict[str, str] = {}
    for name, taxid in db.get("name_to_taxid", {}).items():
        if int(taxid) in host_taxids and len(name.split()) >= 2:
            entries[name.lower()] = t2n.get(int(taxid), name)
    for vern, canonical in _HOST_VERNACULAR.items():
        entries[vern.lower()] = canonical
    return sorted(entries.items(), key=lambda x: -len(x[0]))


def _find_meta_host(text: str,
                    index: list[tuple[str, str]]) -> tuple[str, str]:
    """Return (canonical_name, matched_string) for first host found in text."""
    if not text:
        return "", ""
    tl = text.lower()
    for name_lower, canonical in index:
        pos = tl.find(name_lower)
        while pos != -1:
            before = pos == 0 or not tl[pos - 1].isalpha()
            end    = pos + len(name_lower)
            after  = end >= len(tl) or not tl[end].isalpha()
            if before and after:
                return canonical, name_lower
            pos = tl.find(name_lower, pos + 1)
    return "", ""


# ── Study design inference ─────────────────────────────────────────────────────

_COINF_KEYWORDS = {
    "co-infection", "coinfection", "co infection", "dual infection",
    "mixed infection", "dual inoculation", "co-inoculation", "coinoculation",
    "double infection", "co-inoculated", "co-infected",
}
_ABIOTIC_KEYWORDS = {
    "abiotic stress",
    "drought stress", "water deficit", "water stress",
    "heat stress", "high temperature stress", "thermotolerance",
    "salt stress", "salinity stress", "nacl treatment", "osmotic stress",
    "cold stress", "chilling stress", "freezing stress", "frost stress",
    "flooding stress", "waterlogging", "submergence stress",
    "heavy metal stress", "cadmium stress", "zinc toxicity",
    "uv-b stress", "uv stress",
    "nitrogen starvation", "nutrient deficiency", "phosphorus deficiency",
}
_HOST_STUDY_KEYWORDS = {
    "de novo transcriptome", "transcriptome assembly", "reference transcriptome",
    "seed development", "fruit development", "pod development", "kernel development",
    "grain development", "embryo development", "endosperm development",
    "flower development", "floral development", "anther", "pollen tube",
    "root development", "lateral root development",
    "leaf development", "shoot development",
    "developmental stage", "developmental transcriptome",
    "single-cell rna", "scrna-seq", "cell atlas", "tissue atlas",
    "tissue-specific expression", "gene expression atlas", "transcriptome atlas",
    "photosynthesis", "carbon fixation", "circadian", "vernalization", "flowering time",
}
_FIELD_KEYWORDS = {
    "field sample", "field collection", "field isolate", "field strain",
    "field survey", "field study", "field trial", "field condition",
    "field-collected", "field-grown", "field-infected", "field experiment",
    "field rna-seq", "field transcriptomic", "field populations", "field-sampled",
    "disease survey", "virus survey", "pathogen survey", "population survey",
    "disease surveillance", "pathogen surveillance",
    "surveillance", "epidemiology", "epidemiological",
    "natural infection", "naturally infected", "naturally occurring",
    "wild population", "wild plant", "wild-growing", "wild accession",
    "farm", "orchard", "commercial crop", "commercial orchard", "commercial farm",
    "pathogenomics", "outbreak",
    "growers field", "farmers field", "production field",
    "commercial production", "agronomic field",
}
_LAB_KEYWORDS = {
    "wild type", "wild-type",
    "mutant", "knockout", "knock-out",
    "transgenic", "overexpression", "overexpressing", "overexpressor",
    "rnai", "rna interference", "gene silencing", "gene knockdown",
    "in vitro", "crispr", "t-dna",
    "inoculated", "inoculation",
    "artificial inoculation", "spray inoculation", "wound inoculation",
    "detached leaf", "leaf disc",
    "days post inoculation", "hours post inoculation",
    "growth chamber", "greenhouse", "glasshouse", "growth room",
    "controlled condition", "controlled environment",
    "axenic", "hydroponic", "potted plant", "pot-grown",
}


def _classify(title: str, description: str,
              pub_text: str = "",
              methods_text: str = "",
              has_geoloc: bool = False) -> tuple[str, str, list[str], list[str]]:
    """Return (treatment, study_setting, treat_keywords, set_keywords).

    has_geoloc: if True and no lab/field keywords found, setting → "field".
    """
    no_text = not (title.strip() or description.strip()
                   or pub_text.strip() or methods_text.strip())
    if no_text:
        return "unclear", "no_data", [], []

    text = (title + " " + description + " " + pub_text + " " + methods_text).lower()

    coinf_hits   = [kw for kw in _COINF_KEYWORDS      if kw in text]
    abiotic_hits = [kw for kw in _ABIOTIC_KEYWORDS    if kw in text]
    host_hits    = [kw for kw in _HOST_STUDY_KEYWORDS if kw in text]
    lab_hits     = [kw for kw in _LAB_KEYWORDS        if kw in text]
    field_hits   = [kw for kw in _FIELD_KEYWORDS      if kw in text]

    if coinf_hits:
        treatment  = "coinf_experiment"
        treat_kws  = coinf_hits
    elif abiotic_hits:
        treatment  = "abiotic_stress"
        treat_kws  = abiotic_hits
    elif host_hits:
        treatment  = "host_study"
        treat_kws  = host_hits
    else:
        treatment  = "single"
        treat_kws  = []

    if lab_hits:
        study_setting = "lab"
        set_kws = lab_hits
    elif field_hits:
        study_setting = "field"
        set_kws = field_hits
    elif has_geoloc:
        study_setting = "field"
        set_kws = ["geo_loc_name"]
    else:
        study_setting = "unclear"
        set_kws = []

    return treatment, study_setting, treat_kws, set_kws


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hc", action="store_true",
                    help="restrict to high-confidence runs (same_genus_secondary=False)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    sys.stdout = _Tee(log_dir / "filter.log")
    link_latest(logs_base, log_dir / "filter.log")

    out_tsv = OUT_DIR / "data" / "biosample_kw.tsv"

    for path, hint in [
        (RUNS_TSV,         "run 02_filter_runs.py first"),
        (BIOPROJECTS_PATH, "run 03a_fetch_xml.py first"),
        (BIOSAMPLES_PATH,  "run 03a_fetch_xml.py first"),
        (LITERATURE_PATH,  "run 03b_fetch_literature.py first"),
    ]:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found — {hint}")

    # ── Load runs (biosample_representative only) ─────────────────────────────
    with open(RUNS_TSV, newline="") as f:
        all_rows = list(csv.DictReader(f, delimiter="\t"))
    rep_rows = [r for r in all_rows if r.get("biosample_representative") == "True"]
    if args.hc:
        rep_rows = [r for r in rep_rows if r.get("same_genus_secondary") == "False"]
    print(f"Loaded {len(all_rows):,} runs from {RUNS_TSV}", flush=True)
    print(f"BioSample-representative rows: {len(rep_rows):,}"
          + (" (hc filter)" if args.hc else ""), flush=True)

    # ── Load metadata sources ─────────────────────────────────────────────────
    bp_xml     = load_json(BIOPROJECTS_PATH)
    biosamples = load_json(BIOSAMPLES_PATH)
    lit        = load_json(LITERATURE_PATH)
    print(f"BioProjects: {len(bp_xml):,}  BioSamples: {len(biosamples):,}  "
          f"Literature: {len(lit):,}", flush=True)

    # ── Load LLM tissue ───────────────────────────────────────────────────────
    llm_tissue: dict[str, str] = {}
    if LLM_PATH.exists():
        with open(LLM_PATH, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                llm_tissue[row["BioProject"]] = row.get("llm_tissue", "")
        print(f"LLM tissue: {sum(1 for v in llm_tissue.values() if v):,} BioProjects with value",
              flush=True)

    bp_cache: dict[str, dict] = {}
    for bp in set(bp_xml) | set(lit):
        xml_data = bp_xml.get(bp, {})
        lit_data = lit.get(bp, {})
        entry: dict = {
            "title":               xml_data.get("title", ""),
            "description":         xml_data.get("description", ""),
            "bp_submission_date":  xml_data.get("submission_date", ""),
            **{k: lit_data.get(k, "") for k in [
                "primary_pmid", "primary_pub_date", "primary_publication",
                "abstract", "methods_text", "pmcid", "n_papers_found", "pmid_source",
            ]},
        }
        entry["primary_pub_text"] = (entry["title"] + " " + entry["abstract"]).strip()
        bp_cache[bp] = entry
    print(f"Merged BP metadata: {len(bp_cache):,} entries", flush=True)

    # ── Host name index ───────────────────────────────────────────────────────
    host_index: list[tuple[str, str]] = []
    if PHIBASE_DB.exists():
        with open(PHIBASE_DB) as f:
            db = json.load(f)
        host_index = _build_host_index(db)
        print(f"Host name index: {len(host_index):,} entries", flush=True)

    label = "high-confidence (same_genus_secondary=False)" if args.hc else "all co-infected"
    print(f"Co-infection filter: {label}", flush=True)

    # ── Build output rows ─────────────────────────────────────────────────────
    results = []
    for run in rep_rows:
        bs_acc = run.get("BioSample", "").strip()
        bp     = run.get("BioProject", "").strip()
        bs     = biosamples.get(bs_acc, {})
        meta   = bp_cache.get(bp, {})

        named_host, named_host_source = "", ""
        if host_index:
            raw_bs_host = bs.get("host", "")
            if raw_bs_host:
                h, _ = _find_meta_host(raw_bs_host, host_index)
                if h:
                    named_host        = h
                    named_host_source = "biosample_host"
            if not named_host:
                for src_name, src_text in (
                    ("title",    meta.get("title", "")),
                    ("abstract", meta.get("abstract", "")),
                    ("desc",     meta.get("description", "")),
                    ("methods",  meta.get("methods_text", "")),
                ):
                    h, _ = _find_meta_host(src_text, host_index)
                    if h:
                        named_host        = h
                        named_host_source = src_name
                        break

        has_geoloc = bool(bs.get("geo_loc_name", ""))
        treatment, setting, treat_kws, set_kws = _classify(
            meta.get("title", ""),
            meta.get("description", ""),
            meta.get("primary_pub_text", ""),
            meta.get("methods_text", ""),
            has_geoloc=has_geoloc,
        )

        results.append({
            "BioSample":           bs_acc,
            "BioProject":          bp,
            "SRAStudy":            run.get("SRAStudy", ""),
            "Run":                 run.get("Run", ""),
            "mode":                run.get("mode", ""),
            "tissue":              bs.get("tissue", ""),
            "tissue_category":     (
                # cascade: biosample tissue → isolation_source → llm_tissue
                next(
                    (cat for raw in [
                        bs.get("tissue", ""),
                        bs.get("isolation_source", ""),
                        llm_tissue.get(bp, ""),
                    ] if (cat := normalise_tissue(raw)) != "unknown"),
                    "unknown"
                )
            ),
            "geo_loc_name":        bs.get("geo_loc_name", ""),
            "collection_date":     bs.get("collection_date", ""),
            "lat_lon":             bs.get("lat_lon", ""),
            "bs_host":             bs.get("host", ""),
            "isolation_source":    bs.get("isolation_source", ""),
            "dev_stage":           bs.get("dev_stage", ""),
            "n_runs":              run.get("biosample_n_runs", ""),
            "co_infection_flag":   run.get("co_infection_flag", ""),
            "same_genus_secondary": run.get("same_genus_secondary", ""),
            "library_organism":    run.get("library_organism", ""),
            "library_detected":    run.get("library_detected", ""),
            "stat_pathogens":      run.get("stat_pathogens", ""),
            "stat_hosts":          run.get("stat_hosts", ""),
            "n_pathogens":         run.get("n_pathogens", ""),
            "interaction_status":  run.get("interaction_status", ""),
            "host_pct":            run.get("host_pct", ""),
            "analyzed":            run.get("analyzed", ""),
            "named_host":          named_host,
            "named_host_source":   named_host_source,
            "treatment":           treatment,
            "treatment_keywords":  "|".join(treat_kws),
            "study_setting":       setting,
            "setting_keywords":    "|".join(set_kws),
            "title":               meta.get("title", ""),
            "description":         meta.get("description", ""),
            "abstract":            meta.get("abstract", ""),
            "primary_pmid":        meta.get("primary_pmid", ""),
            "primary_pub_date":    meta.get("primary_pub_date", ""),
            "primary_publication": meta.get("primary_publication", ""),
            "n_papers_found":      meta.get("n_papers_found", 0),
            "pmid_source":         meta.get("pmid_source", "none"),
            "bp_submission_date":  meta.get("bp_submission_date", ""),
        })

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    # ── Summary ───────────────────────────────────────────────────────────────
    from collections import defaultdict

    n_total      = len(results)
    n_coinf      = sum(1 for r in results if r["co_infection_flag"] != "single")
    n_hc         = sum(1 for r in results if r["co_infection_flag"] != "single"
                       and r["same_genus_secondary"] == "False")
    n_geolocated  = sum(1 for r in results if r["geo_loc_name"])
    n_tissue      = sum(1 for r in results if r["tissue"])
    n_aerial      = sum(1 for r in results if r["tissue_category"] in ("leaf", "aerial_other"))
    n_non_aerial  = sum(1 for r in results if r["tissue_category"] in ("root", "seed_fruit", "whole_plant"))
    n_named_host = sum(1 for r in results if r["named_host"])
    n_bs_host    = sum(1 for r in results if r["named_host_source"] == "biosample_host")
    n_with_pmid  = sum(1 for r in results if r["primary_pmid"])
    n_with_title = sum(1 for r in results if r["title"])
    n_with_methods = sum(1 for r in results
                         if bp_cache.get(r["BioProject"], {}).get("methods_text"))

    treat_counts: dict[str, int] = defaultdict(int)
    set_counts:   dict[str, int] = defaultdict(int)
    for r in results:
        treat_counts[r["treatment"]]   += 1
        set_counts[r["study_setting"]] += 1

    summary = (
        f"04_filter_kw summary\n"
        f"Co-infection filter: {label}\n"
        f"BioSamples total:         {n_total:,}\n"
        f"  co-infected:            {n_coinf:,}  ({100*n_coinf/max(n_total,1):.1f}%)\n"
        f"  high-confidence:        {n_hc:,}  ({100*n_hc/max(n_total,1):.1f}%)\n"
        f"BioSample coverage:\n"
        f"  geo_loc_name:           {n_geolocated:,}  ({100*n_geolocated/max(n_total,1):.1f}%)\n"
        f"  tissue:                 {n_tissue:,}  ({100*n_tissue/max(n_total,1):.1f}%)\n"
        f"    aerial:               {n_aerial:,}\n"
        f"    non-aerial:           {n_non_aerial:,}\n"
        f"Named host:\n"
        f"  resolved:               {n_named_host:,}  ({100*n_named_host/max(n_total,1):.1f}%)\n"
        f"  from biosample_host:    {n_bs_host:,}  ({100*n_bs_host/max(n_total,1):.1f}%)\n"
        f"Literature:\n"
        f"  with title:             {n_with_title:,}\n"
        f"  with PMID:              {n_with_pmid:,}\n"
        f"  with PMC methods text:  {n_with_methods:,}\n"
        f"Treatment (auto):\n"
        f"  coinf_experiment:       {treat_counts['coinf_experiment']:,}\n"
        f"  abiotic_stress:         {treat_counts['abiotic_stress']:,}\n"
        f"  host_study:             {treat_counts['host_study']:,}\n"
        f"  single:                 {treat_counts['single']:,}\n"
        f"  unclear (no text):      {treat_counts['unclear']:,}\n"
        f"Study setting (auto):\n"
        f"  field:                  {set_counts['field']:,}\n"
        f"  lab:                    {set_counts['lab']:,}\n"
        f"  unclear:                {set_counts['unclear']:,}\n"
        f"  no_data:                {set_counts['no_data']:,}\n"
    )
    summary_path = log_dir / "filter_summary.txt"
    summary_path.write_text(summary)
    link_latest(logs_base, summary_path)
    print(f"\n{summary}")
    print(f"Written: {out_tsv}  ({n_total:,} rows)")


if __name__ == "__main__":
    main()
