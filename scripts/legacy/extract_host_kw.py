#!/usr/bin/env python3
"""
scripts/extract_host_kw.py — keyword-based host extraction from BioProject metadata.

Searches title + description + abstract + methods_text for plant host names using:
  1. PHI-base host scientific names (from phibase_db.json name_to_taxid, filtered to hosts)
  2. Vernacular crop names (wheat, barley, maize, …)

Outputs:
  output/analysis/host_kw.tsv         one row per BioProject (all 1,750)
  Console: coverage + agreement stats vs STAT host column

Run from crypt/:
  python scripts/extract_host_kw.py
"""

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

CRYPT_TSV  = Path("output/02_filter_runs/data/crypt.tsv")
PHIBASE_DB = Path("output/00_build/data/phibase_db.json")
META_TSV   = Path("output/04_filter_meta/data/bioproject_meta.tsv")
OUT_DIR    = Path("output/analysis")
OUT_TSV    = OUT_DIR / "host_kw.tsv"

# Common crop/plant vernacular names → canonical scientific name
VERNACULAR: dict[str, str] = {
    "wheat":              "Triticum aestivum",
    "bread wheat":        "Triticum aestivum",
    "spring wheat":       "Triticum aestivum",
    "winter wheat":       "Triticum aestivum",
    "durum wheat":        "Triticum durum",
    "durum":              "Triticum durum",
    "emmer":              "Triticum dicoccoides",
    "einkorn":            "Triticum monococcum",
    "barley":             "Hordeum vulgare",
    "maize":              "Zea mays",
    "corn":               "Zea mays",
    "rice":               "Oryza sativa",
    "sorghum":            "Sorghum bicolor",
    "oat":                "Avena sativa",
    "oats":               "Avena sativa",
    "rye":                "Secale cereale",
    "triticale":          "Triticum aestivum",   # hybrid wheat/rye
    "tomato":             "Solanum lycopersicum",
    "potato":             "Solanum tuberosum",
    "soybean":            "Glycine max",
    "soya":               "Glycine max",
    "arabidopsis":        "Arabidopsis thaliana",
    "tobacco":            "Nicotiana tabacum",
    "rapeseed":           "Brassica napus",
    "canola":             "Brassica napus",
    "oilseed rape":       "Brassica napus",
    "brassica":           "Brassica napus",
    "cabbage":            "Brassica oleracea",
    "cotton":             "Gossypium hirsutum",
    "grape":              "Vitis vinifera",
    "grapevine":          "Vitis vinifera",
    "sunflower":          "Helianthus annuus",
    "banana":             "Musa acuminata",
    "cassava":            "Manihot esculenta",
    "sugarcane":          "Saccharum officinarum",
    "sugar cane":         "Saccharum officinarum",
    "strawberry":         "Fragaria ananassa",
    "apple":              "Malus domestica",
    "pear":               "Pyrus communis",
    "citrus":             "Citrus sinensis",
    "orange":             "Citrus sinensis",
    "lemon":              "Citrus limon",
    "peach":              "Prunus persica",
    "cherry":             "Prunus avium",
    "plum":               "Prunus domestica",
    "lettuce":            "Lactuca sativa",
    "cucumber":           "Cucumis sativus",
    "melon":              "Cucumis melo",
    "watermelon":         "Citrullus lanatus",
    "pumpkin":            "Cucurbita maxima",
    "squash":             "Cucurbita pepo",
    "pepper":             "Capsicum annuum",
    "chilli":             "Capsicum annuum",
    "eggplant":           "Solanum melongena",
    "aubergine":          "Solanum melongena",
    "pea":                "Pisum sativum",
    "bean":               "Phaseolus vulgaris",
    "lentil":             "Lens culinaris",
    "chickpea":           "Cicer arietinum",
    "lupin":              "Lupinus angustifolius",
    "clover":             "Trifolium repens",
    "alfalfa":            "Medicago sativa",
    "poplar":             "Populus trichocarpa",
    "eucalyptus":         "Eucalyptus grandis",
    "flax":               "Linum usitatissimum",
    "linseed":            "Linum usitatissimum",
    "coffee":             "Coffea arabica",
    "cacao":              "Theobroma cacao",
    "cocoa":              "Theobroma cacao",
    "brachypodium":       "Brachypodium distachyon",
    "setaria":            "Setaria italica",
}

FSP_RE = re.compile(r'^(\S+\s+\S+\s+f\.\s*sp\.\s*\S+)', re.IGNORECASE)


def _build_host_index(db: dict) -> list[tuple[str, str]]:
    """Return (name_lower, canonical_name) sorted longest-first."""
    host_taxids: set[int] = set(int(t) for t in db.get("host_to_seed", {}))
    taxid_to_name: dict[int, str] = {int(k): v
                                     for k, v in db.get("taxid_to_name", {}).items()}

    entries: dict[str, str] = {}   # lower → canonical

    # Scientific names from name_to_taxid filtered to host taxids
    for name, taxid in db.get("name_to_taxid", {}).items():
        if int(taxid) in host_taxids and len(name.split()) >= 2:
            canonical = taxid_to_name.get(int(taxid), name)
            entries[name.lower()] = canonical

    # Vernacular names
    for vern, canonical in VERNACULAR.items():
        entries[vern.lower()] = canonical

    return sorted(entries.items(), key=lambda x: -len(x[0]))  # longest first


def _search_text(text: str, index: list[tuple[str, str]]) -> tuple[str, str]:
    """Return (canonical_name, matched_string) or ('', '')."""
    if not text:
        return "", ""
    tl = text.lower()
    for name_lower, canonical in index:
        # word-boundary check to avoid "corn" matching "unicorn"
        pattern = r'(?<![a-z])' + re.escape(name_lower) + r'(?![a-z])'
        if re.search(pattern, tl):
            return canonical, name_lower
    return "", ""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(PHIBASE_DB) as f:
        db = json.load(f)
    host_index = _build_host_index(db)
    print(f"Host name index: {len(host_index):,} entries")

    # Load BioProject metadata
    bp_meta: dict[str, dict] = {}
    with open(META_TSV, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            bp_meta[r["BioProject"]] = r

    # Search each BioProject
    results: dict[str, dict] = {}
    source_counts: Counter = Counter()

    for bp, meta in bp_meta.items():
        title   = meta.get("title", "")
        desc    = meta.get("description", "")
        abstract = meta.get("abstract", "")
        methods = meta.get("methods_text", "")

        found, match, source = "", "", ""
        for src_name, src_text in (("title", title), ("abstract", abstract),
                                   ("description", desc), ("methods", methods)):
            found, match = _search_text(src_text, host_index)
            if found:
                source = src_name
                break

        results[bp] = {"named_host": found, "match": match, "source": source}
        if found:
            source_counts[source] += 1

    n_total   = len(bp_meta)
    n_found   = sum(1 for r in results.values() if r["named_host"])
    print(f"\nTotal BioProjects: {n_total:,}")
    print(f"Host identified:   {n_found:,}  ({n_found/n_total*100:.1f}%)")
    print("  by source:")
    for src in ("title", "abstract", "description", "methods"):
        print(f"    {src:<12}: {source_counts[src]:4d}")

    # ── Agreement with STAT host (MAL biosample_rep runs) ────────────────────
    BROAD = {
        "Viridiplantae", "Mesangiospermae", "eudicotyledons", "rosids",
        "Euphyllophyta", "Pentapetalae", "Poales", "asterids", "Gunneridae",
        "Magnoliopsida", "Lamiales", "BOP clade", "PACMAD clade", "Streptophyta",
    }

    agree = disagree = stat_only = kw_only = neither = 0
    disagree_ex: list[str] = []

    with open(CRYPT_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("mode") != "mal":
                continue
            if row.get("biosample_representative") != "True":
                continue
            bp   = row["BioProject"]
            stat = row.get("host", "").strip()
            kw   = results.get(bp, {}).get("named_host", "")

            stat_sp = stat not in BROAD and stat != "Viridiplantae" and len(stat.split()) >= 2
            kw_ok   = bool(kw)

            if stat_sp and kw_ok:
                # Compare genus (first word)
                if stat.split()[0].lower() == kw.split()[0].lower():
                    agree += 1
                else:
                    disagree += 1
                    if len(disagree_ex) < 8:
                        disagree_ex.append(f"  STAT={stat!r}  KW={kw!r}  "
                                           f"BP={bp}")
            elif stat_sp and not kw_ok:
                stat_only += 1
            elif not stat_sp and kw_ok:
                kw_only += 1
            else:
                neither += 1

    comp_total = agree + disagree + stat_only + kw_only + neither
    print(f"\nMAL biosample_rep agreement (genus level, {comp_total:,} runs):")
    print(f"  Both agree      : {agree:5d}  ({agree/comp_total*100:.1f}%)")
    print(f"  Disagree        : {disagree:5d}  ({disagree/comp_total*100:.1f}%)")
    print(f"  STAT only       : {stat_only:5d}  ({stat_only/comp_total*100:.1f}%)")
    print(f"  KW only         : {kw_only:5d}  ({kw_only/comp_total*100:.1f}%)")
    print(f"  Neither         : {neither:5d}  ({neither/comp_total*100:.1f}%)")
    if disagree_ex:
        print("  Disagree examples:")
        for ex in disagree_ex:
            print(ex)

    # ── Write TSV ─────────────────────────────────────────────────────────────
    rows = []
    for bp, meta in bp_meta.items():
        r = results.get(bp, {})
        rows.append({
            "BioProject":    bp,
            "named_host":    r.get("named_host", ""),
            "match":         r.get("match", ""),
            "source":        r.get("source", ""),
            "llm_treatment": meta.get("llm_treatment", ""),
            "llm_setting":   meta.get("llm_study_setting", ""),
            "title":         meta.get("title", ""),
        })

    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWritten: {OUT_TSV}  ({len(rows):,} rows)")

    # Top named hosts
    host_counts: Counter = Counter(
        r["named_host"] for r in results.values() if r["named_host"])
    print("\nTop named hosts:")
    for h, n in host_counts.most_common(20):
        print(f"  {n:4d}  {h}")


if __name__ == "__main__":
    main()
