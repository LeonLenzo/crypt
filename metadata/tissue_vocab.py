#!/usr/bin/env python3
"""
tissue_vocab.py — plant tissue field normalisation + aerial BioSample fetch.

Two modes:

  --survey   Query NCBI BioSample for plant RNA-seq tissue values and output a
             sorted frequency table.  Useful for expanding TISSUE_MAP.
             Output: metadata/output/aerial/tissue_vocab.tsv

  --fetch    Fetch all plant RNA-seq BioSamples that have a tissue attribute,
             classify each with TISSUE_MAP, and save aerial-tissue records.
             Output: metadata/output/aerial/biosamples.tsv
                     metadata/output/aerial/category_counts.tsv

Run from crypt/ root:
    python metadata/tissue_vocab.py --survey
    python metadata/tissue_vocab.py --fetch
"""

import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import http_get

# ── Config ────────────────────────────────────────────────────────────────────

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL       = "leon.lenzo@curtin.edu.au"
TOOL        = "crypt/tissue_vocab"
BATCH       = 500
API_KEY     = os.environ.get("NCBI_API_KEY", "")
RATE        = 0.12 if API_KEY else 0.4
OUT_DIR     = Path("metadata/output/aerial")

# BioSample attribute names that hold tissue information
_TISSUE_ATTRS = {"tissue", "tissue type", "tissue_type", "plant_part", "plant part"}

# ── Normalisation table ───────────────────────────────────────────────────────
# Maps lowercase-stripped tissue string → category.
# Categories: leaf | aerial_other | root | seed_fruit | whole_plant | pathogen

TISSUE_MAP: dict = {
    # ── leaf ──────────────────────────────────────────────────────────────────
    "leaf":                                  "leaf",
    "leaves":                                "leaf",
    "leave":                                 "leaf",
    "young leaf":                            "leaf",
    "young leaves":                          "leaf",
    "seedling leaf":                         "leaf",
    "leaf tissue":                           "leaf",
    "leaf sample":                           "leaf",
    "plant leaf":                            "leaf",
    "blade":                                 "leaf",
    "leaf blade":                            "leaf",
    "leaf or winter bud":                    "leaf",
    "leaf and/or bud":                       "leaf",
    "leaf and/or young floral buds":         "leaf",
    "one bud and two leaves":                "leaf",
    "bulk of young leaves":                  "leaf",
    "true leaf":                             "leaf",
    "herbarium leaf":                        "leaf",
    "infected leaf":                         "leaf",
    "infected leaves":                       "leaf",
    "leafs":                                 "leaf",
    "fresh leaf tissue":                     "leaf",
    "tender leaf":                           "leaf",
    "tender buds and leaves":                "leaf",
    "leaf head":                             "leaf",
    "leafy-head":                            "leaf",
    "leaf material":                         "leaf",
    "flag leaf":                             "leaf",
    "flag leaves":                           "leaf",
    "mature leaf":                           "leaf",
    "mature leaves":                         "leaf",
    "old leaf":                              "leaf",
    "second leaf":                           "leaf",
    "first leaf":                            "leaf",
    "first leaf segment":                    "leaf",
    "first leave segment":                   "leaf",
    "first true leaf":                       "leaf",
    "first true leaves":                     "leaf",
    "juvenile leaf":                         "leaf",
    "leaf disk":                             "leaf",
    "leaf disc":                             "leaf",
    "leaflet":                               "leaf",
    "leaflets":                              "leaf",
    "leaf lamina":                           "leaf",
    "leaf sheath":                           "leaf",
    "leaf petiole":                          "leaf",
    "leaf midrib and petiole":               "leaf",
    "leaf tissues":                          "leaf",
    "leaf and stem":                         "leaf",
    "leaf tissue from v4 stage":             "leaf",
    "the second leaves":                     "leaf",
    "rosette leaves":                        "leaf",
    "rosette leaf":                          "leaf",
    "5-week old leaves":                     "leaf",
    "detached leaf":                         "leaf",
    "inoculated leaf":                       "leaf",
    "secondary vein, leaf 10":               "leaf",
    "shoot apex region containing one unexpanded leaf": "leaf",
    "aerial rosette tissue at 21/22 das (days after sowing)": "leaf",
    "shoots and leaves":                     "leaf",
    "stems and leaves":                      "leaf",
    "shoot/leaves":                          "leaf",
    # ── aerial other ──────────────────────────────────────────────────────────
    "stem":                                  "aerial_other",
    "stems":                                 "aerial_other",
    "shoot":                                 "aerial_other",
    "shoots":                                "aerial_other",
    "young shoot":                           "aerial_other",
    "young shoots":                          "aerial_other",
    "shoot tip":                             "aerial_other",
    "shoot tips":                            "aerial_other",
    "main shoot tip":                        "aerial_other",
    "microshoot":                            "aerial_other",
    "aboveground tissue":                    "aerial_other",
    "aboveground":                           "aerial_other",
    "aerial":                                "aerial_other",
    "aerial parts":                          "aerial_other",
    "aerial part":                           "aerial_other",
    "aerial_part":                           "aerial_other",
    "above-ground":                          "aerial_other",
    "vegetative":                            "aerial_other",
    "apex":                                  "aerial_other",
    "shoot apex":                            "aerial_other",
    "apical shoot":                          "aerial_other",
    "apical meristem":                       "aerial_other",
    "meristem":                              "aerial_other",
    "shoot apical meristem (sam)":           "aerial_other",
    "sam":                                   "aerial_other",
    "sam and surrounding tissue":            "aerial_other",
    "spike":                                 "aerial_other",
    "spikelets":                             "aerial_other",
    "spikelet":                              "aerial_other",
    "rachis":                                "aerial_other",
    "lemma":                                 "aerial_other",
    "glume":                                 "aerial_other",
    "floret":                                "aerial_other",
    "florets":                               "aerial_other",
    "tassel":                                "aerial_other",
    "ear":                                   "aerial_other",
    "ear primordia":                         "aerial_other",
    "ear bract primordium":                  "aerial_other",
    "young ear":                             "aerial_other",
    "panicle":                               "aerial_other",
    "panicles":                              "aerial_other",
    "young panicle":                         "aerial_other",
    "inflorescence":                         "aerial_other",
    "flower":                                "aerial_other",
    "flower bud":                            "aerial_other",
    "floral bud":                            "aerial_other",
    "floral style":                          "aerial_other",
    "terminal flower bud":                   "aerial_other",
    "young flower bud":                      "aerial_other",
    "bud":                                   "aerial_other",
    "buds":                                  "aerial_other",
    "l1 tiller buds":                        "aerial_other",
    "rosette":                               "aerial_other",
    "rosettes":                              "aerial_other",
    "branch":                                "aerial_other",
    "middle branch tip":                     "aerial_other",
    "top branch tip":                        "aerial_other",
    "bottom branch tip":                     "aerial_other",
    "2ndary branch tip":                     "aerial_other",
    "bark":                                  "aerial_other",
    "needle":                                "aerial_other",
    "needles":                               "aerial_other",
    "needle tissue":                         "aerial_other",
    "anther":                                "aerial_other",
    "petal":                                 "aerial_other",
    "sepal":                                 "aerial_other",
    "tepals":                                "aerial_other",
    "corolla limb":                          "aerial_other",
    "pistil":                                "aerial_other",
    "stamen":                                "aerial_other",
    "ovary":                                 "aerial_other",
    "male strobili":                         "aerial_other",
    "female strobili":                       "aerial_other",
    "petiole":                               "aerial_other",
    "internode":                             "aerial_other",
    "pith":                                  "aerial_other",
    "tendril":                               "aerial_other",
    "pollen":                                "aerial_other",
    "pollen tube":                           "aerial_other",
    "pollen tubes":                          "aerial_other",
    "pollen grain":                          "aerial_other",
    "mature pollen":                         "aerial_other",
    "trichome":                              "aerial_other",
    "phyllosphere":                          "aerial_other",
    "hypocotyl":                             "aerial_other",
    "epicotyl":                              "aerial_other",
    "coleoptile":                            "aerial_other",
    "scion":                                 "aerial_other",
    "gall":                                  "aerial_other",
    "pseudobulb":                            "aerial_other",
    "meiocyte":                              "aerial_other",
    "dormant bud":                           "aerial_other",
    "whole head":                            "aerial_other",
    "head":                                  "aerial_other",
    "floral cavity":                         "aerial_other",
    "haustoria":                             "aerial_other",
    "haustorial tissue":                     "aerial_other",
    "haustoria, mixed tissues of plant leaves and fungal mycelia": "aerial_other",
    "urediniospores, soybean leaf":          "aerial_other",
    "trunk":                                 "aerial_other",
    "sheath":                                "aerial_other",
    "leaf sheaths":                          "leaf",
    "stem tissue":                           "aerial_other",
    "stem tissues":                          "aerial_other",
    "wheat spikelet":                        "aerial_other",
    "leaves, pods":                          "aerial_other",
    "leaf, stem":                            "leaf",
    # ── root ──────────────────────────────────────────────────────────────────
    "root":                                  "root",
    "roots":                                 "root",
    "root tissue":                           "root",
    "tap root":                              "root",
    "lateral root":                          "root",
    "primary root":                          "root",
    "root xylem":                            "root",
    "root tip":                              "root",
    "root cap":                              "root",
    "root hair":                             "root",
    "cortex":                                "root",
    "cortex of primary root":                "root",
    "underground tissue":                    "root",
    "rhizome":                               "root",
    "rhizosphere":                           "root",
    "tuber":                                 "root",
    "bulb":                                  "root",
    "corm":                                  "root",
    "stolon":                                "root",
    "rootstock":                             "root",
    "elongation zone":                       "root",
    "whole root":                            "root",
    "whole roots":                           "root",
    "roots, fungal tissue":                  "root",
    "hairy roots":                           "root",
    "maize root":                            "root",
    "plant root":                            "root",
    "plantlet root":                         "root",
    "roots from soybean":                    "root",
    "10 day old root":                       "root",
    "young leaves, roots, and cambial tissues": "root",
    # ── seed / fruit ──────────────────────────────────────────────────────────
    "seed":                                  "seed_fruit",
    "seeds":                                 "seed_fruit",
    "grain":                                 "seed_fruit",
    "kernel":                                "seed_fruit",
    "endosperm":                             "seed_fruit",
    "seed coat":                             "seed_fruit",
    "seed embryo":                           "seed_fruit",
    "fruit":                                 "seed_fruit",
    "fruits":                                "seed_fruit",
    "fruit flesh":                           "seed_fruit",
    "fruit core":                            "seed_fruit",
    "fruit (pericarp)":                      "seed_fruit",
    "fruit exocarp":                         "seed_fruit",
    "fruit of 10 days after pollination":    "seed_fruit",
    "fruit of 18 days after pollination":    "seed_fruit",
    "fruit of 26 days after pollination":    "seed_fruit",
    "fruit of 34 days after pollination":    "seed_fruit",
    "fruit of 42 days after pollination":    "seed_fruit",
    "berry":                                 "seed_fruit",
    "deseeded_berry":                        "seed_fruit",
    "drupe":                                 "seed_fruit",
    "mesocarp":                              "seed_fruit",
    "mesocarp-enriched tissues":             "seed_fruit",
    "flesh":                                 "seed_fruit",
    "peel":                                  "seed_fruit",
    "pod":                                   "seed_fruit",
    "pod wall":                              "seed_fruit",
    "hull":                                  "seed_fruit",
    "pericarp":                              "seed_fruit",
    "ovule":                                 "seed_fruit",
    "ovules":                                "seed_fruit",
    "outer integument":                      "seed_fruit",
    "testa":                                 "seed_fruit",
    "placenta":                              "seed_fruit",
    "silique":                               "seed_fruit",
    "mature dry seeds":                      "seed_fruit",
    "mature embryo":                         "seed_fruit",
    "embryos":                               "seed_fruit",
    "embryonic":                             "seed_fruit",
    "in vitro cotyledon":                    "seed_fruit",
    "kernels":                               "seed_fruit",
    "berries":                               "seed_fruit",
    "juice":                                 "seed_fruit",
    # ── whole plant ───────────────────────────────────────────────────────────
    "whole plant":                           "whole_plant",
    "whole plants":                          "whole_plant",
    "whole body":                            "whole_plant",
    "plants":                                "whole_plant",
    "seedling":                              "whole_plant",
    "seedlings":                             "whole_plant",
    "whole seedling":                        "whole_plant",
    "whole seedlings":                       "whole_plant",
    "fresh whole seedlings":                 "whole_plant",
    "plant":                                 "whole_plant",
    "whole cell":                            "whole_plant",
    "embryo":                                "whole_plant",
    "callus":                                "whole_plant",
    "thallus":                               "whole_plant",
    "cotyledon":                             "whole_plant",
    "cotyledons":                            "whole_plant",
    "root, shoot":                           "whole_plant",
    # ── pathogen material ─────────────────────────────────────────────────────
    "mycelium":                              "pathogen",
    "mycelia":                               "pathogen",
    "hyphae":                                "pathogen",
    "spore":                                 "pathogen",
    "conidia":                               "pathogen",
    "conidium":                              "pathogen",
    "microconidia":                          "pathogen",
    "macroconidia":                          "pathogen",
    "urediniospore":                         "pathogen",
    "urediniospores":                        "pathogen",
    "teliospore":                            "pathogen",
    "spores":                                "pathogen",
}

# Categories considered "aerial" for the fetch mode
AERIAL = {"leaf", "aerial_other"}


_TISSUE_PRIORITY = ["leaf", "aerial_other", "seed_fruit", "root", "whole_plant", "pathogen"]

def normalise_tissue(raw: str) -> str:
    """Map a raw tissue string → category. Returns 'unknown' if not in table.

    Falls back to splitting compound strings on commas/semicolons and returning
    the highest-priority category found across parts.
    """
    key = raw.lower().strip()
    if key in TISSUE_MAP:
        return TISSUE_MAP[key]
    # Compound strings e.g. "leaf, stem", "spike, rachis", "Potato Leaf, Tuber"
    parts = [p.strip() for p in key.replace(";", ",").split(",") if p.strip()]
    cats  = [TISSUE_MAP[p] for p in parts if p in TISSUE_MAP]
    for cat in _TISSUE_PRIORITY:
        if cat in cats:
            return cat
    return "unknown"


# ── Entrez helpers ────────────────────────────────────────────────────────────

def _api_suffix() -> str:
    return f"&api_key={API_KEY}" if API_KEY else ""


def _esearch(query: str, db: str = "biosample", retmax: int = 100_000) -> list:
    url = (f"{ENTREZ_BASE}/esearch.fcgi"
           f"?db={db}&term={urllib.parse.quote(query)}"
           f"&retmax={retmax}&retmode=json&email={EMAIL}&tool={TOOL}"
           f"{_api_suffix()}")
    data = json.loads(http_get(url, headers={"User-Agent": TOOL}))
    ids  = data["esearchresult"]["idlist"]
    total = int(data["esearchresult"]["count"])
    print(f"  Query hits: {total:,}  (fetching {len(ids):,})", flush=True)
    return ids


def _efetch_xml(ids: list, db: str = "biosample") -> str:
    params = {"db": db, "id": ",".join(ids),
              "rettype": "xml", "retmode": "xml",
              "email": EMAIL, "tool": TOOL}
    if API_KEY:
        params["api_key"] = API_KEY
    body = urllib.parse.urlencode(params).encode()
    req  = urllib.request.Request(
        f"{ENTREZ_BASE}/efetch.fcgi", data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode(errors="replace")
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def _tissue_from_bs(bs_elem) -> str:
    for attr in bs_elem.iter("Attribute"):
        name = (attr.get("attribute_name") or attr.get("harmonized_name") or "").lower()
        if name in _TISSUE_ATTRS:
            return (attr.text or "").strip()
    return ""


# ── Survey mode ───────────────────────────────────────────────────────────────

def run_survey(args):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    query = 'txid33090[Organism:exp] AND "tissue"[Attribute Name]'
    if args.rna_only:
        query += ' AND "rna seq"[Strategy]'
    print(f"Survey query: {query}", flush=True)
    ids = _esearch(query, retmax=args.limit)
    time.sleep(RATE)

    counts: Counter = Counter()
    n_batches = (len(ids) + BATCH - 1) // BATCH
    for i in range(0, len(ids), BATCH):
        batch = ids[i:i + BATCH]
        print(f"  batch {i//BATCH+1}/{n_batches}...", end="\r", flush=True)
        xml_text = _efetch_xml(batch)
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            time.sleep(RATE)
            continue
        for bs in root.iter("BioSample"):
            raw = _tissue_from_bs(bs)
            if raw:
                counts[raw] += 1
        time.sleep(RATE)

    out = OUT_DIR / "tissue_vocab.tsv"
    with open(out, "w") as f:
        f.write("count\ttissue\n")
        for val, n in counts.most_common():
            f.write(f"{n}\t{val}\n")
    print(f"\n{sum(counts.values()):,} values, {len(counts):,} unique  →  {out}")


# ── Fetch mode ────────────────────────────────────────────────────────────────

def run_fetch(args):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    query = ('txid33090[Organism:exp] AND "tissue"[Attribute Name]'
             ' AND "rna seq"[Strategy]')
    print(f"Fetch query: {query}", flush=True)
    ids = _esearch(query, retmax=100_000)
    time.sleep(RATE)

    cat_counts: Counter = Counter()
    aerial_rows: list   = []

    n_batches = (len(ids) + BATCH - 1) // BATCH
    for i in range(0, len(ids), BATCH):
        batch = ids[i:i + BATCH]
        print(f"  batch {i//BATCH+1}/{n_batches}...", end="\r", flush=True)
        xml_text = _efetch_xml(batch)
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            time.sleep(RATE)
            continue
        for bs in root.iter("BioSample"):
            acc = bs.get("accession", "")
            raw = _tissue_from_bs(bs)
            cat = normalise_tissue(raw) if raw else "unknown"
            cat_counts[cat] += 1
            if cat in AERIAL:
                aerial_rows.append((acc, raw, cat))
        time.sleep(RATE)

    # Save aerial BioSamples
    bs_out = OUT_DIR / "biosamples.tsv"
    with open(bs_out, "w") as f:
        f.write("accession\ttissue_raw\ttissue_category\n")
        for acc, raw, cat in aerial_rows:
            f.write(f"{acc}\t{raw}\t{cat}\n")

    # Save category breakdown
    tot = sum(cat_counts.values())
    counts_out = OUT_DIR / "category_counts.tsv"
    with open(counts_out, "w") as f:
        f.write("category\tcount\tpct\n")
        for cat, n in cat_counts.most_common():
            f.write(f"{cat}\t{n}\t{100*n/tot:.1f}\n")

    print(f"\nCategory breakdown ({tot:,} total):")
    for cat, n in cat_counts.most_common():
        print(f"  {cat:<15} {n:>6,}  ({100*n/tot:.1f}%)")
    print(f"\nAerial BioSamples: {len(aerial_rows):,}")
    print(f"  → {bs_out}")
    print(f"  → {counts_out}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--survey", action="store_true",
                      help="Enumerate tissue vocabulary (frequency table)")
    mode.add_argument("--fetch",  action="store_true",
                      help="Fetch plant RNA-seq BioSamples and classify by tissue")
    ap.add_argument("--rna-only", action="store_true",
                    help="With --survey: restrict to RNA-Seq library strategy")
    ap.add_argument("--limit", type=int, default=50_000,
                    help="Max UIDs to fetch for --survey (default: 50,000)")
    args = ap.parse_args()

    if args.survey:
        run_survey(args)
    else:
        run_fetch(args)


if __name__ == "__main__":
    main()
