#!/usr/bin/env python3
"""Host grouping, so the cereal priority is mechanical rather than remembered.

Leon's call 2026-10-06: the work is funded for cereals, so cereals are curated first. The
non-cereal data already curated is KEPT and labelled, not retired - it is discovery-set
material in hosts where db_v3's pathogen reference is weakest, which is the least favourable
place to look and therefore worth having.

Three groups:

    cereal   the funded subject. Grain crops: wheat, barley, rice, maize, sorghum, oats,
             rye, triticale, and the millets.
    grass    non-cereal Poaceae (Brachypodium, switchgraass, ryegrass, sugarcane). Kept
             separate because they are cereal-adjacent for reference-database purposes but
             are not the funded subject and are not grain crops.
    other    everything else.

## Why a pathogen name can imply the host

`bioprojects.organism` is the submitter's subject, and for in-planta work it is routinely the
PATHOGEN: Ada21's 903 infected wheat leaves are registered under *Puccinia striiformis*. A
cereal-first ordering that reads only the plant names would therefore sort the single best
cereal co-infection cohort we have into `other`.

So a cereal-specific pathogen maps to `cereal?` - sorted with the cereals, but marked as an
INFERENCE from the pathogen's host range rather than a statement about the library. It must
be confirmed from the paper like anything else. Only pathogens that are effectively
cereal-restricted are listed; generalists (Fusarium oxysporum, Botrytis, Rhizoctonia solani,
Sclerotinia) are deliberately absent because they imply nothing about the host.
"""

from __future__ import annotations

import re

CEREAL = re.compile(
    r"Triticum|Hordeum|Oryza|Zea\s+mays|Zea\s+(nicaraguensis|diploperennis|luxurians|perennis)"
    r"|Sorghum|Secale|Avena|Triticosecale|Triticale|x?Triticosecale"
    r"|Setaria\s+italica|Pennisetum\s+glaucum|Cenchrus\s+americanus|Eleusine\s+coracana"
    r"|Panicum\s+miliaceum|Digitaria\s+exilis|Fagopyrum"      # buckwheat: pseudocereal, grain
    r"|wheat|barley|\brice\b|maize|\bcorn\b|sorghum|\boat\b|\brye\b|millet", re.I)

GRASS = re.compile(
    r"Brachypodium|Panicum\s+virgatum|Setaria\s+viridis|Lolium|Festuca|Saccharum"
    r"|Miscanthus|Agrostis|Poa\s|Bromus|Aegilops|switchgrass|sugarcane", re.I)

# Effectively cereal-restricted pathogens. A hit means the LIBRARY is probably infected cereal
# tissue; it is an inference from host range, never a curated host.
CEREAL_PATHOGEN = re.compile(
    r"Puccinia\s+(striiformis|graminis|triticina|recondita|hordei|coronata|sorghi|polysora)"
    r"|Blumeria\s+graminis|Zymoseptoria|Mycosphaerella\s+graminicola"
    r"|Pyrenophora\s+(teres|tritici|graminea)|Drechslera\s+teres"
    r"|Parastagonospora\s+nodorum|Phaeosphaeria\s+nodorum|Septoria\s+(tritici|nodorum)"
    r"|Pyricularia\s+(oryzae|grisea)|Magnaporthe\s+(oryzae|grisea)"
    r"|Ustilago\s+(maydis|hordei|nuda|tritici)|Ustilaginoidea\s+virens|Sporisorium"
    r"|Tilletia|Claviceps\s+purpurea|Ramularia\s+collo|Rhynchosporium"
    r"|Bipolaris\s+(sorokiniana|oryzae|maydis)|Cochliobolus\s+(sativus|heterostrophus|miyabeanus)"
    r"|Setosphaeria\s+turcica|Exserohilum\s+turcicum|Cercospora\s+zeae"
    r"|Colletotrichum\s+(graminicola|sublineola)|Fusarium\s+(graminearum|pseudograminearum|verticillioides)"
    r"|Gibberella\s+zeae|Microdochium\s+nivale|Xanthomonas\s+oryzae|Rhizoctonia\s+cerealis",
    re.I)

BARLEY = re.compile(r"Hordeum|barley", re.I)
BARLEY_PATHOGEN = re.compile(
    r"Pyrenophora\s+(teres|graminea)|Drechslera\s+teres|Ramularia\s+collo"
    r"|Puccinia\s+hordei|Rhynchosporium|Ustilago\s+(hordei|nuda)"
    r"|Blumeria\s+graminis\s+f\.?\s*sp\.?\s*hordei", re.I)

ORDER = {"cereal": 0, "cereal?": 1, "grass": 2, "other": 3}


def host_group(name: str) -> str:
    """Group a HOST species name. Plant names only; no pathogen inference."""
    n = (name or "").strip()
    if not n:
        return "other"
    if CEREAL.search(n):
        return "cereal"
    if GRASS.search(n):
        return "grass"
    return "other"


def project_group(organism: str, *texts: str) -> tuple[str, str]:
    """Group a PROJECT, which may be registered under its pathogen.

    Returns (group, why). `cereal?` means the group was inferred from a cereal-restricted
    pathogen name and still needs confirming from the paper.
    """
    org = (organism or "").strip()
    g = host_group(org)
    if g != "other":
        return g, f"host name {org!r}"
    blob = " ".join([org, *[t or "" for t in texts]])
    m = CEREAL_PATHOGEN.search(blob)
    if m:
        return "cereal?", f"cereal-restricted pathogen {m.group(0)!r}; host inferred, confirm from the paper"
    m = CEREAL.search(blob)
    if m:
        return "cereal?", f"{m.group(0)!r} appears in the title, not the organism; confirm"
    return "other", f"host name {org!r}" if org else "no organism recorded"
