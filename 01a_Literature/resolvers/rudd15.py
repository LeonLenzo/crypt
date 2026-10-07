#!/usr/bin/env python3
"""Split PRJEB8798 into in-planta glasshouse samples and axenic broth cultures.

    Rudd JJ et al. 2015. "Transcriptome and metabolite profiling of the infection cycle of
    Zymoseptoria tritici on wheat reveals a biphasic interaction with plant immunity..."
    Plant Physiology 167:1158-1185.  doi:10.1104/pp.114.255927   PDF supplied by leon.

51 runs, and they are not all the same kind of sample. Six are the fungus grown in a shake
flask with no plant present at all, which is not a plant RNA-seq sample on any reading:

    "Fungal cultures were propagated in shake flasks (220 rpm) at 18 degC for 3 d (for PDB)
     or 5 d (for CDB) and then harvested via filtration."

The other 45 are an infection time course in a glasshouse:

    "Plant inoculation experiments were done as described previously (Keon et al., 2007) using
     a spore density of 1 x 10^6 spores mL-1 in 0.01% (v/v) Tween 20 in sterile water. Mock
     inoculations of plants were made using just the Tween 20 water solution. Each biological
     replicate plant sample ... was made up of five leaves collected from five independent
     plants randomly distributed in a single walk-in temperature- and humidity-controlled
     GLASSHOUSE."

## Why the archive cannot do this on its own, and where the key is

`scientific_name` splits the project 30 Triticum aestivum / 21 Zymoseptoria tritici IPO323,
which is NOT the split that matters: five of the seven fungal conditions are the same infected
wheat leaves registered under the pathogen. The usable key is the ENA `sample_title`, which is
explicit - `I_9_2`, `M_21_3`, `Z.tritici on wheat leafs (9dpi)`, `Z.tritici grown on potato
dextrose broth (PDB)` - and is cached by `sources.py ena`. So:

    in planta   30 wheat-registered (I_ inoculated and M_ mock, at 1/4/9/14/21 dpi x 3)
              + 15 fungus-registered on wheat leaves at the same five timepoints  = 45
    in vitro     6 broth cultures, CDB x3 and PDB x3

Nothing here is field, so the cohort outcome is the same either way. The split is still worth
writing: an axenic culture run has no host tissue, so it must never sit in a denominator of
plant samples, and `setting` is the only column that can say so.

Run:  python 01a_Literature/resolvers/rudd15.py
"""

from __future__ import annotations

import collections, csv, sys
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _resolver import finish
from _layout import ENA, STUDIES

BIOPROJECT = "PRJEB8798"
REPORT = ENA / f"{BIOPROJECT}.tsv"
OUT = STUDIES / "doi_10.1104_pp.114.255927" / "rudd15_runs.csv"

N_INVITRO, N_INPLANTA = 6, 45


def classify(title: str) -> tuple[str, str, str]:
    """(setting, tissue, arm) from the ENA sample title."""
    t = (title or "").strip()
    low = t.lower()
    if "broth" in low or "culture medium" in low:
        return "in vitro culture", "axenic fungal culture (no plant)", "in vitro"
    if "on wheat leaf" in low or "on wheat leafs" in low:
        return "greenhouse", "leaf", "in planta (fungus-registered)"
    if t[:2] in ("I_", "M_"):
        return "greenhouse", "leaf", "in planta (mock)" if t.startswith("M_") else "in planta (inoculated)"
    return "", "", ""


def main() -> None:
    if not REPORT.exists():
        sys.exit(f"REFUSED: {REPORT} missing. Run: sources.py ena {BIOPROJECT}")
    with open(REPORT) as fh:
        recs = list(csv.DictReader(fh, delimiter="\t"))

    rows = []
    for r in recs:
        setting, tissue, arm = classify(r.get("sample_title", ""))
        if not setting:
            sys.exit(f"REFUSED: unclassifiable sample_title {r.get('sample_title')!r} on "
                     f"{r['run_accession']}. Every title must map, or the in-vitro runs leak "
                     f"into the plant samples silently.")
        rows.append(dict(Run=r["run_accession"], Setting=setting, Tissue=tissue, Arm=arm,
                         SampleTitle=r.get("sample_title", ""),
                         ScientificName=r.get("scientific_name", "")))

    arms = collections.Counter(r["Arm"] for r in rows)
    n_vitro = sum(1 for r in rows if r["Arm"] == "in vitro")
    if n_vitro != N_INVITRO or len(rows) - n_vitro != N_INPLANTA:
        sys.exit(f"REFUSED: expected {N_INVITRO} in vitro and {N_INPLANTA} in planta, "
                 f"got {n_vitro} and {len(rows) - n_vitro}")

    print(f"{len(rows)} runs: " + "  ".join(f"{k}:{v}" for k, v in arms.most_common()),
          file=sys.stderr)
    finish(rows, OUT, key_col="Run", bioprojects=[BIOPROJECT], sort_key=lambda r: r["Run"])


if __name__ == "__main__":
    main()
