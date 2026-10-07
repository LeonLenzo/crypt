#!/usr/bin/env python3
"""Split PRJNA918399 into glasshouse in-planta samples and YPD broth cultures.

    Jones DAB / Rudd JJ et al. 2023. "Combined pangenomics and transcriptomics reveals core
    and redundant virulence processes in a rapidly evolving fungal plant pathogen."
    BMC Biology 21:17.  doi:10.1186/s12915-023-01520-6   PMC9903594   open access.
    GEO GSE222164, SRA PRJNA918399.

The second deposit in this cohort that mixes infected plant tissue with axenic fungal culture,
after PRJEB8798, and the second where the obvious archive field gives the wrong answer. Here
every run is registered under the same organism, so nothing structured separates them; the
split lives in the GEO `source_name`, which is explicit:

    wheat leaves              78
    in-vitro fungal culture   36

Both arms are controlled, so the cohort outcome does not turn on this. It is written anyway
for the same reason as PRJEB8798: an axenic culture run contains no host tissue, and `setting`
is the only column that can say so, which keeps those 36 runs out of any denominator of plant
samples.

The in-planta setting comes from the paper, not the archive:

    "Z. tritici inoculation of wheat was performed in an environment-controlled GLASSHOUSE
     facility (17 degC day/night, 60% relative humidity, 16/8 h light/dark)"
    "A single susceptible wheat cultivar 'Panorama' ... Samples were harvested at 6dpi (within
     the symptomless phase) and 9dpi (transition phase ...). Each sample consisted of 3
     independent pooled seedling leaves"
    "Each isolate was also grown and harvested during log phase growth in YPD broth to
     investigate gene expression away from the wheat plant."

`SampleName` in runs.tsv is the GSM accession, so the GEO record joins with no regex.

Run:  python 01a_Literature/resolvers/ztpan23.py
"""

from __future__ import annotations

import collections, csv, re, sys
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _resolver import finish
from _layout import GEO, STUDIES

BIOPROJECT = "PRJNA918399"
BRIEF = GEO / "GSE222164_brief.txt"
OUT = STUDIES / "doi_10.1186_s12915-023-01520-6" / "ztpan23_runs.csv"

N_PLANTA, N_VITRO = 78, 36

MAP = {
    "wheat leaves": ("greenhouse", "leaf (pooled seedling leaves)", "in planta"),
    "in-vitro fungal culture": ("in vitro culture", "axenic fungal culture (no plant)", "in vitro"),
}


def main() -> None:
    if not BRIEF.exists():
        sys.exit(f"REFUSED: {BRIEF} missing. Run: sources.py geo GSE222164")
    txt = BRIEF.read_text(encoding="utf-8", errors="replace")

    rows = []
    for blk in txt.split("^SAMPLE = ")[1:]:
        gsm = blk.split("\n", 1)[0].strip()
        src = ""
        for line in blk.splitlines():
            m = re.match(r"!Sample_source_name_ch1 = (.*)", line)
            if m:
                src = m.group(1).strip()
        if src not in MAP:
            sys.exit(f"REFUSED: {gsm} has source_name {src!r}, which is not one of {list(MAP)}. "
                     f"An unmapped arm would let in-vitro runs pass as plant samples.")
        setting, tissue, arm = MAP[src]
        rows.append(dict(GSM=gsm, Setting=setting, Tissue=tissue, Arm=arm, SourceName=src))

    n_p = sum(1 for r in rows if r["Arm"] == "in planta")
    n_v = len(rows) - n_p
    if (n_p, n_v) != (N_PLANTA, N_VITRO):
        sys.exit(f"REFUSED: expected {N_PLANTA} in planta and {N_VITRO} in vitro, got {n_p} and {n_v}")

    print(f"{len(rows)} GEO samples: " +
          "  ".join(f"{k}:{v}" for k, v in collections.Counter(r["Arm"] for r in rows).most_common()),
          file=sys.stderr)
    finish(rows, OUT, key_col="GSM", sra_key_spec="SampleName",
           bioprojects=[BIOPROJECT], sort_key=lambda r: r["GSM"])


if __name__ == "__main__":
    main()
