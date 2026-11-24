#!/usr/bin/env python3
"""classify_lab_material.py — what each LAB library was actually made from.

    python release/01_cohort/scripts/classify_lab_material.py

Reads  data/outputs/supp_data_s2_cohort.tsv   which runs are setting=lab
       data/raw/run_titles.tsv                the submitter's description of each library
Writes data/curation/lab_material.tsv         one row per lab run, with the title as evidence

## Why the lab runs need splitting at all

They are the chapter's CONTROLS, and `setting = lab` lumps four different experiments into
one word. A germinated-spore library is pure fungal material and tests the false-positive
spectrum; a 9-dpi infected leaf is a controlled single infection and tests recovery against
a known answer; a haustorium prep is fungal-enriched FROM plant tissue, so host reads are
expected; and an uninoculated plant library is the host-only floor. Reading all 86 as "the
controls" and picking whichever looked cleanest would be choosing the answer.

## The rule, and why the title is evidence rather than a field

`host` in the BioSample record describes the ISOLATE, not the library: every one of
PRJNA396589's 15 runs says `host = wheat`, and six of them are spores that never touched a
plant (crypt-isolate-vs-sample-trap). The run TITLE is the submitter's own words about the
library, so it is a statement about the material. It is still the archive, and the archive
has misled this project in nine documented ways, so a title is a triage signal the paper
overrides, and anything the rule cannot read is written out as `unresolved` rather than
guessed into a bucket.

    spores           germinated or ungerminated urediniospores; no plant in the library
    haustoria        haustorium-enriched from infected tissue; plant reads expected
    infected_tissue  a controlled inoculation sampled in planta
    axenic_culture   grown on medium. NOTE that this can only ever be the Fusarium run:
                     both rust papers state the obvious reason, that rusts "cannot be
                     axenically cultured in the laboratory", so the pure-fungal rust
                     control is a SPORE library and not a culture. The distinction belongs
                     in the caption: spores are a resting stage, not a growing colony.
    host_only        plant material with no pathogen applied
    unresolved       the title does not say; read the paper
"""

import csv
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
COHORT = HERE / "data/outputs/supp_data_s2_cohort.tsv"
TITLES = HERE / "data/raw/run_titles.tsv"
OUT = HERE / "data/curation" / "lab_material.tsv"

# Ordered: the first rule that matches wins, and haustoria is tested BEFORE spores because a
# haustorium prep is often titled with both words.
RULES = [
    ("haustoria",       re.compile(r"haustor", re.I)),
    ("spores",          re.compile(r"\bspores?\b|urediniospore", re.I)),
    ("infected_tissue", re.compile(r"\bdpi\b|infect|plant sample from", re.I)),
]
# Projects whose whole lab arm is settled by something other than the title. Both are read
# off the deposit and both are checked against the paper in declared_pathogens.tsv.
# Runs the title cannot settle and the paper can. Each carries the sentence it rests on.
BY_RUN = {
    **{r: ("infected_tissue",
           "10.1093/gbe/evx241, Methods: \"Infection assays were undertaken with ten "
           "technical replicates on seedlings under controlled environmental conditions.\"")
       for r in ("ERR1607740", "ERR1607741", "ERR1607742", "ERR1607743")},
    "SRR1583659": ("infected_tissue",
                   "10.1186/s13059-015-0590-8, Methods: \"Tests were carried out on "
                   "seedlings under controlled environment conditions\"; the project's "
                   "other six lab runs are titled 'Plant sample from'"),
}
BY_PROJECT = {
    "PRJNA486288": ("axenic_culture",
                    "BioSample host = PDA, i.e. potato dextrose agar, the axenic arm of the "
                    "study's design"),
    "PRJNA722004": ("host_only",
                    "a wheat virome survey's greenhouse wheat; no fungal pathogen applied"),
}


def main():
    with open(COHORT) as fh:
        cohort = {r["run"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(TITLES) as fh:
        titles = {r["run"]: r for r in csv.DictReader(fh, delimiter="\t")}

    rows = []
    for run, c in cohort.items():
        if c["setting"] != "lab":
            continue
        t = titles.get(run, {})
        text = " ".join(x for x in (t.get("experiment_title"), t.get("sample_title")) if x)
        if run in BY_RUN:
            material, note = BY_RUN[run]
        elif c["BioProject"] in BY_PROJECT:
            material, note = BY_PROJECT[c["BioProject"]]
        else:
            material, note = "unresolved", ""
            for name, rx in RULES:
                if rx.search(text):
                    material, note = name, f"title matches /{rx.pattern}/"
                    break
        rows.append({"run": run, "BioProject": c["BioProject"], "material": material,
                     "basis": note, "run_title": text})

    rows.sort(key=lambda r: (r["material"], r["BioProject"], r["run"]))
    with open(OUT, "w", newline="") as fh:
        fh.write("# DERIVED by scripts/classify_lab_material.py from the run titles, then\n"
                 "# checked against the papers. `unresolved` is a worklist, not a category.\n")
        w = csv.DictWriter(fh, ["run", "BioProject", "material", "basis", "run_title"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    c = Counter(r["material"] for r in rows)
    print(f"-> {OUT}  {len(rows)} lab runs")
    for m, n in c.most_common():
        print(f"   {m:16s} {n:3d}  "
              f"{sorted({r['BioProject'] for r in rows if r['material'] == m})}")
    un = [r for r in rows if r["material"] == "unresolved"]
    if un:
        print(f"   UNRESOLVED, read the paper: " +
              ", ".join(f"{r['run']} ({r['BioProject']})" for r in un))


if __name__ == "__main__":
    main()
