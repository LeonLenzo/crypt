#!/usr/bin/env python3
"""audit_setting.py — is a deposit field material, or an isolate grown from it?

    python release/01_cohort/scripts/audit_setting.py

Reads  data/outputs/supp_data_s2_cohort.tsv   the cohort and its curated setting
       data/outputs/supp_data_s1_studies.tsv  which paper generated each deposit
       data/raw/fulltext.tsv + data/raw/fulltext/*   the papers themselves
       data/raw/sample_attrs.jsonl            what the submitter recorded per sample
Writes figures/setting_audit.tsv              one row per BioProject, evidence per column

WHY THIS EXISTS. A rust paper collects a pustule from a field, then multiplies it through
rounds of inoculation on a susceptible variety in a glasshouse and sequences THAT. The
BioSample still says the country the pustule came from, often with a `geo_loc_name` and
sometimes coordinates, so the deposit looks like field sampling at every level except the
sentence in the methods. Curated as `field` it becomes a dot on a collection map that no
one ever collected from, and it enters a field-versus-controlled comparison on the wrong
side. PRJNA749614 was caught this way, by reading one BioSample; PRJEB84390's methods say
it outright.

This does not decide anything. It reads each primary paper and each deposit's BioSamples
and reports the evidence both ways, because the discriminator is a sentence in the methods
and no keyword is reliable on its own: a field survey's methods also say "inoculated" when
describing a glasshouse assay it ran alongside, and an isolate paper still says
"collected". Rank by `review` and read the quoted sentence.
"""

import collections
import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
DATA = MODULE / "data"
OUT = MODULE / "figures" / "setting_audit.tsv"

# Phrases that describe growing the thing that was sequenced. "Multiplied" and
# "successive rounds" are the strongest: nothing in a field survey needs them.
PROPAGATED = [
    r"multipli(ed|cation)", r"successive rounds? of inoculation", r"single pustule",
    r"propagat(ed|ion)", r"increased on", r"maintained on", r"urediniospores were (harvested|collected)",
    r"glasshouse", r"greenhouse", r"growth chamber", r"controlled environment",
    r"seedlings were inoculated", r"axenic", r"in vitro", r"culture[ds]? on",
]
# Phrases that describe taking the thing that was sequenced out of a field.
FIELD = [
    r"field[- ]collected", r"collected (from|in) (the )?field", r"naturally infected",
    r"field surve(y|ys)", r"field sampl(e|es|ing)", r"survey sampl(e|es|ing)",
    r"commercial (crops?|fields?)", r"farmers?.{0,3} fields?",
]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def body_for(doi, manifest):
    row = manifest.get(doi)
    if not row or not row.get("path"):
        return ""
    for base in (MODULE, DATA):
        p = base / row["path"]
        if p.is_file():
            raw = p.read_text(errors="ignore")
            return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", raw))
    return ""


def hits(text, patterns):
    found = []
    for pat in patterns:
        m = re.search(pat, text, re.I)
        if m:
            found.append((pat, text[max(0, m.start() - 110):m.start() + 160].strip()))
    return found


def main():
    cohort = read(DATA / "outputs" / "supp_data_s2_cohort.tsv")
    studies = read(DATA / "outputs" / "supp_data_s1_studies.tsv")
    manifest = {r["doi"]: r for r in read(DATA / "raw" / "fulltext.tsv")}

    owner = {}
    for s in studies:
        if s["verdict"] != "keep":
            continue
        for bp in [a.strip() for a in s["bioprojects_own"].replace(";", ",").split(",")
                   if a.strip()]:
            owner[bp] = s["doi"]

    run_sample = {r["run"]: r["sample"] for r in read(DATA / "raw" / "samples.tsv")}
    attrs = {}
    with open(DATA / "raw" / "sample_attrs.jsonl", encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            attrs[d.get("sample")] = d.get("attrs", {})

    cached, rows = {}, []
    for bp in sorted({r["BioProject"] for r in cohort}):
        runs = [r for r in cohort if r["BioProject"] == bp]
        doi = owner.get(bp, "")
        if doi not in cached:
            cached[doi] = body_for(doi, manifest)
        text = cached[doi]
        grown, wild = hits(text, PROPAGATED), hits(text, FIELD)

        # The submitter's own record. A named strain with no real collection date is the
        # shape of an isolate; a date and no strain is the shape of a field sample.
        strain = dated = 0
        for r in runs:
            a = attrs.get(run_sample.get(r["run"], ""), {})
            if a.get("strain") or a.get("isolate"):
                strain += 1
            date = (a.get("collection_date") or a.get("collection date") or "").lower()
            if date and date not in ("not provided", "missing", "na", "not applicable"):
                dated += 1

        setting = collections.Counter(r["setting"] or "blank" for r in runs)
        says_field = setting.get("field", 0)
        # Flag only what matters: curated as field, with propagation language present
        # and either no field language or the submitter's isolate signature.
        review = bool(says_field and grown and (not wild or strain and not dated))
        rows.append({
            "BioProject": bp, "runs": len(runs), "paper": doi or "(unattributed)",
            "setting": ";".join(f"{k}:{v}" for k, v in setting.most_common()),
            "has_body": int(bool(text)),
            "propagation_phrases": len(grown), "field_phrases": len(wild),
            "samples_with_strain": strain, "samples_with_date": dated,
            "review": "yes" if review else "",
            "strongest_propagation": grown[0][1][:300] if grown else "",
        })

    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    flagged = [r for r in rows if r["review"]]
    nobody = [r for r in rows if not r["has_body"]]
    print(f"-> {OUT.name}  {len(rows)} BioProjects\n")
    print(f"{len(flagged)} curated as field with propagation language and no clear "
          f"field language:\n")
    for r in sorted(flagged, key=lambda r: -r["runs"]):
        print(f"  {r['BioProject']:<14}{r['runs']:>4} runs  {r['paper']}")
        print(f"      strain {r['samples_with_strain']}/{r['runs']}, "
              f"dated {r['samples_with_date']}/{r['runs']}, "
              f"propagation {r['propagation_phrases']}, field {r['field_phrases']}")
        print(f"      \"{r['strongest_propagation'][:190]}\"\n")
    if nobody:
        print(f"{len(nobody)} BioProjects have no paper text to read, so nothing here "
              f"could check them: {', '.join(r['BioProject'] for r in nobody)}")
    print("\nA keyword is not a verdict. Read the quoted sentence before changing a "
          "setting, and record the reason in curated_attributes.tsv.")


if __name__ == "__main__":
    main()
