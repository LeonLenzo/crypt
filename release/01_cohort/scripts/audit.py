#!/usr/bin/env python3
"""audit.py — check the cohort's tables and hand-authored overlays agree with each other.

    python release/01_cohort/scripts/audit.py

Reads every table in data/ and reports only problems. Exits non-zero if any are found,
so it can gate a rebuild.

The overlays (decisions.tsv, relations.tsv, excluded_projects.tsv, curated_attributes.tsv)
are keyed by DOI or accession and nothing enforces that those keys exist. A row naming a
DOI that was never fetched, or an accession that left the cohort, does nothing and says
nothing: it is a decision that looks recorded and is not applied. Every check here exists
because that class of failure is silent in the output.
"""
import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
problems, notes = [], []


def read(path, comments=False):
    if not path.is_file():
        problems.append(f"MISSING FILE: {path.relative_to(DATA.parent)}")
        return []
    text = path.read_text(encoding="utf-8").splitlines()
    if comments:
        text = [l for l in text if l and not l.startswith("#")]
    return list(csv.DictReader(text, delimiter="\t"))


def accessions(field):
    return [a.strip() for a in (field or "").replace(";", ",").split(",") if a.strip()]


dois = {d.strip().lower() for d in (DATA / "inputs" / "pg_dois.txt").read_text().split()
        if d.strip()}
abstracts = {r["doi"].lower(): r for r in read(DATA / "raw" / "abstracts.tsv")}
fulltext = {r["doi"].lower(): r for r in read(DATA / "raw" / "fulltext.tsv")}
screen = {r["doi"].lower(): r for r in read(DATA / "curation" / "screen_fulltext.tsv")}
decisions = read(DATA / "curation" / "decisions.tsv")
relations = read(DATA / "curation" / "relations.tsv", comments=True)
excluded = read(DATA / "curation" / "excluded_projects.tsv", comments=True)
curated = read(DATA / "curation" / "curated_attributes.tsv")
cohort = read(DATA / "outputs" / "supp_data_s2_cohort.tsv")
studies = read(DATA / "outputs" / "supp_data_s1_studies.tsv")

runs = {r["run"] for r in cohort}
projects = {r["BioProject"] for r in cohort}
kept = {r["doi"].lower() for r in studies if r["verdict"] == "keep"}

# 1. corpus completeness
for doi in sorted(dois - set(abstracts)):
    problems.append(f"in pg_dois.txt but never fetched an abstract: {doi}")
bodyless = [d for d in kept if not (fulltext.get(d, {}).get("path") or "").strip()]
if bodyless:
    notes.append(f"{len(bodyless)} kept papers have no body text, so their accessions "
                 f"cannot be scanned: {', '.join(sorted(bodyless))}")

# 2. overlays point at things that exist
for row in decisions:
    if row["doi"].lower() not in dois:
        problems.append(f"decisions.tsv names a DOI not in pg_dois.txt: {row['doi']}")
for row in relations:
    if row["doi"].lower() not in dois:
        problems.append(f"relations.tsv names a DOI not in pg_dois.txt: {row['doi']}")
    elif row["doi"].lower() not in kept and row["relation"] == "generated":
        notes.append(f"relations.tsv credits {row['doi']} with {row['BioProject']}, but "
                     f"that paper is not kept, so the credit is not drawn")
    if row["BioProject"] not in projects and row["relation"] == "generated":
        notes.append(f"relations.tsv credits {row['doi']} with {row['BioProject']}, "
                     f"which is not in the cohort")
for row in excluded:
    if row["BioProject"] in projects:
        problems.append(f"excluded_projects.tsv lists {row['BioProject']} but it is "
                        f"still in the cohort")
stray = [r["run"] for r in curated if r["run"] not in runs]
if stray:
    notes.append(f"{len(stray)} curated_attributes rows are for runs not in the cohort "
                 f"(harmless, but they are curation that is doing nothing)")

# 3. every deposit has exactly one primary, every paper's claims are real
owner = collections.defaultdict(list)
for s in studies:
    if s["verdict"] != "keep":
        continue
    for bp in accessions(s["bioprojects_own"]):
        if bp in projects:
            owner[bp].append(s["doi"])
for bp in sorted(projects - set(owner)):
    problems.append(f"no paper is credited with generating {bp} "
                    f"({sum(1 for r in cohort if r['BioProject'] == bp)} runs)")
for bp, who in sorted(owner.items()):
    if len(who) > 1:
        notes.append(f"{bp} is credited to {len(who)} papers: {', '.join(who)}")

# 4. runs reach the cohort only through a kept paper
named = collections.defaultdict(set)
for s in studies:
    if s["verdict"] == "keep":
        for bp in accessions(s["bioprojects"]) + accessions(s["bioprojects_own"]):
            named[bp].add(s["doi"])
for bp in sorted(projects - set(named)):
    problems.append(f"{bp} is in the cohort but no kept paper names it")

# 5. the figure's own predicate
blank = [r for r in cohort if not r["setting"]]
if blank:
    byproj = collections.Counter(r["BioProject"] for r in blank)
    notes.append(f"{len(blank)} runs have no setting, so they cannot count as field "
                 f"wheat whatever their host: {dict(byproj)}")
wheat_no_setting = [r for r in blank if r["host"].startswith("Triticum")]
if wheat_no_setting:
    problems.append(f"{len(wheat_no_setting)} runs are Triticum with no setting, so they "
                    f"are drawn as not-field-wheat when the truth is unknown: "
                    f"{dict(collections.Counter(r['BioProject'] for r in wheat_no_setting))}")

print(f"cohort: {len(cohort)} runs, {len(projects)} BioProjects, {len(kept)} kept papers")
print(f"overlays: {len(decisions)} decisions, {len(relations)} relations, "
      f"{len(excluded)} excluded projects, {len(curated)} curated runs\n")
for n in notes:
    print(f"  note    {n}")
for p in problems:
    print(f"  PROBLEM {p}")
print(f"\n{len(problems)} problem(s), {len(notes)} note(s)")
sys.exit(1 if problems else 0)
