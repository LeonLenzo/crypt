#!/usr/bin/env python3
"""Add the per-project setting (and any host/location the archive lacks) read from each paper.

    python research/01_literature/retrieve/seed_paper_settings.py

Reads  release/01_cohort/curation/screen_samples.tsv   runs, their project, date, what is filled
Writes release/01_cohort/curation/supplements.tsv       merges `paper:` / `sibling:` rows in

WHY. The single-project papers in the cohort carry no per-accession table to join, so their
one archive-missing fact, `setting` (field vs glasshouse vs lab), is read from the Methods and
recorded here against every run in the project. The quoted sentence that settles each call
lives in SPECS below, so a reviewer sees the basis without re-reading the paper (the audit
chain). Runs whose project the two rust supplements already resolved are left to those rows;
runs those tables missed inherit the setting their project's resolved siblings carry.

This complements seed_supplement_rows.py, which owns the `supp:` rows. Each script rewrites
only its own rows (by source prefix) and preserves the other's, so order does not matter.
"""

import collections
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
COHORT = HERE.parent.parent.parent / "release" / "01_cohort"
SCREEN = COHORT / "curation" / "screen_samples.tsv"
OUT = COHORT / "curation" / "supplements.tsv"
COLS = ["run", "host", "tissue", "setting", "location", "collection_date", "source"]

# One entry per single-project paper. `setting` is the project-wide call; `date_split` overrides
# it per run by collection-year when a project mixes historic isolates with field samples. `host`
# and `location` are filled ONLY where the archive is blank (screen_samples shows the gap).
# `evidence` is the verbatim sentence that settles the setting. `doi` becomes the row's source.
SPECS = {
    "PRJNA918968": dict(
        doi="10.15212/ijafr-2023-0110", setting="field",
        evidence="winter barley crops throughout eight different counties in Ireland in 2021 "
                 "and 2022 ... Three non-overlapping transects were sampled per field"),
    "PRJNA558380": dict(
        doi="10.1094/pdis-08-19-1724-re", setting="field",
        evidence="incidence for the most prevalent viruses across field sites"),
    "PRJNA722004": dict(
        doi="10.3389/fmicb.2021.699078", setting="field", host="Triticum aestivum",
        date_split=(2019, "field", "lab"),  # 2019+ field; earlier = historic lab isolate
        evidence="28 historic and field samples were used for total RNA sequencing"),
    "PRJEB77002": dict(
        doi="10.3389/fpls.2024.1464454", setting="field", host="Triticum aestivum",
        location="Rimski Sancevi, Novi Sad, Serbia",
        evidence="26 Pst-infected wheat samples were collected from the field trials in "
                 "2022-2023 and were subjected to RNA extraction"),
    "PRJNA628830": dict(
        doi="10.1016/j.cropro.2020.105464", setting="field", host="Triticum aestivum",
        evidence="wheat plants exhibiting characteristic symptoms of the disease were "
                 "collected during 2018 from farmers' fields and/or experimental plots"),
    "PRJNA486288": dict(
        doi="10.1080/00275514.2019.1607135", setting="field",
        evidence="we implemented a field pathogenomics strategy where RNA of naturally "
                 "infected samples was collected from multiple locations in Illinois"),
    "PRJNA1250571": dict(
        doi="10.3389/fagro.2025.1721047", setting="field", host="Triticum aestivum",
        location="Alabama, USA",
        evidence="In April 2023, four winter wheat leaf samples were collected from mature "
                 "plants at harvest-ready growth stages in Alabama"),
    "PRJNA749614": dict(
        doi="10.1101/2024.12.10.627719", setting="field", host="Triticum aestivum",
        evidence="hundreds of field isolates ... (companion deposit to PRJEB39201, whose "
                 "runs the s12864 table resolves to field)"),
    "PRJEB47693": dict(
        doi="10.1111/ppa.13532", setting="field", location="United Kingdom",
        # The paper mixes 24 field + 3 controlled-infection samples, but the 4 runs this cohort
        # lacks a setting for are exactly the UK-coded field samples (UK-02..05); the controlled
        # samples are not among our 27, and the other 23 are already curated as field.
        evidence="P. graminis-infected ... samples were collected in Ireland, France, Italy "
                 "and the UK between 2018 to 2020"),
    "PRJDB18040": dict(
        doi="10.3390/genes15060697", setting="field",
        evidence="All plants were sown in a field infected with BaYMV in Takanezawa, Tochigi "
                 "Prefecture, Japan, at the end of October 2019"),
    "PRJNA816710": dict(
        doi="10.1007/s00705-022-05649-7", setting="field", host="Triticum aestivum",
        evidence="During the 2019 winter wheat season, symptoms of severe chlorosis and "
                 "stunting were observed on wheat in the irrigation production areas of South Africa"),
    "PRJNA776691": dict(
        doi="10.3390/plants12122374", setting="field",
        evidence="Wheat (Triticum aestivum L. var. Baekkang) showing yellowing and mosaic "
                 "symptoms were collected from three regions ... in Korea in May 2022"),
}


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def year(s):
    for tok in (s or "").replace("-", " ").replace("/", " ").split():
        if tok.isdigit() and len(tok) == 4:
            return int(tok)
    return None


def main():
    runs = read(SCREEN)
    by_proj = collections.defaultdict(list)
    for r in runs:
        by_proj[r["BioProject"]].append(r)
    # the setting the rust supplements already gave each project, for siblings to inherit
    resolved = {bp: collections.Counter(r["setting"] for r in rs if r.get("setting"))
                for bp, rs in by_proj.items()}

    out = {}
    for r in runs:
        # already resolved by curated or supplement? leave it to that row
        if r.get("setting_src") in ("curated", "supplement"):
            continue
        bp = r["BioProject"]
        spec = SPECS.get(bp)
        row = {"run": r["run"], "host": "", "tissue": "", "setting": "",
               "location": "", "collection_date": "", "source": ""}
        if spec:
            setting = spec["setting"]
            if spec.get("date_split"):
                cut, ge, lt = spec["date_split"]
                y = year(r.get("collection_date"))
                setting = ge if (y and y >= cut) else lt
            row["setting"] = setting
            if spec.get("host") and not r.get("host"):
                row["host"] = spec["host"]
            if spec.get("location") and not r.get("location"):
                row["location"] = spec["location"]
            row["source"] = f"paper:{spec['doi']}"
        elif resolved.get(bp):
            # a run the supplements missed: inherit the project's dominant resolved setting
            row["setting"] = resolved[bp].most_common(1)[0][0]
            row["source"] = f"sibling:{bp}"
        else:
            continue  # no basis yet (a project still to be read)
        out[r["run"]] = row

    kept = [r for r in (read(OUT) if OUT.exists() else [])
            if not r.get("source", "").startswith(("paper:", "sibling:"))]
    rows = kept + list(out.values())
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    tally = collections.Counter(r["source"].split(":")[0] for r in out.values())
    print(f"-> supplements.tsv  {len(rows)} rows total; added {len(out)} "
          f"({dict(tally)}), {len(kept)} supp-rows kept")
    for bp in sorted({r["run"][:0] or SPECS and bp for bp in SPECS}):  # noqa
        pass
    print("   settings written:",
          dict(collections.Counter(r["setting"] for r in out.values())))
    missing = [bp for bp, rs in by_proj.items()
               if any(r.get("setting_src") not in ("curated", "supplement")
                      and r["run"] not in out for r in rs)]
    if missing:
        print("   still no basis (unread projects):", missing)


if __name__ == "__main__":
    main()
