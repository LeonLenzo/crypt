#!/usr/bin/env python3
"""Freeze the research tree's curation into release/01_cohort/inputs/.

    python 01_literature/report/freeze_cohort_inputs.py

This is a RESEARCH-TREE tool. It may import research code; the files it writes may not
depend on any. Its job is to collapse the curation machinery into flat tables so the
release needs none of it:

  - the rule -> evidence_id -> provenance chain resolves to doi/locus/quote on one row
  - the 122 `joins.tsv` rules, which read values out of author spreadsheets, are
    pre-resolved to literal per-run values, so no release script parses xlsx
  - the five BioSample attributes that predicates match on become columns of runs.tsv,
    retiring a 367-file 5.3 MB JSON cache
  - everything keeps the evidence that backs it

Scraped tables go to inputs/ as `raw_*`; tables carrying a human decision go to curation/ as
`curated_*`. The prefix is deliberate even though the directory already says it: these files
get cited and opened on their own, and a reader should never have to recover from context
whether a value was retrieved or decided.

Run it whenever the curation changes. It overwrites; it never merges.
"""

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
ROOT = MODULE.parent
sys.path.insert(0, str(MODULE))

from _resolver import attrs_for  # noqa: E402

SEEDS = MODULE / "seeds"
SILVER = MODULE / "data" / "silver"
REPO = ROOT.parent          # research/ and release/ are siblings
MOD_OUT = REPO / "release" / "01_cohort"
RAW = MOD_OUT / "inputs"          # scraped from the web. No human judgement here.
CUR = MOD_OUT / "curation"        # the module's work: decisions made by hand.
DEPTH = ROOT / "02_classify" / "data" / "silver" / "run_depth.tsv"

# The attribute keys the curation predicates actually match on. Only these travel; the
# rest of the BioSample record does not.
# A value with no origin recorded, or one that was never retrieved. Carries nothing.
NO_VALUE = ("", "absent", "not_fetched")

ATTR_KEYS = ("dev_stage", "treatment", "bs_description", "isolate", "cultivar")

CURATED_FIELDS = ("tissue", "location", "collection_date", "setting",
                  "sampling_selection", "declared_pathogens", "disease_exposure",
                  "tissue_selection")


def read(path, **kw):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", **kw))


def write(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  {path.name:<22} {len(rows):>7,} rows  {len(cols)} cols")


def evidence_index():
    """evidence_id -> (doi, locus, quote), from provenance.tsv."""
    out = {}
    for r in read(SEEDS / "provenance.tsv"):
        eid = (r.get("evidence_id") or "").strip()
        if eid:
            out[eid] = (r.get("doi", ""), r.get("locus", ""), r.get("quote", ""))
    return out


def freeze_curation(ev):
    """Every curation rule, from both rule tables, with its evidence inline.

    Two kinds, distinguished by `scope`:

      project   the rule states one value for every run its predicate matches, and `value`
                holds it
      per_run   the value was read out of an author supplementary table and differs run to
                run, so `value` is empty and the values are in curation_runs.tsv

    Merging the two tables is what keeps the evidence here and only here. Carrying the quote
    on each per-run value instead repeated 46 quotes about 207 times each and cost 7.8 MB.
    """
    rows = []
    for r in read(SEEDS / "curation.tsv"):
        if not (r.get("rule_id") or "").strip():
            continue
        doi, locus, quote = ev.get((r.get("evidence_id") or "").strip(), ("", "", ""))
        rows.append({
            "rule": r["rule_id"], "BioProject": r["BioProject"],
            "predicate": r.get("predicate", ""), "field": r["field"], "scope": "project",
            "value": r.get("value", ""), "override": r.get("override", ""),
            "doi": doi, "locus": locus, "quote": quote, "note": r.get("note", ""),
        })
    for j in read(SEEDS / "joins.tsv"):
        rid = (j.get("rule_id") or "").strip()
        if not rid:
            continue
        doi, locus, quote = ev.get((j.get("evidence_id") or "").strip(), ("", "", ""))
        where = " ".join(x for x in (j.get("file", ""), j.get("sheet", ""),
                                     j.get("value_col", "")) if x)
        rows.append({
            "rule": rid, "BioProject": j["BioProject"],
            "predicate": j.get("predicate", ""), "field": j["field"], "scope": "per_run",
            "value": "", "override": j.get("override", ""),
            "doi": doi, "locus": locus or where, "quote": quote, "note": j.get("note", ""),
        })
    write(CUR / "curated_rules.tsv", rows,
          ["rule", "BioProject", "predicate", "field", "scope", "value", "override",
           "doi", "locus", "quote", "note"])
    return ({r["rule"] for r in rows},
            {r["rule"] for r in rows if r["scope"] == "per_run"})


def freeze_curation_runs(per_run_rules, keep):
    """The per-run values the supplementary-table rules produced.

    Read from silver/runs.tsv, where `<field>_source` names the rule that set each value,
    NOT from cohort_samples.tsv. cohort_samples holds only the runs selected into the
    cohort, so sourcing it dropped per-run values for every run in a cohort BioProject that
    the selection later excluded, and those runs still appear in S2.

    Evidence is not repeated here; join `rule` to curation.tsv.
    """
    rows = []
    for r in read(SILVER / "runs.tsv"):
        if r["BioProject"] not in keep:
            continue
        for field in CURATED_FIELDS:
            src = (r.get(f"{field}_source") or "").strip()
            if field == "collection_date" and not src:
                src = (r.get("date_source") or "").strip()
            if src in per_run_rules:
                rows.append({"Run": r["Run"], "field": field,
                             "value": r.get(field, ""), "rule": src})
    write(CUR / "curated_rules_runs.tsv", rows, ["Run", "field", "value", "rule"])


def freeze_bioprojects():
    """Every BioProject the screen examined, and how far each got.

    This is what narrows 421 examined accessions to the 41 that carry the cohort. Kept
    generous because it is 421 rows: the columns record the screening decision, its route,
    and the umbrella-accession problem (`n_claims`, `primary_conflict`) where one accession
    covers several papers.
    """
    cols = ["BioProject", "organism", "centre", "triage", "accession_source",
            "undermind_decision", "primary_paper", "n_claims", "primary_conflict",
            "sra_runs", "sra_rnaseq", "sra_biosamples", "screened", "gate_passed",
            "in_cohort"]
    write(RAW / "raw_bioprojects.tsv", read(SILVER / "bioprojects.tsv"), cols)


def cohort_projects():
    return {r["BioProject"] for r in read(MODULE / "release" / "cohort_bioprojects.tsv")}


def freeze_runs(keep, project_rules, per_run_rules):
    """One row per archived run in a cohort BioProject: archive-side values only.

    SCOPED to the 41 cohort accessions. silver/runs.tsv holds 33,958 runs across 417
    screened projects, but the 417 -> 41 narrowing is a project-level story and
    bioprojects.tsv tells it. Carrying runs from rejected projects would treble the release
    and answer a question no table here asks.

    Every curated value has an origin recorded in `<field>_source`, and this classifies
    them EXHAUSTIVELY rather than listing the archive ones, so a new origin cannot be
    dropped in silence. Values owned by a rule live in curation.tsv or curation_runs.tsv;
    everything else came from the archive side and travels here with its origin named. The
    counts are printed because an earlier version listed only `archive` and silently lost
    4,055 tissues, 2,968 locations and 1,224 dates to the `biosample` and `import-bs`
    origins.
    """
    runs = [r for r in read(SILVER / "runs.tsv") if r["BioProject"] in keep]
    cache, rows = {}, []
    tally = {"rule": 0, "archive": 0, "none": 0}
    for r in runs:
        bp = r["BioProject"]
        if bp not in cache:
            cache[bp] = attrs_for(bp)
        a = cache[bp].get(r.get("BioSample", ""), {}) or {}
        out = {k: r.get(k, "") for k in
               ("Run", "BioProject", "BioSample", "Experiment", "LibraryName",
                "SampleName", "LibraryStrategy", "LibrarySelection", "spots")}
        for k in ATTR_KEYS:
            v = a.get(k)
            out[k] = v.strip() if isinstance(v, str) else (v or "")
        for f in CURATED_FIELDS:
            src = (r.get(f"{f}_source") or "").strip()
            if f == "collection_date" and not src:
                src = (r.get("date_source") or "").strip()
            if src in project_rules or src in per_run_rules:
                out[f], out[f"{f}_origin"] = "", ""
                tally["rule"] += 1
            elif src in NO_VALUE:
                out[f], out[f"{f}_origin"] = "", ""
                tally["none"] += 1
            else:
                out[f], out[f"{f}_origin"] = r.get(f, ""), src
                tally["archive"] += 1
        rows.append(out)
    cols = (["Run", "BioProject", "BioSample", "Experiment", "LibraryName", "SampleName",
             "LibraryStrategy", "LibrarySelection", "spots"] + list(ATTR_KEYS)
            + [c for f in CURATED_FIELDS for c in (f, f"{f}_origin")])
    write(RAW / "raw_runs.tsv", rows, cols)
    print(f"    value origins: {tally['rule']:,} by rule, "
          f"{tally['archive']:,} archive-side, {tally['none']:,} none")


def freeze_papers():
    """Every screened manuscript and the accessions it reports.

    The host filter that narrows these to candidates is NOT applied here. It is an
    analysis step and belongs in the release, where it can be read.
    """
    pb = read(SILVER / "paper_bioproject.tsv")
    by_ref = {}
    for r in pb:
        by_ref.setdefault((r.get("paper_ref") or "").strip(), []).append(r["BioProject"])
    rows = []
    for p in read(SILVER / "papers.tsv"):
        bps = sorted({b for k in {(p.get("paper_ref") or "").strip(),
                                  (p.get("doi") or "").strip(),
                                  (p.get("paper_key") or "").strip()} if k
                      for b in by_ref.get(k, [])})
        rows.append({
            "paper_key": p.get("paper_key", ""), "doi": p.get("doi", ""),
            "title": " ".join((p.get("title") or "").split()),
            "sources": p.get("sources", ""), "status": p.get("status", ""),
            "n_bioprojects": len(bps), "bioprojects": ";".join(bps),
        })
    write(RAW / "raw_papers.tsv", rows,
          ["paper_key", "doi", "title", "sources", "status",
           "n_bioprojects", "bioprojects"])


def freeze_library_stats():
    """ENA library statistics. Drops the fastq URL list, which is large and transient."""
    if not DEPTH.exists():
        print(f"  SKIPPED library_stats.tsv: no {DEPTH}")
        return
    cols = ["Run", "BioProject", "layout", "strategy", "selection", "reads", "bases",
            "read_len", "platform", "model", "n_fastq", "fastq_bytes"]
    write(RAW / "raw_library_stats.tsv", read(DEPTH), cols)


def freeze_decisions():
    """The scope and exclusion log: what was ruled out of the cohort, and why."""
    rows = read(SEEDS / "decisions.tsv")
    write(CUR / "curated_decisions.tsv", rows,
          ["ts", "level", "key", "decision", "reason", "evidence", "decided_by"])


def main():
    print(f"-> {MOD_OUT.relative_to(ROOT)}")
    ev = evidence_index()
    all_rules, per_run_rules = freeze_curation(ev)
    keep = cohort_projects()
    freeze_curation_runs(per_run_rules, keep)
    freeze_bioprojects()
    freeze_runs(keep, all_rules - per_run_rules, per_run_rules)
    freeze_papers()
    freeze_library_stats()
    freeze_decisions()


if __name__ == "__main__":
    main()
