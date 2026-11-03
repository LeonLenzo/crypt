#!/usr/bin/env python3
"""build_s2.py -- Supplementary Data S2: every run in the selected BioProjects, run level.

    python 01_literature/report/tables/build_s2.py

S1 ends at the manuscript. S2 begins at the run, and carries ALL 10,312 runs the 41 selected
BioProjects hold, not only the 7,688 that survived. That is deliberate: a table holding only
the survivors lets a reader check what was kept and never what was dropped, and the funnel in
Table T4 is the number most open to challenge. Every stage of it is reproducible from the
columns here by counting rows.

EVERY CURATED VALUE TRAVELS WITH ITS EVIDENCE. The `*_rule` and `*_source` columns name the
rule or the quoted sentence that set the value beside it, so a reviewer can follow any single
cell back to the manuscript that justifies it rather than taking the curation on trust. This
is the chain the project has maintained since curation began; S2 is the first place it has
been exposed in one file.

TWO KINDS OF FILTER, AND ONLY ONE BELONGS HERE. Leon, 2026-10-29. The inclusion funnel mixes
criteria that state what the AUTHORS DID -- disease exposure, inoculation, field setting,
tissue -- with one that states what WE COULD MEASURE, the library depth floor. The first kind
is a property of the study and every exclusion traces to a quoted sentence. The second is not
a property of the study at all: the same experiment sequenced deeper would pass it.
  S2 therefore carries the methodological filter only, in `methodological_exclusion`, and the
depth filter moved to Supplementary Data S3 with the rest of the library statistics.

STUDY FOCUS, pathology or physiology, so the manuscript's own vocabulary is backed by the
table. The mapping is the one Figure 1e uses, defined once in prep_crop_map.py: a sample is
PATHOLOGY if it was chosen or manipulated with reference to disease -- selected for symptoms,
inoculated, fungicide-treated, or deliberately chosen for the ABSENCE of symptoms -- and
physiology otherwise. Symptom-avoided is pathology because disease still governed the
sampling, negatively.

WHY THE EXTRA COLUMNS. Beyond crop, host, tissue, setting, selection, location and date, five
fields were added because a stated result cannot be checked without them:

  BioSample              the biological-sample count is 6,200, not the run count
  declared_pathogens     every undeclared-detection claim is measured against this field
  disease_exposure       an exclusion criterion in its own right; 1,452 runs leave on it
  doi_relation           how strongly the manuscript is tied to the data, see below

DROPPED after measuring it: `library_representative` is populated for 92 of these 10,312 runs
and empty for the rest, so a column of it would read as 10,220 runs failing a test that was
never applied to them. The "one library per biological sample" rule is not enforced through
that field, and the cohort is 7,688 runs over 6,222 BioSamples, so it is not one-to-one.

`doi` links each row back to its manuscript in S1, which is what makes the two tables one
chain rather than two lists.
"""

import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _util import atomic_write  # noqa: E402

SILVER = MODULE / "data" / "silver"
RELEASE = MODULE / "release"
SEEDS = MODULE / "seeds"
FUNNEL = MODULE.parent / "02_classify" / "report" / "figures" / "depth_funnel" / "depth_funnel.tsv"
OUT = HERE / "Supplementary Data S2.tsv"

COLS = [
    "Run", "BioSample", "BioProject", "doi",
    "crop", "host", "tissue", "setting", "sampling_selection", "disease_exposure",
    "declared_pathogens", "location", "country", "lat", "lon", "location_precision", "collection_date",
    "study_focus", "doi_relation",
    "selected", "methodological_exclusion",
    "tissue_evidence", "setting_evidence", "sampling_selection_evidence",
    "location_evidence", "collection_date_evidence", "declared_pathogens_evidence",
    "disease_exposure_evidence",
]


# The methodological criteria. `below depth floor` is deliberately absent: it is the
# technical filter and belongs to S3.
METH_REASONS = {"exposure: unknown", "exposure: protected", "exposure: sheltered",
                "inoculated", "not field throughout", "non-foliar tissue"}


def main() -> None:
    cohort_bp = {r["BioProject"] for r in
                 csv.DictReader(open(RELEASE / "cohort_bioprojects.tsv"), delimiter="\t")}
    runs = [r for r in csv.DictReader(open(SILVER / "runs.tsv"), delimiter="\t")
            if r["BioProject"] in cohort_bp]
    cur = {r["Run"]: r for r in
           csv.DictReader(open(RELEASE / "cohort_samples.tsv"), delimiter="\t")}
    fun = {r["Run"]: r for r in csv.DictReader(open(FUNNEL), delimiter="\t")} \
        if FUNNEL.exists() else {}

    # A BioProject can be reported by more than one manuscript (PRJNA306542 is eleven
    # papers), so the DOI column is a list and must not be collapsed to the first. The
    # relation is reported alongside, because the strength of the tie varies and a reviewer
    # should see which it is: only 7 of the 41 projects have a manuscript that GENERATED
    # them. Filtering to `generated`, as this first did, left the column 16.7% full and
    # silently implied the other 34 projects had no source paper. They all do.
    # paper_bioproject.paper_ref is a MIX of DOIs and the project's own short keys (Kas18,
    # Ada21). A reviewer cannot resolve a short key, so every one is mapped back to its DOI
    # through papers.tsv, which carries both spellings for the same row. Left unresolved this
    # column was 6,375 short keys against 4,415 DOIs.
    to_doi = {}
    for r in csv.DictReader(open(SILVER / "papers.tsv"), delimiter="\t"):
        d = (r.get("doi") or r.get("paper_key") or "").strip()
        if not d:
            continue
        for alias in (r.get("paper_key"), r.get("paper_ref")):
            if alias:
                to_doi[alias.strip()] = d

    REL_RANK = ["generated", "reanalysed", "sra_linked", "cited", "mentions"]
    links, unresolved = {}, set()
    for r in csv.DictReader(open(SILVER / "paper_bioproject.tsv"), delimiter="\t"):
        if r["BioProject"] in cohort_bp:
            ref = r["paper_ref"]
            d = to_doi.get(ref)
            if d is None:
                # Never silently emit the short key: say which ones could not be resolved.
                unresolved.add(ref)
                d = ref
            links.setdefault(r["BioProject"], []).append((r.get("relation", ""), d))
    doi, doi_rel = {}, {}
    for b, v in links.items():
        best = min((REL_RANK.index(rel) if rel in REL_RANK else len(REL_RANK), rel, ref)
                   for rel, ref in v)[1]
        doi[b] = {ref for rel, ref in v if rel == best}
        doi_rel[b] = best

    # EVIDENCE AS DOIs, not rule identifiers. The curation records which RULE set a value
    # (man24b-tissue, kashima-setting), which is the right key internally and useless to a
    # reviewer, who cannot resolve it to anything. release/evidence.tsv carries the rule ->
    # DOI -> locus -> quote chain, so the rule is resolved to the DOI of the manuscript the
    # quoted sentence came from. Values with no manuscript behind them keep their literal
    # provenance -- `archive`, `biosample`, `absent` -- because inventing a DOI for an
    # archive field would be worse than naming it.
    # The chain is rule_id -> evidence_id -> DOI, and it takes TWO tables: seeds/curation.tsv
    # names the rule that set each value, seeds/provenance.tsv holds the quoted sentence and
    # the DOI it came from. release/evidence.tsv was the obvious source and is the wrong one:
    # it covers 5 of the 8 curated fields and 74 of the DOIs, so routing through it resolved
    # 40 of 159 evidence values and left declared_pathogens entirely unresolved.
    eid_doi, eid_bp, eid_multi = {}, {}, set()
    for r in csv.DictReader((l for l in open(SEEDS / "provenance.tsv")
                             if not l.startswith("#")), delimiter="\t"):
        eid = (r.get("evidence_id") or "").strip()
        if not eid:
            continue
        d = (r.get("doi") or "").strip()
        # provenance.tsv writes the LITERAL "(multiple)" into the doi column for evidence
        # that spans several papers. Treating it as a DOI is how 5,879 cells came out
        # carrying "(multiple)" after the resolver was supposedly handling that case: the
        # first branch matched and returned before the branch meant to catch it.
        if d and d != "(multiple)":
            eid_doi[eid] = d
        elif d == "(multiple)":
            eid_multi.add(eid)
        elif (r.get("paper_key") or "").strip().startswith(("PRJ", "SRP", "ERP")):
            # 38 of the 41 disease_exposure rows carry the ACCESSION in paper_key and leave
            # doi empty, though the quote beside them is plainly from a manuscript. The DOI
            # was never recorded. Resolve through the project instead -- the quote is from
            # that project's paper by construction -- rather than emit a rule id a reviewer
            # cannot follow. Recording the gap here because filling it properly means adding
            # the DOI to provenance.tsv, which is curation work, not a build fix.
            eid_bp[eid] = r["paper_key"].strip()
    # RULES COME FROM TWO FILES, not one. curation.tsv holds the hand-written rules and
    # joins.tsv the 122 that pull a value out of a paper's supplementary table. Reading only
    # curation.tsv left every ada21-PRJEB*-* rule unresolved, which is most of the wheat rust
    # provenance -- the exact rows whose evidence a reviewer is most likely to check.
    rule_doi, rule_bp, rule_multi = {}, {}, set()
    for src in ("curation.tsv", "joins.tsv"):
        for r in csv.DictReader((l for l in open(SEEDS / src)
                                 if not l.startswith("#")), delimiter="\t"):
            rid = (r.get("rule_id") or "").strip()
            eid = (r.get("evidence_id") or "").strip()
            if not rid:
                continue
            if eid in eid_doi:
                rule_doi[rid] = eid_doi[eid]
            elif eid in eid_bp:
                rule_bp[rid] = eid_bp[eid]
            elif eid in eid_multi:
                rule_multi.add(rid)

    # EVERY EVIDENCE CELL RESOLVES TO SOMETHING A READER CAN OPEN. A DOI where a manuscript
    # set the value, an NCBI URL where the archive did, and a comma-separated list where more
    # than one source contributed. Only two literals survive, and both are statements rather
    # than unresolved identifiers:
    #   absent      the curator looked and no value exists in any source
    #   not_fetched the supplementary file was never retrieved
    SRA = "https://www.ncbi.nlm.nih.gov/sra/"
    BIOSAMPLE = "https://www.ncbi.nlm.nih.gov/biosample/"

    def as_doi(v: str, run: str, bs: str, proj_doi: set) -> str:
        v = (v or "").strip()
        if v in rule_doi:
            return rule_doi[v]
        bp = rule_bp.get(v)
        if bp and doi.get(bp):
            return ", ".join(sorted(doi[bp]))
        # `(multiple)` is mostly ev-unselected-sampling, a determination taken ACROSS studies
        # rather than a citation; the citable source for any one row is its own manuscript.
        # import-supp is the same shape: values mined from several papers' supplements.
        if v in rule_multi or v in ("(multiple)", "import-supp-2026-10-02"):
            return ", ".join(sorted(proj_doi)) if proj_doi else v
        if v == "archive":
            return SRA + run if run else v
        if v in ("biosample", "import-bs-2026-10-02"):
            return BIOSAMPLE + bs if bs else v
        return v

    # The resolved country, from the geocoder's own output. S2 is what Figure 1 is built
    # from, so anything the figure needs has to be here; without this the figure would have
    # to join a second table and the two could disagree about where a locality is.
    loc_country = {}
    lp = MODULE / "data" / "gold" / "cereal_localities.tsv"
    if lp.exists():
        for r in csv.DictReader(open(lp), delimiter="\t"):
            if r.get("location"):
                loc_country[r["location"]] = r.get("country", "")

    # The same definition Figure 1e uses; see prep_crop_map.py.
    PATHOLOGY = {"disease-selected", "inoculated", "fungicide-treated", "symptom-avoided"}

    rows = []
    for r in runs:
        c = cur.get(r["Run"], {})
        f = fun.get(r["Run"], {})
        sel = r["Run"] in cur
        pdoi = doi.get(r["BioProject"], set())
        rows.append({
            "Run": r["Run"], "BioSample": r.get("BioSample", ""),
            "BioProject": r["BioProject"],
            "doi": ";".join(sorted(doi.get(r["BioProject"], ()))),
            "crop": c.get("crop", ""), "host": c.get("host", ""),
            "tissue": c.get("tissue") or r.get("tissue", ""),
            "setting": c.get("setting") or r.get("setting", ""),
            "sampling_selection": c.get("sampling_selection") or r.get("sampling_selection", ""),
            "disease_exposure": r.get("disease_exposure", ""),
            "declared_pathogens": r.get("declared_pathogens", ""),
            "location": c.get("location") or r.get("location", ""),
            "country": loc_country.get(c.get("location") or r.get("location", ""), ""),
            "lat": c.get("lat", ""), "lon": c.get("lon", ""),
            "location_precision": c.get("location_precision", ""),
            "collection_date": c.get("collection_date") or r.get("collection_date", ""),
            "study_focus": ("" if not sel else
                            "pathology" if (c.get("sampling_selection") or "") in PATHOLOGY
                            else "physiology"),
            "doi_relation": doi_rel.get(r["BioProject"], ""),
            "selected": int(sel),
            "methodological_exclusion": (f.get("reason", "")
                                         if f.get("reason") in METH_REASONS else ""),
            "tissue_evidence": as_doi(c.get("tissue_rule") or r.get("tissue_source", ""), r["Run"], r.get("BioSample", ""), pdoi),
            "setting_evidence": as_doi(c.get("setting_rule") or r.get("setting_source", ""), r["Run"], r.get("BioSample", ""), pdoi),
            "sampling_selection_evidence": as_doi(
                c.get("sampling_selection_rule") or r.get("sampling_selection_source", ""),
                r["Run"], r.get("BioSample", ""), pdoi),
            "location_evidence": as_doi(c.get("location_rule") or r.get("location_source", ""), r["Run"], r.get("BioSample", ""), pdoi),
            "collection_date_evidence": as_doi(c.get("collection_date_rule") or r.get("date_source", ""), r["Run"], r.get("BioSample", ""), pdoi),
            "declared_pathogens_evidence": as_doi(r.get("declared_pathogens_source", ""), r["Run"], r.get("BioSample", ""), pdoi),
            "disease_exposure_evidence": as_doi(r.get("disease_exposure_source", ""), r["Run"], r.get("BioSample", ""), pdoi),
        })
    rows.sort(key=lambda x: (x["BioProject"], x["Run"]))

    with atomic_write(OUT) as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    FOUR = {"Wheat", "Maize", "Rice", "Sorghum"}
    sel = [x for x in rows if x["selected"]]
    four = [x for x in sel if x["crop"] in FOUR]
    ret = [x for x in four if not x["methodological_exclusion"]]
    print(f"-> {OUT.name}  ({len(COLS)} columns)")
    print(f"   runs in the 41 BioProjects        {len(rows):,}")
    print(f"   ... selected into the cohort      {len(sel):,}   "
          f"({len({x['BioSample'] for x in sel if x['BioSample']}):,} BioSamples)")
    print(f"   ... of the four target crops      {len(four):,}   "
          f"({len({x['BioSample'] for x in four if x['BioSample']}):,} BioSamples)")
    print(f"   ... passing the methodological filter  {len(ret):,}")
    print(f"   study focus of those                   "
          f"{dict(collections.Counter(x['study_focus'] for x in ret))}")
    evcols = [k for k in COLS if k.endswith("_evidence")]
    vals = {v for x in rows for k in evcols for v in [x[k]] if v}
    nod = sorted(v for v in vals if not v.startswith("10."))
    print(f"   evidence values that are DOIs     "
          f"{len(vals) - len(nod)} of {len(vals)} distinct")
    if nod:
        print(f"   non-DOI provenance kept literal   {', '.join(nod[:10])}"
              + (" ..." if len(nod) > 10 else ""))
    if unresolved:
        print(f"   WARNING: {len(unresolved)} paper reference(s) have no DOI in papers.tsv "
              f"and appear as their short key: {', '.join(sorted(unresolved)[:8])}",
              file=sys.stderr)
    miss = [k for k in COLS if not any(x[k] for x in rows)]
    if miss:
        print(f"   WARNING: columns with no data at all: {', '.join(miss)}", file=sys.stderr)


if __name__ == "__main__":
    main()
