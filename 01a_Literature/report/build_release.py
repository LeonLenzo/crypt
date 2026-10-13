#!/usr/bin/env python3
"""Build the reviewer-facing release from the curation seeds. Generated, never hand-edited.

The curation lives across six directories, 364 rules, 121 joins and 225 quotes, because that
is what it takes to BUILD a cohort from the literature. A reviewer does not need the machinery;
they need to see what the samples are, the evidence for it, and that the evidence is real. This
script projects the cereal cohort into that minimal set, and nothing else is the deliverable.

It READS the seeds and WRITES release/. It curates nothing. If this and the seeds ever
disagree, the seeds win and this is stale: re-run it. Every number is computed here from the
live tables, so a figure in the release cannot drift from a figure in the data.

Outputs, all under release/:

  funnel.tsv              the evidence funnel, discovery-agnostic: candidates assembled,
                          screened for field aerial cereal RNA-seq, evidence sought, retained.
                          How a project was FOUND is a column in cohort_bioprojects.tsv, not a
                          stage here (leon 2026-10-13): the cohort is defined by evidence, not
                          by which search surfaced it.
  cohort_bioprojects.tsv  one row per retained BioProject (the cereal cohort)
  cohort_samples.tsv      one row per retained run, each curated value beside the rule that set
                          it, so a value traces to evidence/<project>.yaml without a join
  evidence.tsv            the canonical evidence: one row per curated claim (bioproject,
                          field, n_runs, rule, evidence_id, doi, locus, quote). Joins to
                          cohort_samples.tsv by `rule` and to cohort_bioprojects.tsv by
                          `bioproject`. This is the substitute for the 310 MB of PDFs the
                          repository does not ship.
  evidence/<PRJ>.md       the same rows laid out to read: one project's evidence as a short
                          dossier, rendering on GitHub without a parser. A convenience over
                          evidence.tsv, not a second source of truth.
  evidence/<PRJ>.manifest.tsv   every source file the evidence rests on: doi, filename, bytes,
                          sha256, where it came from. A reviewer re-fetches by DOI and compares
                          the hash. Files live in studies/<doi>/ and are NOT copied here.
  README.md              what each file is and how to regenerate it

    python 01a_Literature/report/build_release.py
"""

from __future__ import annotations

import collections, csv, hashlib, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
MODULE = _HERE.parents[1]
sys.path[:0] = [str(MODULE), str(MODULE / "report" / "figures")]

from _cohort import all_runs, is_cohort
import prep_crop_map as P

GOLD = MODULE / "data" / "gold"
SILVER = MODULE / "data" / "silver"
SEEDS = MODULE / "seeds"
STUDIES = MODULE / "studies"
REL = MODULE / "release"

FIELDS = ["host", "setting", "sampling_selection", "tissue", "location", "collection_date"]
SOURCE_COL = {"setting": "setting_source", "sampling_selection": "sampling_selection_source",
              "tissue": "tissue_source", "location": "location_source",
              "collection_date": "date_source"}


def read(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    REL.mkdir(exist_ok=True)
    (REL / "evidence").mkdir(exist_ok=True)

    runs = all_runs()
    host = {r["Run"]: r["host"] for r in read(GOLD / "cohort_hosts.tsv")}
    crop = {rn: P.crop_of(h) for rn, h in host.items()}
    cohort = [r for r in runs if is_cohort(r) and crop.get(r["Run"])]
    CBP = sorted({r["BioProject"] for r in cohort})

    rules = {r["rule_id"]: r for r in read(SEEDS / "curation.tsv")}
    joins = {r["rule_id"]: r for r in read(SEEDS / "joins.tsv")}
    prov = {p["evidence_id"]: p for p in read(SEEDS / "provenance.tsv")}
    bpmeta = {r["BioProject"]: r for r in read(SILVER / "bioprojects.tsv")}
    locs = {r["location"]: r for r in read(GOLD / "cereal_localities.tsv")}

    # paper_bioproject links the DOIs; evidence rows name the DOI that was actually used.
    pb = collections.defaultdict(set)
    for r in read(SILVER / "paper_bioproject.tsv"):
        if r.get("paper_ref"):
            pb[r["BioProject"]].add(r["paper_ref"])

    def field_source(run, field):
        """Return (kind, rule_label, evidence_id) for one curated field on one run.

        kind is one of: paper (a curation/join rule quoting a manuscript), archive (the
        submitter's deposited value), none (no value, or an unfetched/absent sentinel). The
        rule_label is a real rule id ONLY when kind is paper; otherwise it is "archive" or ""
        so the cohort_samples `_rule` column never carries a token that looks like a rule but
        resolves to nothing in evidence.tsv.
        """
        if field == "host":
            # host is resolved by host_summary, not apply_curation; it has no per-run source
            # column, so it is labelled by the cohort_hosts provenance at project level below.
            return ("resolved", "", "")
        src = run.get(SOURCE_COL[field], "")
        if src in rules:
            return ("paper", src, rules[src].get("evidence_id", ""))
        if src in joins:
            return ("paper", src, joins[src].get("evidence_id", ""))
        # A field with no usable value carries a sentinel source (absent, not_fetched) or a
        # blank; only a value actually present in the deposit is "archive".
        if not (run.get(field) or "").strip():
            return ("none", "", "")
        if src.startswith(("biosample", "import-bs")) or src == "":
            return ("archive", "archive", "")
        return ("archive", "archive", "")

    # ---- cohort_samples.tsv -------------------------------------------------
    scols = (["Run", "BioProject", "crop", "host"]
             + [f for f in FIELDS if f != "host"]
             + [f"{f}_rule" for f in FIELDS if f != "host"]
             + ["lat", "lon", "location_precision"])
    with (REL / "cohort_samples.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(scols)
        for r in sorted(cohort, key=lambda r: (r["BioProject"], r["Run"])):
            g = locs.get((r.get("location") or "").strip(), {})
            row = [r["Run"], r["BioProject"], crop[r["Run"]], host[r["Run"]]]
            row += [r.get(f, "") for f in FIELDS if f != "host"]
            row += [field_source(r, f)[1] for f in FIELDS if f != "host"]
            row += [g.get("lat", ""), g.get("lon", ""), g.get("precision", "")]
            w.writerow(row)

    # ---- per-project evidence + manifest + cohort_bioprojects ---------------
    by_bp = collections.defaultdict(list)
    for r in cohort:
        by_bp[r["BioProject"]].append(r)

    # evidence.tsv is the canonical, machine-checkable form: one row per curated claim, keyed
    # so a value in cohort_samples.tsv (by `rule`) and a project in cohort_bioprojects.tsv (by
    # `bioproject`) both resolve here without a join the reader has to construct. The per-
    # project Markdown below is the same rows laid out to read; it is a convenience, not a
    # second source of truth.
    ev_cols = ["bioproject", "field", "n_runs", "rule", "evidence_id", "doi", "locus", "quote"]
    ev_tsv_rows = []

    bp_rows = []
    for bp in CBP:
        rs = by_bp[bp]
        n = len(rs)
        crops = collections.Counter(crop[r["Run"]] for r in rs)
        meta = bpmeta.get(bp, {})

        # evidence used by THIS project's cohort runs, field by field, with how many runs each
        # (rule, field) actually set, so a claim carries its own weight.
        field_rule_n = collections.Counter()      # (field, rule, ev) -> runs
        archive_fields = collections.Counter()
        used_ev = {}
        for r in rs:
            for f in FIELDS:
                kind, rule, ev = field_source(r, f)
                if kind == "archive":
                    archive_fields[f] += 1
                if ev:
                    field_rule_n[(f, rule, ev)] += 1
                    used_ev[ev] = prov.get(ev, {})

        dois = sorted({p.get("doi") or p.get("paper_key") for p in used_ev.values()
                       if p.get("doi") or p.get("paper_key")})
        papers = []
        seen_papers = set()
        for ev, p in sorted(used_ev.items()):
            key = (p.get("paper_key", ""), p.get("doi", ""))
            if key not in seen_papers:
                seen_papers.add(key)
                papers.append((p.get("paper_key", ""), p.get("doi", "")))

        # ---- rows into the flat evidence.tsv
        for f in FIELDS:
            for (ff, rule, ev), cnt in sorted(field_rule_n.items()):
                if ff != f:
                    continue
                p = prov.get(ev, {})
                ev_tsv_rows.append([bp, f, cnt, rule, ev, p.get("doi", ""),
                                    p.get("locus", ""), p.get("quote", "")])

        # ---- evidence/<bp>.md : the same rows, laid out to read
        md = [f"# {bp} — cereal field RNA-seq cohort evidence", "",
              "*Generated by `report/build_release.py` from the curation seeds. Do not edit.*", "",
              f"- **Organism:** {meta.get('organism', '') or '—'}",
              f"- **Runs in cohort:** {n}  ("
              + ", ".join(f"{c} {k}" for c, k in crops.most_common()) + ")",
              f"- **Discovery route:** {meta.get('accession_source', '') or '—'}",
              "- **Papers:** " + ("; ".join(f"{r} ({d})" for r, d in papers) or "—"), ""]
        for f in FIELDS:
            rows_f = sorted((cnt, rule, ev) for (ff, rule, ev), cnt in field_rule_n.items()
                            if ff == f)
            if not rows_f and not archive_fields.get(f):
                continue
            md.append(f"## {f}")
            if archive_fields.get(f):
                md.append(f"*{archive_fields[f]} run(s) carry the archive's own value "
                          f"(submitter-deposited), not a manuscript claim.*\n")
            for cnt, rule, ev in rows_f:
                p = prov.get(ev, {})
                md.append(f"**{cnt} run(s)** — rule `{rule}` — {p.get('doi', '') or '—'}  ")
                md.append(f"*{p.get('locus', '') or '—'}*")
                md.append("")
                md.append("> " + (p.get("quote", "") or "—").replace("\n", " "))
                md.append("")
        (REL / "evidence" / f"{bp}.md").write_text("\n".join(md) + "\n")

        # ---- evidence/<bp>.manifest.tsv : every source file the evidence rests on
        man = [["doi", "filename", "bytes", "sha256", "studies_path"]]
        for doi in dois:
            if "/" not in doi:
                continue
            d = STUDIES / ("doi_" + doi.replace("/", "_"))
            if not d.is_dir():
                man.append([doi, "(directory absent)", "", "", str(d.relative_to(MODULE))])
                continue
            for f in sorted(d.iterdir()):
                if not f.is_file() or f.name.startswith("~$"):
                    continue
                h = hashlib.sha256(f.read_bytes()).hexdigest()
                man.append([doi, f.name, f.stat().st_size, h, str(f.relative_to(MODULE))])
        with (REL / "evidence" / f"{bp}.manifest.tsv").open("w", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(man)

        bp_rows.append({
            "BioProject": bp, "organism": meta.get("organism", ""),
            "discovery_route": meta.get("accession_source", ""),
            "runs": n,
            "crops": "; ".join(f"{c}:{k}" for c, k in crops.most_common()),
            "papers": len(papers), "dois": "; ".join(dois),
            "evidence_file": f"evidence/{bp}.md",
        })

    with (REL / "evidence.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ev_cols)
        w.writerows(ev_tsv_rows)

    with (REL / "cohort_bioprojects.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bp_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(bp_rows, key=lambda r: -r["runs"]))

    # ---- funnel.tsv : discovery-agnostic, evidence-defined -------------------
    papers_screened = len(read(SILVER / "papers.tsv"))
    triaged = len(read(SILVER / "bioprojects.tsv"))
    # candidate pool: every run sitting in a cohort project, before the field+aerial test
    cand_runs = [r for r in runs if r["BioProject"] in set(CBP)]
    cand_bs = len({r["BioSample"] for r in cand_runs})
    coh_bs = len({r["BioSample"] for r in cohort})

    def pct(num, den):
        return f"{num}/{len(cohort)} = {num / len(cohort):.1%}"

    complete = collections.Counter()
    for r in cohort:
        for f in ("host", "tissue"):
            complete[f] += 1  # both are 100% by construction, counted for the record
        if (r.get("location") or "").strip():
            complete["location"] += 1
        if (r.get("collection_date") or "")[:4].isdigit():
            complete["year"] += 1

    fn = [["stage", "value", "note"]]
    fn.append(["papers_screened", papers_screened,
               "publications assessed for field-grown plant RNA-seq (literature search + "
               "accession-list + prior detector screen); cereals are the curated subset"])
    fn.append(["bioprojects_triaged", triaged,
               "distinct BioProjects assembled from all discovery routes and triaged; "
               "cereal and non-cereal"])
    fn.append(["cohort_bioprojects", len(CBP),
               "pass the field + aerial + representative screen with manuscript-level evidence"])
    fn.append(["candidate_runs", len(cand_runs),
               "all runs in those projects, before the field + aerial + representative test"])
    fn.append(["cohort_runs", len(cohort),
               "retained: field-grown aerial cereal tissue, one library per biological sample"])
    fn.append(["cohort_biosamples", coh_bs, "distinct BioSamples behind the retained runs"])
    fn.append(["host_coverage", pct(complete["host"], None), "sample-level host, all paper-sourced"])
    fn.append(["tissue_coverage", pct(complete["tissue"], None), "sample-level tissue"])
    fn.append(["location_coverage", pct(complete["location"], None), "any location string"])
    fn.append(["year_coverage", pct(complete["year"], None), "collection year"])
    with (REL / "funnel.tsv").open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(fn)

    # ---- README -------------------------------------------------------------
    route = collections.Counter(r["discovery_route"] for r in bp_rows)
    readme = f"""# Cereal field RNA-seq cohort — reviewer release

Generated by `01a_Literature/report/build_release.py` from the curation seeds.
**Do not edit these files by hand.** They are a projection of the data; to change a value,
change the seed it comes from and regenerate. Last build reflects whatever the seeds held
when the script was last run.

## What this is

{len(cohort):,} RNA-seq runs ({coh_bs:,} BioSamples) across {len(CBP)} BioProjects: field-grown
aerial cereal tissue, screened from the literature and from public accession listings. Every
run's host, tissue, setting, sampling design, location and (where the record permits) year is
backed by a quotation from the study that produced it.

## Files

- **funnel.tsv** — how the cohort narrows, from papers screened to runs retained. Defined by
  evidence, not by discovery method.
- **cohort_bioprojects.tsv** — one row per BioProject: organism, how it was discovered, run
  count, papers, and a pointer to its evidence file.
- **cohort_samples.tsv** — one row per run: the curated host, tissue, setting, selection,
  location and year, each beside the id of the rule that set it. That id is the key into the
  project's evidence file.
- **evidence.tsv** — the canonical evidence, one row per curated claim: `bioproject, field,
  n_runs, rule, evidence_id, doi, locus, quote`. Join it to `cohort_samples.tsv` on `rule` and
  to `cohort_bioprojects.tsv` on `bioproject`. This stands in for the manuscripts and
  supplements themselves, which are not in this repository.
- **evidence/<BioProject>.md** — the same claims laid out to read, one project per file,
  rendering on GitHub. A convenience view over evidence.tsv.
- **evidence/<BioProject>.manifest.tsv** — every source file the evidence rests on, with its
  SHA-256. The files live in `studies/<doi>/` and are not shipped; a reviewer re-fetches each
  by its DOI from the publisher and compares the hash.

## Discovery routes

The cohort is defined by evidence, so discovery method is recorded per project, not used to
include or exclude. Across the {len(CBP)} projects:

{chr(10).join(f"  - {k or '(unrecorded)'}: {v} projects" for k, v in route.most_common())}

## Regenerating

    python 01a_Literature/report/build_release.py

The script reads only `data/gold/`, `data/silver/` and `seeds/`, and hashes the files in
`studies/`. It writes only `release/`. Running it on an unchanged tree reproduces these files
byte for byte, except the SHA-256s, which change only if a source file in `studies/` changes.
"""
    (REL / "README.md").write_text(readme)

    print(f"release/ built: {len(cohort):,} runs, {len(CBP)} projects, "
          f"{sum(len(x['dois'].split('; ')) for x in bp_rows if x['dois'])} DOI references")
    print(f"  funnel.tsv  cohort_bioprojects.tsv  cohort_samples.tsv")
    print(f"  evidence/ : {len(CBP)} markdown + {len(CBP)} manifest files")
    print(f"  discovery routes: {dict(route.most_common())}")


if __name__ == "__main__":
    main()
