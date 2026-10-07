#!/usr/bin/env python3
"""Everything known locally about one BioProject, in a single call.

Written 2026-10-06. Curating a project was opening with five or six exploratory calls -
the bioprojects row, a grep across the tables for candidate papers, the BioSample attribute
cache, the run field distributions, the offtarget verdict - and one of those (grepping the
attribute JSON) routinely dumped tens of kilobytes of raw records to find a single key.

This prints the same information as value-counts, capped. The rule it exists to enforce:
**an empty `primary_paper` does not mean no paper is known.** `accession_papers.tsv` is
populated by an earlier accession sweep and already held two open-access papers for
PRJNA1119650 while its primary_paper sat blank.

Run:  python 01a_Literature/project_brief.py PRJNA1119650
"""

from __future__ import annotations

import collections, csv, json, sys
from pathlib import Path

from _layout import GOLD, SEEDS

HERE = Path(__file__).resolve().parent
D = HERE / "data"
CAP = 6          # distinct values printed before collapsing to a count plus examples


def read(name: str) -> list[dict]:
    """Read a table by bare filename, from whichever layer holds it (see _paths)."""
    p = next((d / name for d in (D, SEEDS, GOLD) if (d / name).exists()), D / name)
    if not p.exists():
        return []
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def counts(vals, cap: int = CAP) -> str:
    c = collections.Counter("" if v is None else v for v in vals)
    if len(c) == 1:
        k = next(iter(c))
        return f"{k!r} (all)" if k else "(all blank)"
    if len(c) <= cap:
        return "  ".join(f"{k or '(blank)'}:{n}" for k, n in c.most_common())
    top = "  ".join(f"{k or '(blank)'}:{n}" for k, n in c.most_common(3))
    return f"{len(c)} distinct, top: {top}"


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: project_brief.py <BioProject>")
    bp = sys.argv[1].strip()

    row = next((r for r in read("bioprojects.tsv") if r["BioProject"] == bp), None)
    if not row:
        sys.exit(f"{bp} is not in bioprojects.tsv")
    print(f"=== {bp}  {row.get('organism','')}  [{row.get('centre','')}]")
    print(f"    runs {row.get('sra_runs')}  rnaseq {row.get('sra_rnaseq')}  "
          f"triage {row.get('triage')}  geo% {row.get('geo_pct')}  date% {row.get('date_pct')}")
    print(f"    primary_paper {row.get('primary_paper') or '(unset)'}  "
          f"[{row.get('primary_source')}]  claims {row.get('n_claims')}")

    ot = next((r for r in read("offtarget.tsv") if r["BioProject"] == bp), None)
    if ot:
        print(f"    screen: {ot['recommend'].upper()} - {ot['reason'][:150]}")

    # Candidate papers. This block is the reason the script exists.
    cands = [r for r in read("accession_papers.tsv") if r["accession"] == bp]
    print(f"\n--- candidate papers from the accession sweep: {len(cands)}")
    for r in cands:
        slug = "doi_" + (r["doi"] or "").replace("/", "_")
        d = HERE / "studies" / slug
        filed = "filed" if d.is_dir() and any(d.rglob("*")) else "NOT FILED"
        print(f"    {r['doi']}  PMID {r.get('pmid','')}  OA={r.get('is_oa','')}  {filed}")
        print(f"      {(r.get('title') or '')[:150]}")

    claims = [r for r in read("paper_bioproject.tsv") if r["BioProject"] == bp]
    if claims:
        print(f"\n--- existing claims: {len(claims)}")
        for r in claims:
            print(f"    {r['relation']:<16} {r['paper_ref']}  ({r.get('run_scope') or 'all'})")

    rules = [r for r in read("curation.tsv") if r["BioProject"] == bp]
    joins = [r for r in read("joins.tsv") if r["BioProject"] == bp]
    if rules or joins:
        print(f"\n--- ALREADY CURATED: {len(rules)} rule(s), {len(joins)} join(s)")
        for r in rules:
            print(f"    rule {r['rule_id']:<28} {r['field']} = {r['value'][:60]}")
        for r in joins:
            print(f"    join {r['rule_id']:<28} {r['field']} <- {r['file']}:{r['value_col']}")

    runs = [r for r in read("runs.tsv") if r["BioProject"] == bp]
    if runs:
        print(f"\n--- runs.tsv: {len(runs)}")
        for c in ("LibraryStrategy", "LibrarySelection", "SampleName", "LibraryName",
                  "tissue", "location", "collection_date", "setting", "sampling_selection"):
            if c in runs[0]:
                print(f"    {c:<20} {counts(r[c] for r in runs)}")

    p = D / "biosample_attrs" / f"{bp}.json"
    if not p.exists():
        print("\n--- no BioSample attribute cache")
    else:
        a = json.load(p.open())
        print(f"\n--- BioSample attrs: {a.get('sampled')}/{a.get('total')} "
              f"complete={a.get('complete')}")
        attrs = a.get("attrs", {})
        keys = sorted({k for v in attrs.values() for k in v})
        for k in keys:
            print(f"    {k:<24} {counts(v.get(k, '') for v in attrs.values())}")
        # A design code hiding in an attribute is the cheapest possible join key.
        for k in keys:
            vals = [v.get(k, "") for v in attrs.values()]
            if len(set(vals)) == len(vals) and all("_" in x or "-" in x for x in vals if x):
                print(f"    ^ {k} is UNIQUE per sample and looks like a design code: "
                      f"a zero-regex join key (e.g. {vals[0]!r})")


if __name__ == "__main__":
    main()
