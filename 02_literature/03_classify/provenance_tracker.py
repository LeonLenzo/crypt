#!/usr/bin/env python3
"""provenance_tracker.py — which BioProjects still need provenance, and which are done.

A worklist, not an analysis. Ranks every project in the analysed cohort by how much provenance
is still missing, so chasing supplementary tables can be prioritised and crossed off. Reruns
cheaply, so it stays current as supplements land in
02_literature/02_text/data/supp_data/.

"Spatial resolution" is the thing that matters and is easy to mistake for coverage: a project
labelled `Ethiopia` for all 153 samples is nominally resolved while carrying no within-project
structure at all. `localities` counts distinct location strings, so 1 means a single label.

**A single location is often the correct answer, not a gap.** For a controlled inoculation
experiment every sample really did come from one glasshouse, and chasing its supplement cannot
help: PRJNA1217477's supplement turned out to tabulate the inoculum isolates, not the plants.
Projects whose title indicates indoor or in-vitro work are therefore marked SETTING? rather than
TARGET. They are not provenance problems, they are candidates for exclusion from a field cohort,
which is a bigger issue: up to ~13% of the "field" cohort looks like it is not field, and that
propagates into the field-versus-greenhouse comparison the chapter rests on.

    python 02_literature/03_classify/provenance_tracker.py
    python 02_literature/03_classify/provenance_tracker.py --min-samples 10 --md
"""
import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _paths import ROOT, SAMPLES  # noqa: E402

# Title words that place the work indoors or in vitro, so a single location is expected.
NONFIELD = re.compile(r"\b(in vitro|artificial|detached|leaf sheath|axenic|appressorium|"
                      r"growth chamber|greenhouse|glasshouse|compatible and incompatible|"
                      r"wild[- ]type and|mutant|inoculated with|transcriptomic atlas|"
                      r"time series|time ?course)\b", re.I)

PROV = Path(__file__).resolve().parent / "data/cohort_provenance.tsv"
SUPP = ROOT / "02_literature/02_text/data/supp_provenance.tsv"
TEXTCACHE = ROOT / "02_literature/02_text/data/text_cache.jsonl"
OUT = Path(__file__).resolve().parent / "data/provenance_tracker.tsv"


def pmcids():
    """DOI -> PMC id, from the full-text cache's pdf_url field (recorded as epmc:PMCnnn)."""
    out = {}
    if not TEXTCACHE.exists():
        return out
    for line in open(TEXTCACHE):
        k, _, payload = line.partition("\t")
        try:
            d = json.loads(payload)
        except ValueError:
            continue
        u = str(d.get("pdf_url") or "")
        if "epmc:" in u:
            out[k.strip()] = u.split("epmc:")[1].strip()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-samples", type=int, default=10,
                    help="only list projects with at least this many samples (default 10)")
    ap.add_argument("--md", action="store_true", help="markdown table instead of fixed width")
    args = ap.parse_args()

    prov = list(csv.DictReader(open(PROV), delimiter="\t"))
    samples = {r["BioSample"]: r for r in csv.DictReader(open(SAMPLES), delimiter="\t")}
    supp_files = collections.defaultdict(set)
    if SUPP.exists():
        bs_file = {r["BioSample"]: r["source_file"] for r in
                   csv.DictReader(open(SUPP), delimiter="\t")}
    else:
        bs_file = {}
    pmc = pmcids()

    byp = collections.defaultdict(list)
    for r in prov:
        byp[r["BioProject"]].append(r)

    rows = []
    for p, rs in byp.items():
        n = len(rs)
        locs = {r["location"] for r in rs if r["location"]}
        yrs = sum(1 for r in rs if r["collection_year"])
        src = collections.Counter(r["location_source"] for r in rs).most_common(1)[0][0]
        files = {bs_file[r["BioSample"]] for r in rs if r["BioSample"] in bs_file}
        meta = samples.get(rs[0]["BioSample"], {})
        doi = next((r["doi"] for r in rs if r["doi"]), "")
        # the gap: samples lacking a year, plus samples with no spatial structure at all
        gap = (n - yrs) + (n if len(locs) <= 1 else 0)
        title = meta.get("title") or ""
        if files:
            status = "DONE"
        elif len(locs) > 1 and yrs == n:
            status = "ok"
        elif NONFIELD.search(title):
            # one location is expected here; the open question is whether it belongs in a
            # field cohort at all, not whether its geography can be improved
            status = "SETTING?"
        elif not doi:
            status = "NO DOI"
        else:
            status = "TARGET"
        rows.append({
            "status": status, "BioProject": p, "n": n,
            "localities": len(locs), "years_known": yrs,
            "source": src, "gap_score": gap,
            "pmc": pmc.get(doi, ""), "doi": doi,
            "supp_file": ";".join(sorted(files)),
            "title": (meta.get("title") or "")[:60],
        })

    order = {"TARGET": 0, "NO DOI": 1, "SETTING?": 2, "ok": 3, "DONE": 4}
    rows.sort(key=lambda r: (order.get(r["status"], 9), -r["gap_score"]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["status", "BioProject", "n", "localities", "years_known", "source", "gap_score",
            "pmc", "doi", "supp_file", "title"]
    with open(OUT, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    st = collections.Counter(r["status"] for r in rows)
    tot = sum(r["n"] for r in rows)
    print(f"{len(rows)} projects, {tot} BioSamples")
    for k in ("TARGET", "NO DOI", "SETTING?", "ok", "DONE"):
        if st[k]:
            ns = sum(r["n"] for r in rows if r["status"] == k)
            print(f"  {k:7s} {st[k]:4d} projects {ns:5d} samples")
    print(f"\n{OUT}\n")

    show = [r for r in rows if r["n"] >= args.min_samples]
    if args.md:
        print("| | project | n | locs | yrs | PMC | DOI | supp / title |")
        print("|---|---|---|---|---|---|---|---|")
        for r in show:
            box = {"DONE": "[x]", "ok": "[~]", "SETTING?": "[!]"}.get(r["status"], "[ ]")
            last = r["supp_file"] or r["title"]
            print(f"| {box} | {r['BioProject']} | {r['n']} | {r['localities']} | "
                  f"{r['years_known']} | {r['pmc'] or '-'} | `{r['doi'] or '-'}` | {last} |")
    else:
        print(f"{'':3s} {'project':14s} {'n':>4s} {'locs':>4s} {'yrs':>4s} {'PMC':13s} DOI")
        for r in show:
            box = {"DONE": "[x]", "ok": "[~]", "SETTING?": "[!]"}.get(r["status"], "[ ]")
            print(f"{box:3s} {r['BioProject']:14s} {r['n']:4d} {r['localities']:4d} "
                  f"{r['years_known']:4d} {r['pmc'] or '-':13s} {r['doi'] or '-'}")


if __name__ == "__main__":
    main()
