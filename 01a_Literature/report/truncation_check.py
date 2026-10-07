#!/usr/bin/env python3
"""truncation_check.py — how much of each study the STAT gate actually kept.

A BioProject contributing one or two samples to the cohort looks like a tiny study. It
usually is not. Nobody sequences one sample: what happened is that the study was screened in
full and the STAT gate kept the handful of runs that carried enough pathogen signal. The rest
are still in SRA, unexamined.

This measures that, per BioProject, at three points:

    sra       runs in SRA, from NCBI esearch. The external truth. Needs --verify.
    screened  runs present in stat_cache.jsonl, so runs our SRA search found AND profiled.
    passed    runs in runs.tsv, so runs that cleared the kingdom thresholds.

`sra` versus `screened` is what the SRA *search* missed. `screened` versus `passed` is what
the *gate* removed. They are different failures and conflating them misattributes the bias;
measured on 2026-10-05 the search missed almost nothing (1,954 of 1,960 runs for the worst
case) and the gate did effectively all of the truncation.

The local pass needs no network and takes seconds. `--verify` adds two esearch calls per
BioProject, so pass a short list.

Usage:

    python 01a_Literature/report/truncation_check.py                      # local, all cohort projects
    python 01a_Literature/report/truncation_check.py --min-screened 100   # only the big ones
    python 01a_Literature/report/truncation_check.py --verify PRJNA383416 PRJNA306542
    python 01a_Literature/report/truncation_check.py --verify-worst 12    # verify the worst offenders

Reads `stat_cache.jsonl` (2.4 GB) by regex over the first bytes of each line. Do NOT
json.loads the whole line: the `_stat` payload is ~4 KB per run and parsing all 607,197 of
them is minutes of work for two fields that sit in the first 200 characters.
"""
import argparse, collections, csv, json, os, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _paths import ROOT

STAT_CACHE = ROOT / "01_stat/02_fetch/data/stat_cache.jsonl"
RUNS       = ROOT / "01_stat/03_filter/data/runs.tsv"
SAMPLES    = ROOT / "02_literature/03_classify/data/samples.tsv"
OUT        = Path(__file__).resolve().parent / "data/truncation.tsv"

_BP = re.compile(r'"BioProject":"([^"]+)"')
# BioProject sits within the first ~200 bytes of every line; widen only if that changes.
_HEAD = 220


def screened_counts() -> collections.Counter:
    """run count per BioProject across everything STAT profiled."""
    if not STAT_CACHE.exists():
        sys.exit(f"stat_cache.jsonl not found at {STAT_CACHE}")
    n = collections.Counter()
    with open(STAT_CACHE) as fh:
        for line in fh:
            m = _BP.search(line[:_HEAD])
            if m:
                n[m.group(1)] += 1
    return n


def passed_counts() -> tuple[collections.Counter, collections.Counter]:
    """(runs, distinct biosamples) per BioProject among gate survivors."""
    runs, bs = collections.Counter(), collections.defaultdict(set)
    with open(RUNS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            runs[r["BioProject"]] += 1
            bs[r["BioProject"]].add(r["BioSample"])
    return runs, collections.Counter({k: len(v) for k, v in bs.items()})


def cohort_projects() -> dict:
    """BioProject -> (biosamples in the classified cohort, host, title)."""
    out = {}
    with open(SAMPLES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bp = r["BioProject"]
            e = out.setdefault(bp, dict(n=0, host="", title=r.get("title", "")))
            e["n"] += 1
            e["host"] = e["host"] or r.get("llm_host_resolved") or r.get("stat_host") or ""
    return out


def _esearch_count(term: str, key: str) -> int:
    q = {"db": "sra", "term": term, "retmode": "json", "rettype": "count"}
    if key:
        q["api_key"] = key
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(url, timeout=60) as r:
        return int(json.loads(r.read())["esearchresult"]["count"])


def verify(bioprojects: list[str]) -> dict:
    """SRA run counts per BioProject, total and RNA-Seq-only.

    ENA's filereport endpoint was returning an HTML error page on 2026-10-05, so this goes to
    NCBI, which is authoritative for PRJNA accessions anyway. 9 req/s is the ceiling with a
    key; the sleep below is deliberately well under it.
    """
    key = os.environ.get("NCBI_API_KEY", "")
    if not key:
        print("warning: no NCBI_API_KEY; rate limit is 2.5 req/s. "
              "Use `source ~/.profile` first.", file=sys.stderr)
    out = {}
    for bp in bioprojects:
        try:
            n_all = _esearch_count(f"{bp}[BioProject]", key)
            time.sleep(0.15)
            n_rna = _esearch_count(f"{bp}[BioProject] AND RNA-Seq[Strategy]", key)
        except Exception as e:
            print(f"  {bp}: {type(e).__name__}: {e}", file=sys.stderr)
            time.sleep(0.4)
            continue
        out[bp] = (n_all, n_rna)
        time.sleep(0.15)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-screened", type=int, default=0,
                    help="only report BioProjects with at least this many runs screened")
    ap.add_argument("--verify", nargs="*", metavar="PRJ",
                    help="also query NCBI for these BioProjects' true run counts")
    ap.add_argument("--verify-worst", type=int, metavar="N",
                    help="query NCBI for the N BioProjects that lost the most runs")
    args = ap.parse_args()

    print("reading stat_cache.jsonl ...", file=sys.stderr)
    screened = screened_counts()
    pruns, pbs = passed_counts()
    cohort = cohort_projects()
    print(f"  {sum(screened.values()):,} runs screened across {len(screened):,} BioProjects",
          file=sys.stderr)

    recs = []
    for bp, e in cohort.items():
        s, p = screened.get(bp, 0), pruns.get(bp, 0)
        if s < args.min_screened:
            continue
        recs.append(dict(BioProject=bp, screened=s, gate_passed=p,
                         gate_biosamples=pbs.get(bp, 0), cohort_biosamples=e["n"],
                         gate_kept_pct=round(100 * p / s, 2) if s else "",
                         runs_unexamined=s - p, host=e["host"],
                         title=e["title"].replace("\t", " ")))
    recs.sort(key=lambda r: -r["runs_unexamined"])

    targets = list(args.verify or [])
    if args.verify_worst:
        targets += [r["BioProject"] for r in recs[:args.verify_worst]
                    if r["BioProject"] not in targets]
    if targets:
        print(f"verifying {len(targets)} BioProjects against NCBI ...", file=sys.stderr)
        got = verify(targets)
        for r in recs:
            if r["BioProject"] in got:
                r["sra_runs"], r["sra_rnaseq"] = got[r["BioProject"]]

    fields = ["BioProject", "screened", "gate_passed", "gate_biosamples", "cohort_biosamples",
              "gate_kept_pct", "runs_unexamined", "sra_runs", "sra_rnaseq", "host", "title"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(recs)

    tot_s = sum(r["screened"] for r in recs)
    tot_p = sum(r["gate_passed"] for r in recs)
    print(f"\n{len(recs)} BioProjects: {tot_s:,} runs screened, {tot_p:,} passed the gate "
          f"({100*tot_p/tot_s:.1f}%), {tot_s-tot_p:,} never examined")
    print(f"\n{'BioProject':<14}{'screened':>9}{'passed':>7}{'cohort':>7}{'kept':>7}  "
          f"{'host':<20}title")
    for r in recs[:25]:
        kept = f"{r['gate_kept_pct']}%" if r["gate_kept_pct"] != "" else "-"
        print(f"{r['BioProject']:<14}{r['screened']:>9}{r['gate_passed']:>7}"
              f"{r['cohort_biosamples']:>7}{kept:>7}  {r['host'][:19]:<20}{r['title'][:44]}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
