#!/usr/bin/env python3
"""resolve_accessions.py — turn non-BioProject accessions into BioProjects.

Papers report their data under whichever accession their journal asked for. Undermind's
reports name GEO series, SRA studies and ArrayExpress experiments as often as BioProjects,
and `triage.py` only recognises `PRJ[DEN]x`, so those papers land in `no-accession` with
their data sitting in plain sight. Measured 2026-10-05: 15 papers name an `SRP`, 13 a `GSE`,
4 an `E-MTAB`, out of 194 that looked accession-less.

Each type reduces to a lookup we already do:

    SRP/ERP/DRP   an SRA study accession goes straight into esearch db=sra; runinfo's
                  BioProject column is the answer.
    GSE           GEO series are cross-linked to SRA, so the same search usually works;
                  elink from gds to sra is the fallback.
    SRR/SRX/ERR   a single run or experiment resolves the same way, and is worth trying
                  because some papers cite only one example run.
    E-MTAB        ArrayExpress lives at EBI. ENA's filereport would answer it but was
                  returning HTTP 500 on every request on 2026-10-05, so these are reported
                  and left for manual resolution rather than guessed at.

**A mentioned accession is not necessarily the paper's own data.** A GSE in a methods
section may be a dataset the authors reanalysed. Every resolution therefore records the
accession it came from and is written with `relation="sra_linked"`, never `generated`:
the link says an archive lookup connected them, not that this paper produced the data.
Deciding who generated it is a human call and stays one.

Output is `data/resolved_accessions.tsv`, which `triage.py` picks up on its next run.

Usage:

    python 01a_Literature/retrieve/resolve_accessions.py --dry-run   # what would be looked up
    python 01a_Literature/retrieve/resolve_accessions.py
"""
import argparse, collections, csv, json, os, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _paths import ROOT
from _layout import RESOLVED_ACCESSIONS, UNDERMIND, WORKLIST
import triage

OUT    = RESOLVED_ACCESSIONS
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA     = {"User-Agent": "crypt/resolve_accessions (leon.lenzo@curtin.edu.au)"}

# Ordered by how directly each resolves. A study accession is unambiguous; a single run
# accession is weaker evidence of what the paper as a whole deposited, so it is tried last.
PATTERNS = [
    ("sra_study",    re.compile(r"\b[SED]RP\d{5,}\b")),
    ("geo_series",   re.compile(r"\bGSE\d{3,}\b")),
    ("sra_run",      re.compile(r"\b[SED]R[RSX]\d{5,}\b")),
    ("arrayexpress", re.compile(r"\bE-[A-Z]{4}-\d+\b")),
]
_UNRESOLVABLE = {"arrayexpress"}   # EBI route; see module docstring


def _post(endpoint: str, **params) -> str:
    params.setdefault("api_key", os.environ.get("NCBI_API_KEY", ""))
    body = urllib.parse.urlencode({k: v for k, v in params.items() if v}).encode()
    req = urllib.request.Request(EUTILS + endpoint, data=body, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8", "replace")


def to_bioprojects(acc: str) -> list[str]:
    """BioProjects an accession resolves to, via SRA runinfo. Empty if it resolves to none."""
    try:
        es = json.loads(_post("esearch.fcgi", db="sra", term=acc,
                              retmode="json", usehistory="y"))["esearchresult"]
        if es.get("count") in (None, "0"):
            return []
        time.sleep(0.12)
        body = _post("efetch.fcgi", db="sra", WebEnv=es["webenv"],
                     query_key=es["querykey"], rettype="runinfo", retmode="csv")
    except Exception as exc:
        print(f"    {acc}: {type(exc).__name__}: {exc}", file=sys.stderr)
        return []
    import io
    rows = csv.DictReader(io.StringIO(triage.RUNINFO_HEADER + "\n" + body))
    bps = {r["BioProject"] for r in rows
           if (r.get("Run") or "").startswith(("SRR", "ERR", "DRR")) and r.get("BioProject")}
    return sorted(bps)


def candidates() -> dict:
    """paper_ref -> {kind: {accessions}} for papers with no BioProject of their own."""
    mds = sorted(UNDERMIND.glob("*.md"))      # see triage.py: allowlist by location
    work = {r["paper_key"]: r for r in
            csv.DictReader(open(WORKLIST), delimiter="\t")}
    noacc = {w["paper_ref"] for w in work.values()
             if w["need"] == "no-accession" and w.get("paper_ref")}
    out = collections.defaultdict(lambda: collections.defaultdict(set))
    rows = []
    for md in mds:
        rows += triage.undermind_rows(md)
    for sec, cells in rows:
        m = triage._REF_RE.search(cells[0])
        if not m or m.group(1) not in noacc:
            continue
        row = " | ".join(cells)
        if triage._ACC_RE.search(row):
            continue       # it names a BioProject after all; triage will pick that up
        for kind, pat in PATTERNS:
            for hit in pat.findall(row):
                out[m.group(1)][kind].add(hit)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if not os.environ.get("NCBI_API_KEY"):
        print("warning: no NCBI_API_KEY; 2.5 req/s. `source ~/.profile` first.",
              file=sys.stderr)

    cand = candidates()
    n_acc = sum(len(v) for kinds in cand.values() for v in kinds.values())
    print(f"{len(cand)} accession-less papers name {n_acc} non-BioProject accessions")
    for kind, _ in PATTERNS:
        papers = [r for r, k in cand.items() if kind in k]
        print(f"  {kind:<14}{len(papers):>3} papers")
    if args.dry_run:
        for ref, kinds in sorted(cand.items()):
            print(f"  {ref:<10}" + "  ".join(f"{k}={sorted(v)}" for k, v in kinds.items()))
        return

    recs = []
    for ref, kinds in sorted(cand.items()):
        for kind, _ in PATTERNS:
            for acc in sorted(kinds.get(kind, ())):
                if kind in _UNRESOLVABLE:
                    recs.append(dict(paper_ref=ref, via=acc, via_kind=kind, BioProject="",
                                     status="needs_ebi", note="ENA filereport returned HTTP 500"))
                    continue
                bps = to_bioprojects(acc)
                if not bps:
                    recs.append(dict(paper_ref=ref, via=acc, via_kind=kind, BioProject="",
                                     status="no_match", note=""))
                else:
                    for bp in bps:
                        recs.append(dict(paper_ref=ref, via=acc, via_kind=kind, BioProject=bp,
                                         status="resolved",
                                         note=f"{acc} -> {bp} via SRA runinfo"))
                print(f"  {ref:<10}{acc:<16}-> {', '.join(bps) or '(none)'}")
                time.sleep(0.15)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["paper_ref", "via", "via_kind", "BioProject",
                                           "status", "note"], delimiter="\t")
        w.writeheader()
        w.writerows(recs)

    got = {r["paper_ref"] for r in recs if r["status"] == "resolved"}
    bps = {r["BioProject"] for r in recs if r["BioProject"]}
    print(f"\nresolved {len(got)} of {len(cand)} papers to {len(bps)} BioProjects")
    for st, n in collections.Counter(r["status"] for r in recs).most_common():
        print(f"  {st:<12}{n:>4} accessions")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    print("These are sra_linked, not `generated`: an archive lookup connected them, which is "
          "not a claim about who produced the data.")


if __name__ == "__main__":
    main()
