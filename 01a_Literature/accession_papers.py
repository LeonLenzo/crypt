#!/usr/bin/env python3
"""accession_papers.py — which papers mention a given accession.

The reverse of everything else here. `meta_search.py` and `resolve_accessions.py` go from a
paper to its data; this goes from the data back to every paper that cites it, by full-text
search over Europe PMC.

It exists because submitters share BioProjects. Wei et al. 2016 states that its 12 RNA-seq
pools sit inside PRJNA306542, which holds 433 runs; a full-text search on that accession
returns **eleven** papers, and their titles account for metadata that was otherwise
inexplicable. The `ecotype: upland/lowland` field belongs to a transcriptomic-divergence
paper, the nine `shallow root tips` samples to a root-transcriptomics paper, neither of which
is Wei et al.

Two things this fixes:

1. **Attribution.** A single link saying "this paper relates to this BioProject" lets one
   paper claim runs belonging to ten others. Finding the siblings is what makes a correct
   `run_scope` possible.
2. **Coverage.** An umbrella project may be mostly in scope even when the paper we found is
   not. Judging the project by the one paper we happened to retrieve throws away the rest.

Europe PMC is the right instrument because it matches the accession STRING in full text, so
it finds papers no semantic search would connect. It is not a claim about who generated the
data: a paper mentioning an accession may be reanalysing it, which is the same caution that
applies to `relation` in the link table.

EBI shares a 2 request/second limiter across its services, hence the sleep.

Usage:

    python 01a_Literature/accession_papers.py --accessions PRJNA306542
    python 01a_Literature/accession_papers.py --umbrella-candidates
    python 01a_Literature/accession_papers.py --umbrella-candidates --min-runs 100
"""
import argparse, collections, csv, json, sys, time, urllib.parse, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _paths import ROOT
from _layout import ACCESSION_PAPERS, BIOPROJECTS, PAPERS, PAPER_BIOPROJECT

OUT  = ACCESSION_PAPERS
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA   = {"User-Agent": "crypt/accession_papers (leon.lenzo@curtin.edu.au)"}


def papers_for(acc: str, page_size: int = 100) -> list[dict]:
    """Every Europe PMC record whose text contains this accession."""
    q = urllib.parse.urlencode({"query": f'"{acc}"', "format": "json",
                                "pageSize": page_size, "resultType": "core"})
    req = urllib.request.Request(f"{EPMC}?{q}", headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=90) as r:
            d = json.loads(r.read())
    except Exception as exc:
        print(f"    {acc}: {type(exc).__name__}: {exc}", file=sys.stderr)
        return []
    out = []
    for r in d.get("resultList", {}).get("result", []):
        out.append(dict(
            accession=acc,
            pmid=r.get("pmid", ""), doi=(r.get("doi") or "").lower(),
            year=r.get("pubYear", ""), journal=r.get("journalTitle", ""),
            # Europe PMC returns HTML entities and italics markup in titles
            title=(r.get("title", "") or "").replace("&lt;i&gt;", "")
                                            .replace("&lt;/i&gt;", "").strip(),
            is_oa=r.get("isOpenAccess", ""), source=r.get("source", ""),
        ))
    return out


def umbrella_candidates(min_runs: int) -> list[str]:
    """Projects big enough, and thinly enough claimed, to be shared between papers.

    The signal is a mismatch between size and attention: many runs but one or two papers
    linked. It is a heuristic, not a test, and it is cheap to over-include because the
    lookup costs one request.
    """
    reg = {r["BioProject"]: r for r in
           csv.DictReader(open(BIOPROJECTS), delimiter="\t")}
    claims = collections.Counter()
    conflicts = set()
    with open(PAPER_BIOPROJECT) as fh:
        for l in csv.DictReader(fh, delimiter="\t"):
            claims[l["BioProject"]] += 1
    for bp, r in reg.items():
        if r.get("primary_conflict"):
            conflicts.add(bp)
    out = [bp for bp, r in reg.items()
           if int(r["sra_runs"] or 0) >= min_runs and claims[bp] <= 2]
    return sorted(set(out) | conflicts)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", nargs="*", default=[], metavar="PRJ")
    ap.add_argument("--umbrella-candidates", action="store_true")
    ap.add_argument("--min-runs", type=int, default=50)
    args = ap.parse_args()

    accs = list(dict.fromkeys(args.accessions))
    if args.umbrella_candidates:
        accs = list(dict.fromkeys(accs + umbrella_candidates(args.min_runs)))
    if not accs:
        sys.exit("no accessions; pass --accessions or --umbrella-candidates")

    known = {r["doi"].lower() for r in
             csv.DictReader(open(PAPERS), delimiter="\t") if r.get("doi")}
    print(f"{len(accs)} accessions to look up; {len(known)} DOIs already in papers.tsv\n")

    recs = []
    for i, acc in enumerate(accs, 1):
        got = papers_for(acc)
        new = [p for p in got if p["doi"] and p["doi"] not in known]
        for p in got:
            p["already_known"] = "" if p["doi"] and p["doi"] not in known else "yes"
            recs.append(p)
        flag = "  <-- SHARED" if len(got) > 1 else ""
        print(f"  [{i}/{len(accs)}] {acc:<14}{len(got):>3} papers, {len(new):>3} new{flag}")
        time.sleep(0.5)

    # MERGE, never replace. Running this for two accessions once wiped the 69 rows an earlier
    # umbrella sweep had found, and the loss was silent: the file still looked valid, just
    # smaller, and the link table built from it reported PRJNA306542 as a one-paper project.
    prior = []
    if OUT.exists():
        with open(OUT) as fh:
            prior = list(csv.DictReader(fh, delimiter="\t"))
    fresh = {(r["accession"], r["doi"]) for r in recs}
    kept = [r for r in prior if (r["accession"], r.get("doi", "")) not in fresh]
    if kept:
        print(f"  merging with {len(kept)} existing rows for other accessions")
    recs = recs + kept

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "doi", "pmid", "year", "journal",
                                            "title", "is_oa", "source", "already_known"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(recs)

    by_acc = collections.Counter(r["accession"] for r in recs)
    shared = {a: n for a, n in by_acc.items() if n > 1}
    print(f"\n{len(recs)} paper-accession mentions across {len(by_acc)} accessions")
    print(f"  accessions mentioned by more than one paper: {len(shared)}")
    newdois = {r["doi"] for r in recs if not r["already_known"] and r["doi"]}
    print(f"  DOIs not yet in papers.tsv: {len(newdois)}")
    print("\nmost shared:")
    for a, n in sorted(shared.items(), key=lambda kv: -kv[1])[:10]:
        print(f"  {a:<14}{n:>3} papers")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    print("A mention is not authorship: some of these are reanalyses. The run_scope and "
          "relation columns are where that gets decided, by a human.")


if __name__ == "__main__":
    main()
