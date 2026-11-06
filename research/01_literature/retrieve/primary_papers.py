#!/usr/bin/env python3
"""primary_papers.py — for a BioProject with no credited paper, who actually made it?

    python 01_literature/retrieve/primary_papers.py --accessions PRJEB15280
    python 01_literature/retrieve/primary_papers.py --orphans

`accession_papers.py` answers "which papers MENTION this accession", which is the right
question for finding umbrella siblings and the wrong one for attribution: a reanalysis
mentions an accession exactly as its depositor does. This adds the two things that
separate them.

**The archive record says which GROUP deposited.** ENA serves INSDC project XML for PRJEB,
PRJNA and PRJDB alike, carrying `center_name`, the study title the submitter wrote, and
ENA-FIRST-PUBLIC. Note what it does NOT reliably carry: of 23 cohort BioProjects checked
on 2026-11-05, exactly one had a PUBMED link in either ENA PROJECT_LINKS or NCBI elink.
Do not expect the archive to name the paper; expect it to name the lab and the date.

**The release date bounds the claim.** Data is deposited at or before publication, so the
submitter's paper sits near ENA-FIRST-PUBLIC and a paper years later is reanalysing. This
is what exposes a browser paper claiming nine BioProjects as its own: PRJEB31334 went
public in 2019 and PRJNA256347 in 2014, against a 2021 paper.

Neither is proof, so this writes candidates ranked earliest first with the centre and the
date beside them, and the pick stays a curation call. See CURATING.md.
"""
import argparse, csv, re, sys, time, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
ROOT = MODULE.parent
sys.path.insert(0, str(ROOT))

MENTIONS = MODULE / "data" / "silver" / "accession_papers.tsv"
OUT = MODULE / "data" / "silver" / "primary_papers_worklist.tsv"
COHORT = ROOT.parent / "release" / "01_cohort" / "data" / "outputs" / "supp_data_s2_cohort.tsv"
STUDIES = ROOT.parent / "release" / "01_cohort" / "data" / "outputs" / "supp_data_s1_studies.tsv"
UA = {"User-Agent": "crypt-cohort/1.0 (leon.lenzo@curtin.edu.au)"}
COLUMNS = ["accession", "runs", "centre", "first_public", "study", "candidates", "rank",
           "year", "doi", "journal", "in_screen", "paper"]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def split(field):
    return [a.strip() for a in field.replace(";", ",").split(",") if a.strip()]


def deposit(acc):
    """Submitting centre, study title and release date, from the INSDC project record."""
    time.sleep(0.5)   # EBI shares a 2 req/s limiter across its services
    try:
        body = urllib.request.urlopen(urllib.request.Request(
            f"https://www.ebi.ac.uk/ena/browser/api/xml/{acc}", headers=UA),
            timeout=60).read()
        project = ET.fromstring(body).find("PROJECT")
    except Exception as exc:
        return {"centre": f"<{type(exc).__name__}>", "study": "", "first_public": ""}
    if project is None:
        return {"centre": "", "study": "", "first_public": ""}
    public = next((a.findtext("VALUE") or "" for a in project.iter("PROJECT_ATTRIBUTE")
                   if a.findtext("TAG") == "ENA-FIRST-PUBLIC"), "")
    return {"centre": project.get("center_name", ""),
            "study": (project.findtext("TITLE") or "").strip(),
            "first_public": public}


def uncredited():
    """Cohort BioProjects with no paper in the screen that credibly generated them.

    `bioprojects_own` is taken at face value except where one paper claims several
    deposits, which is the shape of a browser or meta-analysis paper rather than a
    submitter, and is the case this script exists to resolve.
    """
    studies, cohort = read(STUDIES), read(COHORT)
    runs = {}
    for row in cohort:
        runs[row["BioProject"]] = runs.get(row["BioProject"], 0) + 1
    credited = set()
    for s in studies:
        own = [b for b in split(s["bioprojects_own"]) if b in runs]
        if len(own) == 1:
            credited.update(own)
    return [b for b in sorted(runs, key=lambda b: -runs[b]) if b not in credited], runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", nargs="*", default=[], metavar="PRJ")
    ap.add_argument("--orphans", action="store_true",
                    help="every cohort BioProject with no single-deposit paper crediting it")
    args = ap.parse_args()

    accessions, runs = (args.accessions, {})
    if args.orphans or not accessions:
        accessions, runs = uncredited()
    if not accessions:
        sys.exit("no accessions; pass --accessions or --orphans")

    screened = {s["doi"].lower() for s in read(STUDIES) if s["verdict"] == "keep"}
    mentions = {}
    for row in read(MENTIONS):
        mentions.setdefault(row["accession"], []).append(row)

    rows = []
    print(f"{len(accessions)} BioProjects to attribute\n")
    for n, acc in enumerate(accessions, 1):
        record = deposit(acc)
        hits = sorted(mentions.get(acc, []), key=lambda r: r.get("year") or "9999")
        base = dict(accession=acc, runs=runs.get(acc, ""), candidates=len(hits), **record)
        if not hits:
            rows.append({**base, "rank": "", "year": "", "doi": "", "journal": "",
                         "in_screen": "", "paper": ""})
        for rank, hit in enumerate(hits, 1):
            rows.append({**base, "rank": rank, "year": hit.get("year", ""),
                         "doi": hit.get("doi", ""), "journal": hit.get("journal", ""),
                         "in_screen": "yes" if hit.get("doi", "").lower() in screened else "",
                         "paper": hit.get("title", "")})
        earliest = hits[0].get("year", "?") if hits else "none"
        print(f"  [{n}/{len(accessions)}] {acc:<14}{len(hits):>3} candidates, "
              f"earliest {earliest:<6} public {record['first_public']:<11}"
              f"{record['centre'][:34]}")

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    attributed = {r["accession"] for r in rows if r["doi"]}
    covered = sum(runs.get(a, 0) for a in attributed)
    print(f"\nwrote {OUT}")
    print(f"{len(attributed)}/{len(accessions)} have at least one candidate, "
          f"{covered:,} of {sum(runs.get(a, 0) for a in accessions):,} runs")
    blank = [a for a in accessions if a not in attributed]
    if blank:
        print(f"no paper names these at all, so treat as unpublished deposits until "
              f"shown otherwise:\n  {', '.join(blank)}")
    print("\nA centre and a date are not an author list. Rank 1 is a candidate, not a "
          "verdict.")


if __name__ == "__main__":
    main()
