#!/usr/bin/env python3
"""Build the paper and accession curation worklists for release/01_cohort/curation/.

    python 01_literature/report/make_curation_worklists.py

Two decisions in this project were never written down: which screened manuscripts are
candidates, and which examined accessions carry the cohort. S1's 141 papers came from a host
filter, not from a decision, and it disagrees with the nine `undermind_decision` values that
were recorded. `in_cohort` in bioprojects.tsv says yes for 317 of 421 accessions while the
cohort is 41, so it is not that decision either.

This writes one row per paper and per accession with `decision` and `reason` blank, and every
scraped fact that bears on the call beside it. A reason is pre-filled ONLY where it follows
from scraped data and cannot be wrong: no runs, no RNA-seq runs, no accessions reported.

NOT pre-filled: anything resting on the archive's `organism` field. It names the host for only
27 of the 41 cohort accessions and misses 14 rust surveys deposited under *Puccinia*, so an
automatic "host is not a target crop" would be wrong on a third of the real cohort. The
organism travels as a column to read, not as a verdict.

RE-RUNNING IS SAFE. Existing decisions and reasons are read back and preserved; only the
scraped columns are refreshed and new rows appended. `_layout.py` records what happens
otherwise: hand annotation here was once discarded by the next triage run.
"""

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
ROOT = MODULE.parent
SILVER = MODULE / "data" / "silver"
REPO = ROOT.parent          # research/ and release/ are siblings
CURATION = REPO / "release" / "01_cohort" / "curation"


def read(path):
    if not Path(path).exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    done = sum(1 for r in rows if (r.get("decision") or "").strip())
    pre = sum(1 for r in rows if (r.get("reason") or "").strip())
    print(f"  {path.name:<16} {len(rows):>5} rows   {done:>4} decided   "
          f"{pre:>4} with a reason   {len(rows) - done:>4} TO DO")


def keep_decisions(path, key):
    """Existing hand decisions, so a re-run cannot discard them."""
    return {r[key]: (r.get("decision", ""), r.get("reason", ""))
            for r in read(path) if r.get(key) and (r.get("decision") or "").strip()}


def papers_worklist():
    out = CURATION / "curated_papers.tsv"
    prior = keep_decisions(out, "paper_key")
    bp = {r["BioProject"]: r for r in read(SILVER / "bioprojects.tsv")}
    cohort = {r["BioProject"] for r in read(MODULE / "release" / "cohort_bioprojects.tsv")}
    links = {}
    for r in read(SILVER / "paper_bioproject.tsv"):
        links.setdefault((r.get("paper_ref") or "").strip(), []).append(r["BioProject"])

    rows = []
    for p in read(SILVER / "papers.tsv"):
        key = p.get("paper_key", "")
        accs = sorted({a for k in {key, (p.get("doi") or "").strip(),
                                   (p.get("paper_ref") or "").strip()} if k
                       for a in links.get(k, [])})
        orgs = sorted({(bp.get(a, {}).get("organism") or "").strip() for a in accs} - {""})
        decision, reason = prior.get(key, ("", ""))
        if not reason and not accs:
            reason = "reports no accessions"
        rows.append({
            "paper_key": key, "doi": p.get("doi", ""),
            "title": " ".join((p.get("title") or "").split())[:200],
            "decision": decision, "reason": reason,
            "found_by": p.get("sources", ""),
            "undermind_view": p.get("undermind_decision", ""),
            "n_accessions": len(accs), "accessions": ";".join(accs),
            "accession_organisms": ";".join(orgs),
            "accessions_in_cohort": sum(1 for a in accs if a in cohort),
        })
    rows.sort(key=lambda r: (-r["accessions_in_cohort"], -r["n_accessions"], r["paper_key"]))
    write(out, rows, ["paper_key", "doi", "title", "decision", "reason", "found_by",
                      "undermind_view", "n_accessions", "accessions",
                      "accession_organisms", "accessions_in_cohort"])
    return rows


# Hosts and scope, used ONLY to attach a stated-criterion reason where the archive's own
# organism field settles it. A microbe organism settles nothing, because the host cannot be
# read off it: that is how 14 rust surveys deposited under *Puccinia* stayed in the cohort.
TARGET = re.compile(r"triticum|zea mays|oryza sativa|sorghum", re.I)
OTHER_CEREAL = re.compile(r"hordeum|avena|secale|triticale|triticosecale|aegilops|brachypodium|"
                          r"setaria|panicum|pennisetum|eleusine", re.I)
MICROBE = re.compile(r"puccinia|blumeria|zymoseptoria|pyrenophora|parastagonospora|fusarium|"
                     r"magnaporthe|pyricularia|ustilago|rhizoctonia|bipolaris|alternaria|"
                     r"colletotrichum|phytophthora|botrytis|verticillium|sclerotinia|claviceps|"
                     r"tilletia|phakopsora|melampsora|ramularia|cochliobolus|metagenome|virus|"
                     r"viroid|bacteri|xanthomonas|pseudomonas|aspergillus|penicillium|"
                     r"trichoderma|mycosphaerella|septoria|"
                     # Viral taxonomy names a pathogen, never the host. The realm and the
                     # standard rank suffixes are included because PRJNA796936 is deposited
                     # as `Riboviria` and a bare `virus` pattern missed it, handing a real
                     # cohort accession the reason "host is not a cereal".
                     r"riboviria|orthornavirae|monodnaviria|viria$|virae\b|virales|viridae|"
                     r"virinae|phage", re.I)


def criterion_reason(organism: str) -> str:
    """A reason that follows from a criterion already stated in the Methods, or "".

    Returns nothing for a microbe organism or a blank one: neither settles the host, and an
    automatic host exclusion there would be wrong on exactly the accessions that matter.
    """
    org = (organism or "").strip()
    if not org or MICROBE.search(org) or TARGET.search(org):
        return ""
    if OTHER_CEREAL.search(org):
        return "host is a cereal but not one of the four target crops"
    return "host is not a cereal"


def projects_worklist():
    out = CURATION / "curated_projects.tsv"
    prior = keep_decisions(out, "BioProject")
    cohort = {r["BioProject"] for r in read(MODULE / "release" / "cohort_bioprojects.tsv")}
    papers = {}
    for r in read(SILVER / "paper_bioproject.tsv"):
        papers.setdefault(r["BioProject"], []).append((r.get("paper_ref") or "").strip())

    def n(r, k):
        try:
            return int(r.get(k) or 0)
        except ValueError:
            return 0

    rows = []
    for b in read(SILVER / "bioprojects.tsv"):
        acc = b["BioProject"]
        decision, reason = prior.get(acc, ("", ""))
        source = ""
        # A criterion can never explain an accession that IS in the cohort. Suppressing it
        # here is belt and braces beside the host patterns, and it is what caught
        # PRJNA796936 (deposited as `Riboviria`) being handed a host-based exclusion.
        in_cohort = acc in cohort
        if not reason and not in_cohort:
            if n(b, "sra_runs") == 0:
                reason, source = "no runs in SRA", "criterion"
            elif n(b, "sra_rnaseq") == 0:
                reason, source = "no RNA-seq runs", "criterion"
            else:
                reason = criterion_reason(b.get("organism", ""))
                source = "criterion" if reason else ""
        elif acc in prior:
            source = "hand"
        rows.append({
            "BioProject": acc, "decision": decision, "reason": reason,
            "reason_source": source,
            "organism": b.get("organism", ""), "centre": b.get("centre", ""),
            "found_by": b.get("accession_source", ""),
            "sra_runs": b.get("sra_runs", ""), "sra_rnaseq": b.get("sra_rnaseq", ""),
            "sra_biosamples": b.get("sra_biosamples", ""),
            "triage": b.get("triage", ""),
            "undermind_view": b.get("undermind_decision", ""),
            "papers": ";".join(sorted({x for x in papers.get(acc, []) if x})),
            "in_cohort_now": "yes" if acc in cohort else "",
        })
    rows.sort(key=lambda r: (r["in_cohort_now"] != "yes", -n(r, "sra_rnaseq"), r["BioProject"]))
    write(out, rows, ["BioProject", "decision", "reason", "reason_source", "organism",
                      "centre", "found_by",
                      "sra_runs", "sra_rnaseq", "sra_biosamples", "triage",
                      "undermind_view", "papers", "in_cohort_now"])
    return rows


def main():
    print(f"-> {CURATION.relative_to(ROOT)}")
    papers_worklist()
    projects_worklist()


if __name__ == "__main__":
    main()
