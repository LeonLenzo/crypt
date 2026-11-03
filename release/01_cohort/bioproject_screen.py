#!/usr/bin/env python3
"""Find the data accessions each paper reports, in its full text and abstract.

    python release/01_cohort/bioproject_screen.py

Reads  raw/fulltext/*.xml           where full text was retrieved
       raw/raw_papers.tsv           the abstract, for papers without it
Writes curation/bioproject_screen_papers.tsv   one row per (paper, accession)

WHERE a paper names an accession decides how much it is worth, so the row records it. An
accession in a data-availability statement is the paper's own deposit; the same accession in
the introduction is usually somebody else's work being cited. The enclosing section title and
the surrounding text travel with every hit so the call can be checked rather than trusted.

Accession types are kept apart rather than merged. Only PRJ* is a BioProject; GSE, SRP and
E-MTAB are study accessions in other registries that have to be resolved to one, and SRR/ERR
and SAM* are runs and samples that imply a project without naming it. Collapsing them would
silently invent BioProjects that the paper never gave.
"""

import collections
import csv
import glob
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
FULLTEXT = HERE / "raw" / "fulltext"
PAPERS = HERE / "raw" / "raw_papers.tsv"
OUT = HERE / "curation" / "bioproject_screen_papers.tsv"

PATTERNS = {
    "bioproject": r"PRJ(?:NA|EB|DB)\d{4,}",
    "sra_study": r"\b[SED]RP\d{5,}",
    "geo_series": r"\bGSE\d{4,}",
    "arrayexpress": r"\bE-[A-Z]{4}-\d+",
    "run": r"\b[SED]RR\d{5,}",
    "biosample": r"\bSAM[NED][A-Z]?\d{5,}",
}
ACCESSION_RE = re.compile("|".join(f"(?P<{k}>{v})" for k, v in PATTERNS.items()))

# Section headings that mean "this is our deposit". Matched against the enclosing <title>.
DEPOSIT_SEC = re.compile(r"data availab|availability of data|accession|data deposition|"
                         r"supporting information|data access", re.I)

TAG = re.compile(r"<[^>]+>")
SEC = re.compile(r"<sec\b[^>]*>|</sec>|<title[^>]*>(.*?)</title>", re.S | re.I)


def sections(xml):
    """Yield (section title, plain text) for each <sec>, plus the rest of the article."""
    out, stack, pos = [], [], 0
    for m in re.finditer(r"<sec\b[^>]*>(.*?)</sec>", xml, re.S | re.I):
        block = m.group(1)
        t = re.search(r"<title[^>]*>(.*?)</title>", block, re.S | re.I)
        title = TAG.sub(" ", t.group(1)) if t else ""
        out.append((" ".join(title.split()), TAG.sub(" ", block)))
    out.append(("", TAG.sub(" ", xml)))
    return out


def scan(text, where, doi, seen, rows):
    text = " ".join(text.split())
    for m in ACCESSION_RE.finditer(text):
        kind = m.lastgroup
        acc = m.group(0)
        if (doi, acc) in seen:
            continue
        seen.add((doi, acc))
        a, b = max(0, m.start() - 90), min(len(text), m.end() + 90)
        rows.append({
            "doi": doi, "accession": acc, "type": kind,
            "section": where,
            "in_deposit_section": int(bool(DEPOSIT_SEC.search(where))),
            "context": text[a:b].strip(),
        })


def main():
    abstracts = {}
    with open(PAPERS, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            abstracts[r["doi"]] = r.get("abstract", "")

    files = {Path(p).stem.replace("_", "/", 1): p
             for p in glob.glob(str(FULLTEXT / "*.xml"))}
    rows, searched = [], collections.Counter()
    for doi in abstracts:
        seen = set()
        path = next((p for s, p in files.items() if Path(p).stem == re.sub(
            r"[^A-Za-z0-9._-]", "_", doi)), None)
        if path:
            xml = Path(path).read_text(errors="ignore")
            for title, text in sections(xml):
                scan(text, title, doi, seen, rows)
            searched["fulltext"] += 1
        elif abstracts[doi]:
            scan(abstracts[doi], "abstract", doi, seen, rows)
            searched["abstract only"] += 1
        else:
            searched["nothing to search"] += 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["doi", "accession", "type", "section", "in_deposit_section", "context"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    kinds = collections.Counter(r["type"] for r in rows)
    papers = len({r["doi"] for r in rows})
    bp = len({r["accession"] for r in rows if r["type"] == "bioproject"})
    print(f"-> {OUT.name}  {len(rows)} hits across {papers} papers")
    print(f"   searched: {dict(searched)}")
    print(f"   by type:  {dict(kinds.most_common())}")
    print(f"   distinct BioProjects named: {bp}")
    print(f"   hits in a data-availability section: "
          f"{sum(1 for r in rows if r['in_deposit_section'])}")


if __name__ == "__main__":
    main()
