#!/usr/bin/env python3
"""Turn the PDFs we hold into text the release pipeline can read.

    python research/01_literature/retrieve/pdf_to_text.py
    python research/01_literature/retrieve/pdf_to_text.py --drop ~/Downloads

Reads  release/01_cohort/curation/fulltext_wanted.tsv   the papers still missing text
       studies/doi_<slug>/*.pdf                         PDFs already on disk
       --drop DIR                                       anything newly downloaded
Writes release/01_cohort/inputs/fulltext_manual/<slug>.txt

WHY THIS LIVES IN research/. The release tree holds what a reviewer can run. These PDFs are
paywalled, obtained through an institutional subscription, and no reviewer can fetch them, so
converting them is our step and not part of the reproducible pipeline. The release receives
the text as an input and records `source=manual`, which is the honest limit.

WHY NOT GROBID. GROBID returns TEI XML with sections marked up, and needs a Java service and
a ~500 MB container. The only thing its structure buys here is telling a data-availability
statement apart from a citation, and pymupdf plus the heading line that precedes a hit
recovers most of that for one pip install. Measured on nature25966: 19 pages, 57,586
characters, and both PRJNA383416 and SRP115041 found along with the "Data availability"
heading. If the section attribution turns out to matter more than that, GROBID is the upgrade.

A DOI becomes a filename by replacing every character outside [A-Za-z0-9._-] with "_", the
same rule the release scripts use, so a converted file lands where they look for it.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
ROOT = MODULE.parent
REPO = ROOT.parent
WANTED = REPO / "release" / "01_cohort" / "curation" / "fulltext_wanted.tsv"
DST = REPO / "release" / "01_cohort" / "inputs" / "fulltext_manual"
STUDIES = MODULE / "studies"
MIN_CHARS = 4000


def slug(doi):
    return re.sub(r"[^A-Za-z0-9._-]", "_", doi.strip())


def extract(path):
    import pymupdf
    with pymupdf.open(path) as doc:
        return "\n".join(page.get_text() for page in doc)


def norm(s):
    """Lowercase alphanumerics only, for comparing a title against a filename."""
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


def pdfs_for(doi, drop, title=""):
    """Every PDF on disk that plausibly belongs to this DOI.

    Matched three ways, because publisher downloads are not named after the DOI: the full
    slug, the DOI's suffix (`science.aba5435`), and that suffix's last dotted part
    (`aba5435`), which is what Science, Wiley and Nature typically hand you. Saves renaming
    twenty files by hand.
    """
    out = list((STUDIES / f"doi_{doi.strip().lower().replace('/', '_')}").glob("*.pdf"))
    if not drop:
        return out
    full = slug(doi).lower()
    suffix = doi.strip().lower().split("/", 1)[-1]
    tail = suffix.rsplit(".", 1)[-1]
    keys = [k for k in (full, suffix, tail) if len(k) > 5]
    # Wiley, Elsevier and the societies name a download after the article, not the DOI:
    # "JIPB - 2026 - Kang - A conserved fungal Egh16-like effector...". The first 40
    # characters of the normalised title are enough to place those, and nothing else.
    tkey = norm(title)[:40]
    for f in drop.rglob("*.pdf"):
        stem = f.stem.lower().replace("/", "_")
        if any(k in stem for k in keys) or (len(tkey) >= 24 and tkey in norm(f.stem)):
            out.append(f)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drop", type=Path,
                    help="a directory of freshly downloaded PDFs, matched to DOIs by filename")
    ap.add_argument("--force", action="store_true", help="re-convert even if text exists")
    args = ap.parse_args()

    try:
        import pymupdf  # noqa: F401
    except ImportError:
        sys.exit("pymupdf is needed: pip install --user pymupdf")
    if not WANTED.exists():
        sys.exit(f"no {WANTED}; run release/01_cohort/fetch_fulltext.py first")

    with open(WANTED, newline="", encoding="utf-8") as fh:
        wanted = [(r["doi"], r.get("title", ""))
                  for r in csv.DictReader(fh, delimiter="\t") if r["doi"]]
    DST.mkdir(parents=True, exist_ok=True)
    print(f"{len(wanted)} papers wanted")

    done = skipped = short = nopdf = 0
    for doi, title in wanted:
        out = DST / f"{slug(doi)}.txt"
        if out.exists() and not args.force:
            skipped += 1
            continue
        found = pdfs_for(doi, args.drop, title)
        if not found:
            nopdf += 1
            continue
        # The largest, where a study holds both the article and its supplement.
        text = extract(max(found, key=lambda p: p.stat().st_size))
        # A PDF of scanned pages yields almost nothing. Better to leave it missing than to
        # record an empty file as full text and have the gate count it as retrieved.
        if len(text) < MIN_CHARS:
            short += 1
            print(f"  {doi}: only {len(text)} chars, probably not text-layer. Left alone")
            continue
        out.write_text(text, encoding="utf-8")
        done += 1
        print(f"  {doi} -> {out.name}  {len(text):,} chars")

    print(f"\nconverted {done}, already had text {skipped}, no PDF found {nopdf}, "
          f"too little text {short}")
    if args.drop:
        # Name every PDF in the drop that matched nothing, so a mis-named file is visible
        # rather than silently ignored.
        claimed = {p.resolve() for doi, ti in wanted for p in pdfs_for(doi, args.drop, ti)}
        loose = [p for p in args.drop.rglob("*.pdf") if p.resolve() not in claimed]
        if loose:
            print(f"\n{len(loose)} PDF(s) in {args.drop} matched no wanted DOI:")
            for p in loose[:12]:
                print(f"  {p.name}")
            print("  rename one to its DOI with / replaced by _, e.g. "
                  "10.1126_science.aba5435.pdf")
    if done:
        print("Now run release/01_cohort/fetch_fulltext.py to take them into the manifest.")


if __name__ == "__main__":
    main()
