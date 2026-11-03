#!/usr/bin/env python3
"""scan_exposure.py — pull the sentences that bear on disease_exposure, per BioProject.

Does NOT assign a level. It extracts candidate evidence so a human assignment rests on quoted
sentences rather than on recall, and so a NEGATIVE scan is recorded as evidence in its own
right rather than inferred from silence. See curate/DISEASE-EXPOSURE.md for the levels.

The scan is deliberately generous and will over-match. A paper describing a glasshouse arm that
did not produce the sequenced samples still matches `indoor`, and rule 1 of the curation
document exists for exactly that case. Reading the quote is the point.

    python 01_literature/curate/scan_exposure.py            # all cohort projects
    python 01_literature/curate/scan_exposure.py PRJNA297424
"""

import csv
import html
import unicodedata
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
ROOT = MODULE.parent
sys.path.insert(0, str(MODULE))

COHORT = MODULE / "release" / "cohort_samples.tsv"
EVIDENCE = MODULE / "release" / "evidence"
STUDIES = MODULE / "studies"
OUT = MODULE / "data" / "silver" / "exposure_scan.tsv"

# Ordered by the level each supports. A hit does not assign the level; it surfaces a sentence.
CUES = [
    ("sheltered", r"rain[- ]?(out|shelter)|closed? on rainy|protect\w*\s+(?:the\s+)?field\s+from\s+rain|"
                  r"excluded? (?:the )?rain|rain exclusion|canopy (?:was |is |could be )?clos"),
    ("protected", r"fungicid|seed (?:was |were )?treat|plant protection|pesticid|"
                  r"azoxystrobin|tebuconazole|propiconazole|chlorothalonil|mancozeb"),
    ("indoor",    r"glasshouse|greenhouse|growth chamber|growth cabinet|polytunnel|"
                  r"controlled[- ]environment|phytotron"),
    ("open",      r"\bfield (?:experiment|trial|condition|site|plot)|open field|"
                  r"randomi[sz]ed complete block|research farm|experimental station"),
    ("tissue:representative", r"representative plant|avoid\w* edge|healthy[- ]looking|"
                              r"disease[- ]free plant|asymptomatic"),
    ("tissue:positional", r"flag leaf|days after anthesis|\bDAA\b|fully expanded leaf|"
                          r"\b(?:fourth|third|second|uppermost) leaf"),
    ("tissue:symptomatic", r"lesion|symptomatic|infected leaf|diseased (?:leaf|plant|tissue)"),
]
# Words whose ABSENCE is itself evidence. Counted per paper and reported, because
# "fungicide 0 across the full text" supports `open` where silence does not.
MIN_PDF_CHARS = 4000   # below this a PDF is scanned images, not text
NEGATIVE_SCAN = ["fungicide", "pesticide", "insecticide", "inoculat", "infection",
                 "lesion", "greenhouse", "growth chamber", "rain"]


# PDF typesetting renders fi, fl, ff as single ligature glyphs, and pypdf returns them as
# such. "field" comes back as "\ufb01eld" and "flag leaf" as "\ufb02ag leaf", so a regex
# containing either pair silently fails on every PDF while working on every XML. This cost the
# rain-exclusion evidence for PRJNA510920, whose methods say "protecting the \ufb01eld from
# rain". NFKD decomposes most of them; the explicit map covers what it leaves.
LIGATURES = {"\ufb00": "ff", "\ufb01": "fi", "\ufb02": "fl", "\ufb03": "ffi",
             "\ufb04": "ffl", "\ufb05": "st", "\ufb06": "st"}


def defang(s: str) -> str:
    for k, v in LIGATURES.items():
        s = s.replace(k, v)
    return unicodedata.normalize("NFKD", s)


def text_of(path: Path) -> str:
    raw = path.read_bytes()
    try:
        s = raw.decode("utf-8", "ignore")
    except Exception:                                   # noqa: BLE001
        return ""
    if path.suffix.lower() == ".pdf":
        # pypdf, with a guard. A scanned PDF with no text layer returns a few stray characters
        # per page, and quoting from that produces mangled sentences that look like evidence
        # and are not. Anything under MIN_PDF_CHARS is treated as unread so a human opens it.
        try:
            import pypdf
            r = pypdf.PdfReader(str(path))
            s = " ".join((pg.extract_text() or "") for pg in r.pages)
        except Exception as e:                          # noqa: BLE001
            print(f"  (unreadable pdf {path.name}: {type(e).__name__})", file=sys.stderr)
            return ""
        s = re.sub(r"\s+", " ", defang(s))
        if len(s) < MIN_PDF_CHARS:
            print(f"  (pdf has no usable text layer: {path.name}, {len(s)} chars)",
                  file=sys.stderr)
            return ""
        return s
    if path.suffix.lower() in (".xml", ".html", ".htm"):
        s = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", s, flags=re.S | re.I)
        s = re.sub(r"<[^>]+>", " ", s)
    return defang(html.unescape(re.sub(r"\s+", " ", s)))


MAX_QUOTE = 700


def sentences(t: str, pat: str, window: int = 3):
    """Quotes around each match, TRIMMED to length rather than discarded for being long.

    This dropped matches silently until 2026-10-21. The window was built out to sentence
    boundaries and then filtered on `40 < len(frag) < 700`, so a real cue inside a long
    sentence vanished with no warning. The rain-exclusion sentence in PRJNA609211, which is the
    single most important piece of evidence in this scan, built a 779-character window and was
    thrown away; the project reported no shelter cue while its own methods section said the
    canopy "could be closed on rainy days".

    A length limit on a QUOTE is a formatting preference. It must never decide whether evidence
    exists. Trim around the match instead, and keep the match in view.
    """
    out = []
    for m in re.finditer(pat, t, re.I):
        lo = t.rfind(".", 0, max(0, m.start() - 300))
        hi = t.find(".", m.end() + 300)
        a = lo + 1 if lo > 0 else max(0, m.start() - 300)
        b = hi + 1 if hi > 0 else min(len(t), m.end() + 300)
        frag = t[a:b].strip()
        if len(frag) > MAX_QUOTE:                       # trim AROUND the match, keep it centred
            half = MAX_QUOTE // 2
            c = (m.start() + m.end()) // 2
            frag = ("..." + t[max(a, c - half):min(b, c + half)].strip() + "...")
        if len(frag) > 20 and frag not in out:
            out.append(frag)
        if len(out) >= window:
            break
    return out


def main() -> None:
    want = sys.argv[1:] or None
    rows = list(csv.DictReader(open(COHORT), delimiter="\t"))
    projects = sorted({r["BioProject"] for r in rows})
    if want:
        projects = [p for p in projects if p in want]

    # Which study directories back each project, from the generated evidence files.
    dois = defaultdict(set)
    for p in projects:
        f = EVIDENCE / f"{p}.md"
        if f.exists():
            for m in re.finditer(r"10\.\d{4,9}/[^\s);]+", f.read_text(errors="ignore")):
                dois[p].add(m.group(0).rstrip(".,)"))

    out, summary = [], []
    for p in projects:
        texts, pdfs = [], []
        for d in dois[p]:
            key = "doi_" + d.replace("/", "_")
            sd = STUDIES / key
            if not sd.is_dir():
                continue
            for f in sd.iterdir():
                # EVERY text-bearing format, not just the two I first thought of. Leon
                # pushed back on a "no readable text" report for 10 projects and was right:
                # the studies tree holds 82 _fulltext.xml, 40 pdf, 22 _cachedtext.txt and 4
                # html, and only the first two were being read. Five dirs genuinely hold no
                # text; everything else was a file-extension oversight on my part.
                if f.name in ("SOURCE.txt", "NEEDED.md") or f.suffix.lower() in (
                        ".xlsx", ".csv", ".tsv", ".yaml", ".zip", ".docx", ".md"):
                    continue
                if f.suffix.lower() == ".pdf":
                    pdfs.append(f.name)
                if f.suffix.lower() in (".pdf", ".xml", ".txt", ".html", ".htm"):
                    t = text_of(f)
                    if t:
                        texts.append((d, f.name, t))
        if not texts:
            summary.append((p, "NO READABLE TEXT", len(pdfs), ""))
            out.append([p, "", "", "no_text", "", f"{len(pdfs)} pdf(s) present, not parsed"])
            continue

        hits = defaultdict(list)
        negs = []
        for doi, fn, t in texts:
            for level, pat in CUES:
                for frag in sentences(t, pat):
                    hits[level].append((doi, fn, frag))
            negs.append((doi, {w: len(re.findall(w, t, re.I)) for w in NEGATIVE_SCAN}))
        for level, v in hits.items():
            for doi, fn, frag in v:
                out.append([p, doi, fn, level, frag, ""])
        for doi, counts in negs:
            out.append([p, doi, "", "negative_scan", "",
                        ", ".join(f"{k} {n}" for k, n in counts.items())])
        summary.append((p, ",".join(sorted(hits)) or "none", len(pdfs), ""))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["BioProject", "doi", "file", "cue", "quote", "note"])
        w.writerows(out)

    print(f"-> {OUT.relative_to(ROOT)}  ({len(out)} rows over {len(projects)} projects)\n")
    print(f"  {'BioProject':<14}{'cues found':<52}{'pdfs'}")
    for p, c, n, _ in summary:
        print(f"  {p:<14}{c[:50]:<52}{n}")


if __name__ == "__main__":
    main()
