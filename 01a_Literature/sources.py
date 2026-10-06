#!/usr/bin/env python3
"""Fetch and cache the sources a curation decision rests on, using the right endpoint.

Written 2026-10-06. Every endpoint below was learned the hard way over the preceding days,
and the knowledge kept living only in a memory note, so each session re-derived the URL
patterns and re-made the same mistakes. This file is where that knowledge lives now.

Everything fetched is cached under `01a_Literature/studies/doi_<slug>/`, which is also what
`apply_curation.py --check-sources` looks for, so filing happens by construction rather than
as a step somebody remembers.

Subcommands, in the order a project usually needs them:

    papers   PRJNAxxxxxx            which papers name this accession
    fulltext PMCID                  fetch the full text, file it, list its sections
    sections PMCID|path  PATTERN    print only the sections whose title matches
    scan     PMCID|path             count the disease/treatment vocabulary; no context dump
    geo      GSExxxxx               per-sample characteristics, as value-counts
    ena      PRJxxxxxx              per-run fields, ONE request for the whole project

## What each source is good for, and how it fails

Europe PMC search  the only reliable accession->paper route. Note the quotes around the
                   accession: an unquoted query matches loosely.
fullTextXML        open-access full text. Reliable.
fullTextPDF        404s often, including for papers whose XML is served fine. Do not treat
                   a PDF 404 as the paper being unavailable.
supplementaryFiles serves TRUNCATED zips (seen four times in one day: PNAS, Nat Genet, BMC,
                   OUP). The truncated zip still lists member filenames, which is how to
                   learn what to ask the publisher host for. Not wrapped here because it
                   fails more often than it works.
Springer/BMC/Nature supplements, directly and reliably:
                   https://static-content.springer.com/esm/art%3A<DOI urlencoded>/MediaObjects/<file>
GEO                `?acc=GSE...&targ=gsm&form=text&view=brief` carries tissue, dates,
                   treatment and the `molecule` confirmation that SRA leaves blank. GEO
                   deposits are the EMPTIEST in SRA and the CHEAPEST to resolve. Serves a
                   reCAPTCHA if hit repeatedly in one session, so this caches hard.
ENA portal         one request per PROJECT, never per sample (~150x faster). Does not mirror
                   every NCBI project - all eight Sato Arabidopsis projects came back blank -
                   and `scientific_name` is the submitter's subject, not necessarily the host.
Paywalled          Leon can get them. Ask, naming the specific missing fact.
"""

from __future__ import annotations

import collections, json, re, sys, urllib.parse, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent / "01a_Literature"
STUDIES = HERE / "studies"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"

# Words whose ABSENCE is what licenses `sampling_selection = unselected`, and whose presence
# changes the stratum. Counted, never dumped: a context dump of 17 terms was the single most
# wasteful print in the 2026-10-06 session.
VOCAB = ["fungicide", "pesticide", "insecticide", "inocul", "disease", "pathogen", "infect",
         "symptom", "asymptomatic", "lesion", "sterili", "surface-steril", "greenhouse",
         "growth chamber", "glasshouse", "field", "plot", "pool"]


def get(url: str, timeout: int = 120) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "crypt-curation/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return fh.read().decode("utf-8", errors="replace")


def flatten(xml: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml))


def resolve(arg: str) -> tuple[str, Path | None]:
    """Accept a PMCID or a path. Returns (text, path)."""
    p = Path(arg)
    if p.exists():
        return p.read_text(encoding="utf-8", errors="replace"), p
    hits = sorted(STUDIES.rglob(f"*{arg}*")) if arg.startswith("PMC") else []
    if hits:
        return hits[0].read_text(encoding="utf-8", errors="replace"), hits[0]
    return get(f"{EPMC}/{arg}/fullTextXML"), None


def cmd_papers(acc: str) -> None:
    q = urllib.parse.quote(f'"{acc}"')
    d = json.loads(get(f"{EPMC}/search?query={q}&format=json&pageSize=25&resultType=core"))
    print(f"{acc}: {d['hitCount']} hit(s)")
    for r in d["resultList"]["result"]:
        print(f"  {r.get('doi','(no doi)')}  {r.get('pmcid','(no pmcid)')}  "
              f"OA={r.get('isOpenAccess')}  {r.get('pubYear')}")
        print(f"    {(r.get('title') or '')[:150]}")
        print(f"    {r.get('journalInfo',{}).get('journal',{}).get('title','')}")


def cmd_fulltext(pmcid: str) -> None:
    xml = get(f"{EPMC}/{pmcid}/fullTextXML")
    m = re.search(r'<article-id pub-id-type="doi">(.*?)</article-id>', xml)
    if not m:
        sys.exit(f"REFUSED: no DOI in {pmcid}'s XML; cannot decide where to file it")
    doi = m.group(1).strip()
    out = STUDIES / ("doi_" + doi.replace("/", "_"))
    out.mkdir(parents=True, exist_ok=True)
    f = out / f"{doi.split('/')[-1]}_fulltext.xml"
    f.write_text(xml)
    print(f"{doi}  {pmcid}  {len(xml):,} bytes -> {f.relative_to(HERE.parent)}")
    for t in re.findall(r"<title>(.*?)</title>", xml, re.S):
        t = re.sub(r"<[^>]+>", "", t).strip()
        if t:
            print(f"  SEC: {t}")


def cmd_sections(arg: str, pattern: str, cap: int = 2500) -> None:
    xml, _ = resolve(arg)
    rx = re.compile(pattern, re.I)
    for blk in xml.split("<title>")[1:]:
        title = re.sub(r"<[^>]+>", "", blk.split("</title>")[0]).strip()
        if not rx.search(title) or "</title>" not in blk:
            continue
        body = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", blk.split("</title>", 1)[1])).strip()
        print(f"{'='*20} {title}\n{body[:cap]}\n")


def cmd_scan(arg: str) -> None:
    xml, p = resolve(arg)
    t = flatten(xml).lower()
    print(f"scan {p or arg}: {len(t):,} chars")
    for w in VOCAB:
        n = t.count(w)
        if n:
            print(f"  {w:<16} {n}")
    absent = [w for w in VOCAB if w not in t]
    print(f"  ABSENT: {', '.join(absent) if absent else '(none)'}")


def cmd_geo(gse: str) -> None:
    cache = HERE / "data" / "geo" / f"{gse}_brief.txt"
    if cache.exists():
        txt = cache.read_text(encoding="utf-8", errors="replace")
        print(f"(cached {cache.name})")
    else:
        txt = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                  f"?acc={gse}&targ=gsm&form=text&view=brief", timeout=180)
        if "recaptcha" in txt.lower():
            sys.exit("REFUSED: GEO served a reCAPTCHA; wait, then retry")
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(txt)
        print(f"cached -> {cache.relative_to(HERE.parent)}")
    recs = []
    for blk in txt.split("^SAMPLE = ")[1:]:
        d = {"GSM": blk.split("\n", 1)[0].strip()}
        for line in blk.splitlines():
            m = re.match(r"!Sample_(\w+) = (.*)$", line)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if k == "characteristics_ch1":
                kk, _, vv = v.partition(":")
                d[kk.strip()] = vv.strip()
            elif k in ("source_name_ch1", "molecule_ch1", "title"):
                d[k] = v
        recs.append(d)
    print(f"{gse}: {len(recs)} samples")
    for k in sorted({k for r in recs for k in r} - {"GSM"}):
        c = collections.Counter(r.get(k, "") for r in recs)
        if len(c) <= 6:
            print(f"  {k:<22} " + "  ".join(f"{a or '(blank)'}:{n}" for a, n in c.most_common()))
        else:
            top = "  ".join(f"{a}:{n}" for a, n in c.most_common(3))
            print(f"  {k:<22} {len(c)} distinct, top: {top}")


def cmd_ena(acc: str) -> None:
    fields = ("run_accession,sample_accession,scientific_name,collection_date,country,"
              "location,sample_title,library_strategy,library_selection")
    txt = get("https://www.ebi.ac.uk/ena/portal/api/filereport"
              f"?accession={acc}&result=read_run&format=tsv&fields={fields}")
    lines = [l for l in txt.splitlines() if l.strip()]
    if len(lines) < 2:
        print(f"{acc}: ENA returned no rows. It does not mirror every NCBI project; "
              f"use the BioSample attribute cache or GEO instead")
        return
    hdr = lines[0].split("\t")
    rows = [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]
    cache = HERE / "data" / "ena" / f"{acc}.tsv"
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(txt)
    print(f"{acc}: {len(rows)} runs -> {cache.relative_to(HERE.parent)}")
    for k in hdr[1:]:
        c = collections.Counter(r.get(k, "") for r in rows)
        if len(c) <= 6:
            print(f"  {k:<20} " + "  ".join(f"{a or '(blank)'}:{n}" for a, n in c.most_common()))
        else:
            print(f"  {k:<20} {len(c)} distinct, top: "
                  + "  ".join(f"{a}:{n}" for a, n in c.most_common(3)))


CMDS = {"papers": (cmd_papers, 1), "fulltext": (cmd_fulltext, 1), "sections": (cmd_sections, 2),
        "scan": (cmd_scan, 1), "geo": (cmd_geo, 1), "ena": (cmd_ena, 1)}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in CMDS:
        sys.exit(__doc__)
    fn, n = CMDS[sys.argv[1]]
    args = sys.argv[2:]
    if len(args) < n:
        sys.exit(f"usage: sources.py {sys.argv[1]} " + " ".join(["<arg>"] * n))
    fn(*args)
