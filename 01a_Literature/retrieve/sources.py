#!/usr/bin/env python3
"""Fetch and cache the sources a curation decision rests on, using the right endpoint.

Written 2026-10-06. Every endpoint below was learned the hard way over the preceding days,
and the knowledge kept living only in a memory note, so each session re-derived the URL
patterns and re-made the same mistakes. This file is where that knowledge lives now.

Everything fetched is cached under `01a_Literature/studies/doi_<slug>/`, which is also what
`apply_curation.py --check-sources` looks for, so filing happens by construction rather than
as a step somebody remembers.

Subcommands, in the order a project usually needs them:

    cached   DOI                    file the full text 02_literature already downloaded
    grep     DOI|path  PATTERN      contexts around a regex in a filed paper
    bioproject PRJNAxxxxxx [...]      the submitters' own description, GEO id and PMIDs
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

import collections, csv, json, os, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import BIOPROJECT_XML, ENA, GEO, MODULE, STUDIES, TEXT_CACHE

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


def cache_text(doi: str) -> str | None:
    """Full text for a DOI from 02_literature's own cache, or None.

    The cache holds 711 full texts already. Before 2026-10-08 `resolve()` went straight to
    Europe PMC and passed a DOI where the endpoint expects a PMCID, so every DOI lookup
    raised an HTTPError and the cache was never consulted. Reading a paper cost a network
    round trip for a file that was on disk.
    """
    if not TEXT_CACHE.exists():
        return None
    want = doi.strip().lower()
    with TEXT_CACHE.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            key, _, blob = line.partition("\t")
            if key.strip().lower() != want:
                continue
            try:
                d = json.loads(blob)
            except Exception:
                return None
            if isinstance(d, str):
                return d
            for k in ("fulltext", "full_text", "text", "body", "xml"):
                if isinstance(d.get(k), str) and d[k].strip():
                    return d[k]
            return json.dumps(d)
    return None


def resolve(arg: str) -> tuple[str, Path | None]:
    """Accept a DOI, a PMCID or a path. Returns (text, path).

    Order matters: local first, network last. A DOI goes to the 02_literature text cache,
    which is where the papers for this project actually live.
    """
    p = Path(arg)
    if p.exists():
        return p.read_text(encoding="utf-8", errors="replace"), p
    if arg.startswith("10."):
        t = cache_text(arg)
        if t:
            return t, None
        hits = sorted(STUDIES.rglob(f"*{arg.replace('/', '_')}*"))
        hits = [h for h in hits if h.is_file()]
        if hits:
            return hits[0].read_text(encoding="utf-8", errors="replace"), hits[0]
        sys.exit(f"REFUSED: {arg} is not in the text cache and not under studies/. "
                 f"Fetch it with `sources.py fulltext <PMCID>` first.")
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
    print(f"{doi}  {pmcid}  {len(xml):,} bytes -> {f.relative_to(MODULE)}")
    for t in re.findall(r"<title>(.*?)</title>", xml, re.S):
        t = re.sub(r"<[^>]+>", "", t).strip()
        if t:
            print(f"  SEC: {t}")


def cmd_sections(arg: str, pattern: str, cap=2500) -> None:
    cap = int(cap)          # argv strings
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
    cache = GEO / f"{gse}_brief.txt"
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
        print(f"cached -> {cache.relative_to(MODULE)}")
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
            elif k in ("growth_protocol_ch1", "extract_protocol_ch1"):
                # Where SETTING lives for a GEO deposit: the per-sample growth protocol is
                # routinely the only statement that the plants were in a field at all.
                d.setdefault(k, v[:300])
        recs.append(d)
    print(f"{gse}: {len(recs)} samples")
    for k in sorted({k for r in recs for k in r} - {"GSM"}):
        c = collections.Counter(r.get(k, "") for r in recs)
        if len(c) <= 6:
            print(f"  {k:<22} " + "  ".join(f"{a or '(blank)'}:{n}" for a, n in c.most_common()))
        else:
            top = "  ".join(f"{a}:{n}" for a, n in c.most_common(3))
            print(f"  {k:<22} {len(c)} distinct, top: {top}")


def _cache_index() -> dict:
    """{doi -> record} from 02_literature's text cache. Lines are `doi<TAB>json`."""
    out = {}
    if not TEXT_CACHE.exists():
        return out
    with TEXT_CACHE.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "\t" not in line:
                continue
            doi, js = line.split("\t", 1)
            try:
                out[doi.strip().lower()] = json.loads(js)
            except ValueError:
                continue
    return out


def cmd_cached(doi: str) -> None:
    """File a paper's full text from 02_literature's cache instead of refetching it.

    `02_literature/02_text/data/text_cache.jsonl` holds 711 full texts keyed by DOI, built
    while mining supplements for the Kraken cohort. It covers 80 of the 81 STAT-frame cereal
    projects, so for anything that came in through that frame this is the first thing to try:
    no network, no rate limit, and the paper is already on disk.

    The cached text is PLAIN TEXT extracted from a PDF or from Europe PMC, not JATS XML, so
    `sections` cannot split it on <title>. Use `grep` on the filed file instead.
    """
    rec = _cache_index().get(doi.strip().lower())
    if rec is None:
        sys.exit(f"REFUSED: {doi} is not in the cache. Try `sources.py papers` / `fulltext`.")
    txt = (rec.get("full_text") or "").strip()
    if not txt:
        sys.exit(f"REFUSED: {doi} is cached but empty "
                 f"(oa_status={rec.get('oa_status')}, error={rec.get('error')})")
    out = STUDIES / ("doi_" + doi.replace("/", "_"))
    out.mkdir(parents=True, exist_ok=True)
    f = out / f"{doi.split('/')[-1]}_cachedtext.txt"
    f.write_text(txt)
    print(f"{doi}  {len(txt):,} chars  oa={rec.get('oa_status')}  "
          f"src={rec.get('pdf_url') or '?'}  -> {f.relative_to(MODULE)}")
    cmd_scan(str(f))


def cmd_grep(arg: str, pattern: str, width: str = "260", hits: str = "3") -> None:
    """Contexts around a regex in a filed paper. Works on cached plain text and on JATS XML.

    The replacement for `sections` when the text has no markup, and the honest tool for
    checking a single claim: it shows the sentence a curated value will rest on.
    """
    width, hits = int(width), int(hits)
    text, _ = resolve(arg)
    t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", text))
    rx = re.compile(pattern, re.I)
    n = 0
    for m in rx.finditer(t):
        print(f"  ...{t[max(0, m.start() - width):m.start() + width].strip()}\n")
        n += 1
        if n >= hits:
            break
    if not n:
        print(f"  (no match for {pattern!r})")


def cmd_bioproject(*accs: str) -> None:
    """The NCBI BioProject records for one or more accessions, as title/GEO/PMIDs/description.

    The most under-used source in this project. A BioProject `Description` is free text the
    SUBMITTERS wrote about their own experiment, and it routinely contains the methods
    outright - "This dataset comprises 3'RNAseq data from a greenhouse experiment... Plants
    were grown under controlled conditions". Thirty records cost two requests, which makes it
    the cheapest setting evidence available anywhere.

    It also yields the GEO series id and any linked PMIDs, so it resolves papers for projects
    the Europe PMC accession sweep cannot reach.

    Remember the asymmetry: archive metadata may rule a project OUT of the cohort but not IN,
    because a deposit still needs an available manuscript to enter (leon 2026-10-06).
    """
    accs = [a.strip() for a in accs if a.strip()]
    BIOPROJECT_XML.mkdir(parents=True, exist_ok=True)
    key = os.environ.get("NCBI_API_KEY", "")
    kq = f"&api_key={key}" if key else ""
    for i in range(0, len(accs), 30):
        chunk = accs[i:i + 30]
        term = " OR ".join(f"{a}[Project Accession]" for a in chunk)
        u = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject"
             f"&retmax=200&retmode=json&term={urllib.parse.quote(term)}{kq}")
        ids = json.loads(get(u, 90))["esearchresult"]["idlist"]
        if not ids:
            print(f"  (no UIDs for {chunk[0]}..{chunk[-1]})", file=sys.stderr)
            continue
        x = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=bioproject"
                f"&id={','.join(ids)}&retmode=xml{kq}", 180)
        for blk in re.split(r"(?=<Project>)", x)[1:]:
            m = re.search(r'accession="(PRJ[A-Z]{2}\d+)"', blk)
            if not m:
                continue
            acc = m.group(1)
            (BIOPROJECT_XML / f"{acc}.xml").write_text(blk)

            def g(tag: str) -> str:
                mm = re.search(rf"<{tag}>(.*?)</{tag}>", blk, re.S)
                return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", mm.group(1))).strip() if mm else ""

            geo = re.search(r'<CenterID center="GEO"[^>]*>(GSE\d+)</CenterID>', blk)
            pubs = re.findall(r'<Publication[^>]*id="([^"]+)"', blk)
            print(f"=== {acc}   GEO {geo.group(1) if geo else '-'}   PMID {','.join(pubs) or '-'}")
            print(f"    TITLE: {g('Title')[:150]}")
            d = g("Description")
            print(f"    DESC : {d[:520] if d else '(none)'}")
        time.sleep(0.4)


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
    cache = ENA / f"{acc}.tsv"
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(txt)
    print(f"{acc}: {len(rows)} runs -> {cache.relative_to(MODULE)}")
    for k in hdr[1:]:
        c = collections.Counter(r.get(k, "") for r in rows)
        if len(c) <= 6:
            print(f"  {k:<20} " + "  ".join(f"{a or '(blank)'}:{n}" for a, n in c.most_common()))
        else:
            print(f"  {k:<20} {len(c)} distinct, top: "
                  + "  ".join(f"{a}:{n}" for a, n in c.most_common(3)))


def cmd_find(query: str, rows: str = "8") -> None:
    """Europe PMC topic search, for the deposit whose paper nothing links to.

    `papers` searches on the accession STRING and is the right instrument when the paper
    cites its data. It returns nothing for a deposit the authors never cited, or cited only
    in a supplement EPMC did not index - 4 of 8 field candidates on 2026-10-08. This is the
    fallback: search the distinctive nouns from the BioProject title (a gene symbol, a
    cultivar, an institute) and check the hits against the deposit by hand.

    A hit is a CANDIDATE, never a resolution. The paper has to state the accession, or
    describe the same samples, before anything is curated from it.
    """
    q = urllib.parse.quote(query)
    d = json.loads(get(f"{EPMC}/search?query={q}&format=json&pageSize={int(rows)}"
                       f"&resultType=core"))
    print(f"{d['hitCount']} hit(s) for {query!r}\n")
    for r in d["resultList"]["result"]:
        cached = "cached" if TEXT_CACHE.exists() and cache_text(r.get("doi") or "~") else ""
        print(f"  {r.get('doi') or '(no doi)':<36} {r.get('pmcid') or '-':<12} "
              f"OA={r.get('isOpenAccess')} {r.get('pubYear')} {cached}")
        print(f"    {(r.get('title') or '')[:140]}")


CMDS = {"papers": (cmd_papers, 1), "fulltext": (cmd_fulltext, 1), "sections": (cmd_sections, 2),
        "scan": (cmd_scan, 1), "geo": (cmd_geo, 1), "ena": (cmd_ena, 1),
        "bioproject": (cmd_bioproject, 1), "cached": (cmd_cached, 1), "grep": (cmd_grep, 2), "find": (cmd_find, 1)}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in CMDS:
        sys.exit(__doc__)
    fn, n = CMDS[sys.argv[1]]
    args = sys.argv[2:]
    if len(args) < n:
        sys.exit(f"usage: sources.py {sys.argv[1]} " + " ".join(["<arg>"] * n))
    fn(*args)
