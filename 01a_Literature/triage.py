#!/usr/bin/env python3
"""triage.py — one row per candidate study, carrying the whole funnel.

Before any manual work, every candidate BioProject should answer two questions without a
human: **how much of it exists**, and **what will it take to resolve per-run metadata**.
This produces both, normalised across four tables.

    papers.tsv              one row per paper
    bioprojects.tsv         one row per BioProject: the funnel and the triage state
    paper_bioproject.tsv    one row per claim, carrying `relation`
    runs.tsv                one row per run: the per-sample table the adapters fill

The link table is not optional. paper and BioProject are many-to-many in both directions here:
11 of 78 papers name more than one project and 8 of 96 projects are named by more than one
paper. Ada21 names ten BioProjects and generated one of them, having reprocessed 486 older
datasets, so a foreign key on either table would let it claim nine datasets belonging to the
papers that actually produced them. `relation` defaults to `cited` and is promoted to
`generated` or `reanalysed` only on words in the text.

The funnel, widest to narrowest:

    sra_runs        runs in SRA (NCBI esearch). The external truth.
    sra_rnaseq      of those, RNA-Seq strategy.
    screened        runs present in stat_cache.jsonl: our SRA search found AND profiled them.
    gate_passed     runs in runs.tsv: cleared the STAT kingdom thresholds.
    runs_unexamined sra_runs - gate_passed. Reads that exist and were never looked at.

Measured 2026-10-05, the gap that matters is `screened` to `gate_passed`, not `sra_runs` to
`screened`: the SRA search found 1,954 of PRJNA383416's 1,960 runs and the gate kept 1. A
study contributing one sample to the cohort is usually a large study the gate truncated, not
a small study. Without these columns a registry row reads as "1 sample" and the 1,959
unexamined runs are invisible.

The triage state says what resolving per-run metadata will cost:

    sra-complete    per-run BioSamples carrying geo AND date. Nothing to do.
    sra-partial     per-run BioSamples, but geo or date missing. A supplement may fill it.
    one-biosample   every run points at a single BioSample, so SRA has no per-run metadata
                    at all and a supplement is mandatory. PRJNA383416 is the type specimen:
                    1,960 runs, one BioSample titled "RNA isolated from multiple tissues and
                    conditions". Its per-run detail was recoverable only by joining the
                    paper's Supplementary Table 1 on SRA's LibraryName.
    few-biosamples  fewer BioSamples than runs, so some pooling or the same pattern partially.
    no-runs         the accession resolves to nothing, or to something that is not plant reads.

Runinfo is cached per accession under `data/runinfo/`, so reruns cost nothing and the
expensive pass happens once.

Usage:

    python 01a_Literature/triage.py --from-undermind       # scrape accessions from the .md files
    python 01a_Literature/triage.py --accessions PRJNA383416 PRJEB39201
    python 01a_Literature/triage.py --from-undermind --refresh    # ignore the cache
"""
import argparse, collections, csv, io, json, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _paths import ROOT

STAT_CACHE = ROOT / "01_stat/02_fetch/data/stat_cache.jsonl"
RUNS       = ROOT / "01_stat/03_filter/data/runs.tsv"
CACHE      = HERE / "data/runinfo"
BS_CACHE   = HERE / "data/biosample_attrs"
BS_BATCH   = 300    # efetch accepts far more, but a failed batch is cheaper to retry small
BS_SAMPLE  = 20     # BioSamples sampled per project to estimate field coverage.
                    # Coverage is near-bimodal per project (measured 2026-10-05: projects sit
                    # at 0% or 100% far more often than between), so 20 settles the triage
                    # state as reliably as 60 at a third of the cost. Studies that survive
                    # screening get a full fetch via --full-biosamples; over-sampling here is
                    # work paid for twice.
EUTILS     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA         = {"User-Agent": "crypt/triage (leon.lenzo@curtin.edu.au)"}

# efetch rettype=runinfo omits the header row when driven from WebEnv, so it is supplied here.
# Changing this list silently misaligns every column; verify against a known run if NCBI ever
# alters the runinfo schema.
RUNINFO_HEADER = (
    "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,"
    "download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,"
    "LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,"
    "ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,"
    "source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,"
    "Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,"
    "RunHash,ReadHash")

_BP_RE  = re.compile(r'"BioProject":"([^"]+)"')
_ACC_RE = re.compile(r"\bPRJ[DEN][A-Z]\d+\b")
_REF_RE = re.compile(r"\\\[([A-Z][a-z]{2}\d{2}[a-z]?)\\\]")
_HEAD   = 220   # BioProject sits in the first ~200 bytes of a stat_cache line

# Values NCBI stores to mean "no answer". Counting these as populated metadata overstates
# coverage badly; the same list is applied in 02_literature.
_PLACEHOLDER = {"missing", "not applicable", "not collected", "na", "n/a", "none", "unknown",
                "not provided", "restricted access", "not available", "unspecified", ""}

# Undermind's decision strings are free text ("Meets, field subset only", "Does not meet as an
# independent field collection"). Normalise to a base verdict for filtering, and keep the
# original verbatim, because the qualifier is often the useful part.
def _base_decision(text: str) -> str:
    t = text.lower()
    if t.startswith("meets"):          return "meets"
    if t.startswith("does not meet"):  return "does_not_meet"
    if t.startswith("unresolved"):     return "unresolved"
    if "not a new collection" in t or "reanalysis" in t or "same dataset" in t:
        return "duplicate_dataset"
    if "same study" in t:              return "duplicate_dataset"
    return "other"


def undermind_rows(md: Path) -> list[tuple[str, list[str]]]:
    """(section title, cells) for every data row of every table in an Undermind report.

    Section titles are what carry provenance: "Accession-list screening batch*" rows came from
    the accessions we supplied out of STAT, "Search-result screening batch*" rows are what the
    search found on its own. The table before the first batch heading has no label and is
    reported as `lead-table` rather than guessed at.
    """
    out, sec, started = [], "lead-table", False
    for i, line in enumerate(md.read_text().splitlines()):
        # the document repeats its own H1 after the table of contents; content starts there
        if line.startswith("# ") and i > 10:
            started = True
            continue
        if not started:
            continue
        if line.startswith("## "):
            sec = line[3:].strip()
            continue
        if not line.startswith("|") or re.match(r"^\|[:\s|-]+\|$", line):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) >= 3 and "Decision" not in cells[1]:
            out.append((sec, cells))
    return out


def undermind_provenance(md_files: list[Path]) -> dict:
    """BioProject -> dict(source, undermind_decision, undermind_decision_raw, undermind_ref).

    `source` is the union over every mention, so a project named both in the accessions we
    supplied out of STAT and in the search's own results is `both`. That distinction measures
    the search's independent recall and must not be collapsed: `accession_list` alone means
    the search did not find it, `undermind_search` alone means our SRA queries never reached it.
    """
    srcs = collections.defaultdict(set)
    dec, refs = {}, {}
    for md in md_files:
        for sec, cells in undermind_rows(md):
            s = sec.lower()
            src = ("accession_list" if "accession-list" in s
                   else "undermind_search" if "search-result" in s
                   else "undermind_report")
            row = " | ".join(cells)
            accs = set(_ACC_RE.findall(row))
            for acc in accs:
                srcs[acc].add(src)
            # The first cell names the row's subject. A decision attaches to that subject, so
            # only take it when the subject is this accession, or when the subject is a paper
            # reference (then the decision covers every accession the row names).
            subject_acc = set(_ACC_RE.findall(cells[0]))
            ref = _REF_RE.search(cells[0])
            if len(cells) > 1:
                for acc in (subject_acc or (accs if ref else set())):
                    dec[acc] = cells[1]
                    if ref:
                        refs.setdefault(acc, ref.group(1))

    def collapse(names: set) -> str:
        real = names - {"undermind_report"}
        if len(real) > 1:
            return "both"
        return (real or names).pop()

    return {a: dict(source=collapse(srcs[a]),
                    undermind_decision=_base_decision(dec.get(a, "")),
                    undermind_decision_raw=dec.get(a, ""),
                    undermind_ref=refs.get(a, ""))
            for a in srcs}


def _post(endpoint: str, **params) -> str:
    params.setdefault("api_key", os.environ.get("NCBI_API_KEY", ""))
    body = urllib.parse.urlencode({k: v for k, v in params.items() if v}).encode()
    req = urllib.request.Request(EUTILS + endpoint, data=body, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read().decode("utf-8", "replace")


_DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\)\]]+")

# What a row has to say for a paper's claim on a BioProject to be more than "mentioned it".
# Default is `cited`, never `generated`: Ada21 names ten BioProjects and generated one, so
# guessing `generated` from a mention would double-count nine datasets against the papers that
# actually produced them. Promotion needs words on the page.
_REANALYSIS = ("not a new collection", "reanalys", "reprocess", "reused", "same dataset",
               "already counted", "same field dataset")
_GENERATED  = ("generated", "deposited", "new field", "this study produced", "raw reads reported")


def undermind_references(md_files: list[Path]) -> dict:
    """Undermind citation key -> dict(doi, citation), from the References section."""
    out = {}
    for md in md_files:
        for line in md.read_text().splitlines():
            m = _REF_RE.match(line.strip())
            if not m:
                continue
            ref, rest = m.group(1), line.strip()[m.end():].strip()
            doi = _DOI_RE.search(rest)
            out.setdefault(ref, dict(doi=doi.group(0).rstrip(').,') if doi else "",
                                     citation=re.sub(r"\s+", " ", rest)[:400]))
    return out


def paper_links(md_files: list[Path]) -> list[dict]:
    """One row per (paper, BioProject) claim, with the relation and the text it rests on.

    paper <-> BioProject is many-to-many in both directions in this corpus: 11 of 78 papers
    name more than one project, and 8 of 96 projects are named by more than one paper. So the
    link cannot be a column on either table without losing claims or inventing them.
    """
    links, seen = [], set()
    for md in md_files:
        for sec, cells in undermind_rows(md):
            row = " | ".join(cells)
            low = row.lower()
            refs = set(_REF_RE.findall(row))
            accs = set(_ACC_RE.findall(row))
            if not refs or not accs:
                continue
            relation = ("reanalysed" if any(k in low for k in _REANALYSIS)
                        else "generated" if any(k in low for k in _GENERATED)
                        else "cited")
            for ref in sorted(refs):
                for acc in sorted(accs):
                    if (ref, acc) in seen:
                        continue
                    seen.add((ref, acc))
                    links.append(dict(paper_ref=ref, BioProject=acc, relation=relation,
                                      section=sec, evidence=re.sub(r"\s+", " ", row)[:300]))
    return links


def runinfo(acc: str, refresh: bool = False) -> list[dict]:
    """Every run in a BioProject, from cache when possible."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{acc}.json"
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    es = json.loads(_post("esearch.fcgi", db="sra", term=f"{acc}[BioProject]",
                          retmode="json", usehistory="y"))["esearchresult"]
    if es.get("count") in (None, "0"):
        path.write_text("[]")
        return []
    time.sleep(0.12)
    body = _post("efetch.fcgi", db="sra", WebEnv=es["webenv"], query_key=es["querykey"],
                 rettype="runinfo", retmode="csv")
    rows = [r for r in csv.DictReader(io.StringIO(RUNINFO_HEADER + "\n" + body))
            if (r.get("Run") or "").startswith(("SRR", "ERR", "DRR"))]
    path.write_text(json.dumps(rows))
    return rows


def local_counts() -> tuple[collections.Counter, collections.Counter]:
    """(runs screened by STAT, runs that passed the gate), per BioProject."""
    screened = collections.Counter()
    if STAT_CACHE.exists():
        with open(STAT_CACHE) as fh:
            for line in fh:
                m = _BP_RE.search(line[:_HEAD])
                if m:
                    screened[m.group(1)] += 1
    else:
        print(f"warning: {STAT_CACHE} missing; screened counts will be 0", file=sys.stderr)
    passed = collections.Counter()
    with open(RUNS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            passed[r["BioProject"]] += 1
    return screened, passed


# ENA/DDBJ-submitted BioSamples (SAME*/SAMD*) return EMPTY from NCBI efetch, silently, with
# no error. Every PRJEB project in this corpus is affected: 16 projects and 3,191 runs read as
# "no metadata" until routed here instead. The same trap is handled in
# 02_literature/01_search/meta_search.py; this mirrors its attribute mapping.
#
# ENA's own portal API (filereport) would answer in one request per project rather than one per
# sample, but it returned HTTP 500 on every field list on 2026-10-05. Prefer it if it recovers.
_EBI_BASE     = "https://www.ebi.ac.uk/biosamples/samples"
_EBI_PREFIXES = ("SAME", "SAMD")
_EBI_ATTR_MAP = {
    "geographic location (country and/or sea)": "geo_loc_name",
    "collection date":                          "collection_date",
    "organism part":                            "tissue",
    "host scientific name":                     "host",
    "description":                              "bs_description",
}
_EBI_LOCALITY = "geographic location (region and locality)"


# ENA's portal answers for a WHOLE PROJECT in one request, and returns more than the
# per-sample route does. Measured 2026-10-05 on PRJDB7234: 1,557 runs in 3.6 s against an
# estimated 8.5 minutes for 1,017 sequential EBI BioSamples calls, and the portal additionally
# carries `country` ("Japan:Osaka, Takatsuki") which the per-sample route never returned.
#
# It was returning HTTP 500 on every field list earlier the same day, which is why the slow
# path was built at all. Treat portal failure as expected and fall back rather than trusting it.
_ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
_ENA_PORTAL_FIELDS = ("run_accession,sample_accession,collection_date,country,location,"
                      "tissue_type,cultivar,description,sample_title,host,strain,isolate,"
                      "dev_stage,sample_alias")
# portal field -> the name used everywhere else here
_ENA_PORTAL_MAP = {"collection_date": "collection_date", "country": "geo_loc_name",
                   "location": "lat_lon", "tissue_type": "tissue", "cultivar": "cultivar",
                   "host": "host", "strain": "strain", "isolate": "isolate",
                   "dev_stage": "dev_stage", "sample_title": "bs_description"}


def ena_portal_attrs(bp: str) -> dict:
    """BioSample accession -> attributes, for a whole ENA/DDBJ project in one request.

    Returns {} on any failure so the caller falls back to the per-sample route. Placeholder
    values are dropped here as everywhere else.
    """
    import io
    q = urllib.parse.urlencode({"accession": bp, "result": "read_run",
                                "fields": _ENA_PORTAL_FIELDS, "format": "tsv"})
    try:
        req = urllib.request.Request(f"{_ENA_PORTAL}?{q}", headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            body = r.read().decode("utf-8", "replace")
        if not body.startswith("run_accession"):
            return {}
    except Exception as exc:
        print(f"    ENA portal unavailable for {bp} ({type(exc).__name__}); "
              f"falling back to per-sample", file=sys.stderr)
        return {}
    out: dict = {}
    for row in csv.DictReader(io.StringIO(body), delimiter="\t"):
        acc = row.get("sample_accession", "")
        if not acc:
            continue
        d = out.setdefault(acc, {})
        for src, dst in _ENA_PORTAL_MAP.items():
            v = (row.get(src) or "").strip()
            if v and v.lower() not in _PLACEHOLDER and dst not in d:
                d[dst] = v
    return out


def _ebi_biosample(acc: str) -> dict:
    req = urllib.request.Request(f"{_EBI_BASE}/{acc}",
                                 headers={**UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            chars = json.loads(r.read()).get("characteristics", {})
    except Exception:
        return {}
    out = {}
    for ena_key, our_key in _EBI_ATTR_MAP.items():
        items = chars.get(ena_key) or []
        val = (items[0].get("text") or "").strip() if items else ""
        if val and val.lower() not in _PLACEHOLDER:
            out[our_key] = val
    loc = chars.get(_EBI_LOCALITY) or []
    locality = (loc[0].get("text") or "").strip() if loc else ""
    if locality and locality.lower() not in _PLACEHOLDER:
        out["geo_loc_name"] = (f"{out['geo_loc_name']}: {locality}"
                               if out.get("geo_loc_name") else locality)
    return out


def biosample_attrs(accs: list[str], refresh: bool = False,
                    project: str | None = None) -> dict:
    """BioSample accession -> attribute dict, placeholders dropped, cached per BioProject batch.

    Values NCBI uses to mean "no answer" are removed here rather than downstream, because a
    populated-looking `geo_loc_name` of "missing" inflates every coverage number computed from
    this table. Note the converse is not handled and cannot be: "USA: California, Davis" is a
    real string naming an institution, not a collection site, and only the study design
    distinguishes the two.
    """
    out = {}
    ebi = [a for a in accs if a[:4] in _EBI_PREFIXES]
    ncbi = [a for a in accs if a[:4] not in _EBI_PREFIXES]
    if ebi and project:
        got = ena_portal_attrs(project)
        if got:
            hit = {a: got[a] for a in ebi if a in got}
            print(f"    {len(hit)} of {len(ebi)} ENA/DDBJ BioSamples from the portal "
                  f"in ONE request", file=sys.stderr)
            out.update(hit)
            ebi = [a for a in ebi if a not in hit]
    if ebi:
        # EBI BioSamples serves one sample per request and EBI shares a 2 req/s limiter.
        print(f"    {len(ebi)} ENA/DDBJ BioSamples via EBI", file=sys.stderr)
        for a in ebi:
            out[a] = _ebi_biosample(a)
            time.sleep(0.5)
    for i in range(0, len(ncbi), BS_BATCH):
        batch = ncbi[i: i + BS_BATCH]
        try:
            raw = _post("efetch.fcgi", db="biosample", id=",".join(batch),
                        rettype="xml", retmode="xml")
            root = ET.fromstring(raw)
        except Exception as exc:
            print(f"    biosample batch error: {type(exc).__name__}: {exc}", file=sys.stderr)
            time.sleep(0.5)
            continue
        for bs in root.findall(".//BioSample"):
            d = {}
            for a in bs.iter("Attribute"):
                k = (a.get("harmonized_name") or a.get("attribute_name") or "").lower().strip()
                v = (a.text or "").strip()
                if k and v.lower() not in _PLACEHOLDER:
                    d[k] = v
            out[bs.get("accession", "")] = d
        time.sleep(0.15)
    return out


def load_attr_cache(path: Path, bs_all: list[str]) -> tuple[dict, bool, int]:
    """(attrs, complete, sampled) from a cache file, tolerating the old bare-dict format.

    The cache records its own scope because it is usually a SAMPLE, not the whole project.
    Without that, 44 of 181 caches held 60 of up to 1,017 BioSamples while looking complete,
    and every run outside the sample would have been written to runs.tsv with blank metadata
    indistinguishable from metadata that genuinely does not exist.
    """
    if not path.exists():
        return {}, False, 0
    obj = json.loads(path.read_text())
    if isinstance(obj, dict) and isinstance(obj.get("attrs"), dict):
        return obj["attrs"], bool(obj.get("complete")), int(obj.get("sampled", 0))
    # legacy: a bare {accession: attrs} map written before the scope was recorded
    return obj, len(obj) >= len(bs_all) and bool(obj), len(obj)


def save_attr_cache(path: Path, attrs: dict, bs_all: list[str]) -> None:
    path.write_text(json.dumps(dict(sampled=len(attrs),
                                    total=len(bs_all),
                                    complete=len(attrs) >= len(bs_all),
                                    attrs=attrs)))


def triage_state(rows: list[dict], attrs: dict | None = None) -> tuple[str, str, str]:
    """(state, geo coverage %, date coverage %) for one BioProject.

    Cardinality first, because it decides whether per-run metadata can exist at all. Only when
    there is a BioSample per run does field completeness become the question, which is why
    `sra-complete` and `sra-partial` are reachable only from `per-run-biosamples`.
    """
    if not rows:
        return "no-runs", "", ""
    n_runs = len(rows)
    bs = {r["BioSample"] for r in rows if r.get("BioSample")}
    if len(bs) <= 1:
        return "one-biosample", "", ""
    if len(bs) < n_runs:
        return "few-biosamples", "", ""
    if not attrs:
        return "per-run-biosamples", "", ""
    seen = [attrs[a] for a in bs if a in attrs]
    if not seen:
        return "per-run-biosamples", "", ""
    geo = sum(1 for d in seen if d.get("geo_loc_name"))
    date = sum(1 for d in seen if d.get("collection_date"))
    both = sum(1 for d in seen if d.get("geo_loc_name") and d.get("collection_date"))
    g, t = round(100 * geo / len(seen)), round(100 * date / len(seen))
    # A study is only "complete" when nearly every sample carries both. The threshold is 90%
    # rather than 100% because a handful of blank rows in an otherwise documented study is a
    # submission slip, not a reason to queue it for manual work.
    state = "sra-complete" if 100 * both / len(seen) >= 90 else "sra-partial"
    return state, g, t


COHORT   = HERE / "data/kraken_cohort_studies.tsv"
RESOLVED  = HERE / "data/resolved_accessions.tsv"
RUN_SCOPE = HERE / "data/run_scope.tsv"
MENTIONS  = HERE / "data/accession_papers.tsv"
FOUND     = HERE / "data/found.tsv"


def found_accessions() -> list[dict]:
    """BioProjects recorded in found.tsv during review, as links.

    This is the feedback loop that makes found.tsv worth keeping. An accession read out of a
    paper's data statement is a fact the pipeline can act on, but only if something reads it:
    PRJNA609211 and PRJNA825139 sat in bioprojects.tsv for an hour with no rows in runs.tsv,
    because the full rebuild only knew the accessions Undermind and the cohort had supplied.

    relation is `cited`: the paper named the accession, which is stronger than a full-text
    `mentions` hit and weaker than `generated`. Who produced the data stays a human call.
    """
    if not FOUND.exists():
        return []
    out, seen = [], set()
    with open(FOUND) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("kind") != "bioproject" or not r.get("value"):
                continue
            key = (r["paper_key"], r["value"])
            if key in seen:
                continue
            seen.add(key)
            out.append(dict(paper_ref=r["paper_key"], BioProject=r["value"],
                            relation="cited", section="found_during_review",
                            evidence=f"recorded in found.tsv by {r.get('found_by','')}: "
                                     f"{r.get('note','')[:200]}"))
    return out



def mention_links() -> list[dict]:
    """Links from accession_papers.py: papers whose full text names an accession.

    relation is `mentions`, which is weaker than every other relation here and deliberately
    so. A Europe PMC full-text hit proves the accession appears in the paper and nothing
    more: of the eleven papers naming PRJNA306542, five contributed no data to it and were
    reanalysing. Promoting a mention to `generated` is the error the relation column exists
    to prevent.

    Without these the link table knows only what Undermind and our own cohort told it, so
    PRJNA306542 reads as a one-paper project when it is in fact six experiments from six
    papers, and `studies/by-project/` built from it shows one symlink instead of eleven.
    """
    if not MENTIONS.exists():
        return []
    out, seen = [], set()
    with open(MENTIONS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = ((r.get("doi") or "").strip().lower(), r["accession"])
            if not all(key) or key in seen:
                continue
            seen.add(key)
            out.append(dict(paper_ref=key[0], BioProject=r["accession"],
                            relation="mentions", section="europepmc_fulltext",
                            evidence=f"Europe PMC full text names {r['accession']}; "
                                     f"{r.get('year','')} {r.get('title','')[:140]}"))
    return out


_RANGE = re.compile(r"^([A-Z]+)(\d+)\s*-\s*([A-Z]+)?(\d+)$", re.I)


def expand_spec(spec: str, rows: list[dict]) -> set[str]:
    """Run accessions a scope spec covers, resolved against this project's actual runs.

    A paper often deposits into a shared umbrella BioProject and accounts for only part of it.
    Wei et al. 2016 names 12 RNA-seq pools inside PRJNA306542, which holds 433 runs from the
    same centre; without a scope the link claims all 433 and overstates that paper's
    contribution 36-fold. Any count built on the link table inherits that error.

    Grammar, deliberately small:
        all                      every run in the project (the default when unspecified)
        SRX1521275-SRX1521286    inclusive range over Experiment or Run accessions
        SRR1,SRR2,...            an explicit list
    Ranges expand against the accessions that EXIST rather than by arithmetic, so a gap in the
    numbering cannot invent a run that was never deposited.
    """
    spec = (spec or "all").strip()
    if spec.lower() in ("", "all", "*"):
        return {r["Run"] for r in rows}
    out: set[str] = set()
    for part in (p.strip() for p in spec.split(",") if p.strip()):
        m = _RANGE.match(part)
        if m:
            pre, lo, _, hi = m.group(1).upper(), int(m.group(2)), m.group(3), int(m.group(4))
            for r in rows:
                for col in ("Experiment", "Run"):
                    v = r.get(col, "")
                    mm = re.match(r"^([A-Z]+)(\d+)$", v, re.I)
                    if mm and mm.group(1).upper() == pre and lo <= int(mm.group(2)) <= hi:
                        out.add(r["Run"])
            continue
        for r in rows:
            if part.upper() in (r.get("Run", "").upper(), r.get("Experiment", "").upper()):
                out.add(r["Run"])
    return out


def run_scopes() -> dict:
    """(paper_key, BioProject) -> dict(spec, note, evidence). Curated, one row per claim."""
    if not RUN_SCOPE.exists():
        return {}
    with open(RUN_SCOPE) as fh:
        return {(r["paper_key"], r["BioProject"]): r
                for r in csv.DictReader(fh, delimiter="\t") if r.get("paper_key")}



def resolved_links() -> list[dict]:
    """Links recovered by resolve_accessions.py from GEO/SRA-study/run accessions.

    Written as `sra_linked`, never `generated`: an archive lookup connected the paper to the
    BioProject, which is not a claim about who produced the data. The accession the link came
    through is kept in `evidence` so the chain is auditable, which matters because a single
    cited run is weak evidence: Kim22c's five example runs resolve to five different projects.
    """
    if not RESOLVED.exists():
        return []
    with open(RESOLVED) as fh:
        return [dict(paper_ref=r["paper_ref"], BioProject=r["BioProject"],
                     relation="sra_linked", section="resolved_accession",
                     evidence=f"{r['via']} ({r['via_kind']}) -> {r['BioProject']} via SRA runinfo")
                for r in csv.DictReader(fh, delimiter="\t")
                if r["status"] == "resolved" and r["BioProject"]]



def _norm_doi(d: str) -> str:
    return (d or "").strip().lower().rstrip(").,;")


def cohort_studies() -> list[dict]:
    """Papers already in the Kraken cohort: one row per BioProject, carrying its resolved DOI.

    These are candidates for the new frame too, and they must be IN the registry rather than
    alongside it, or "what is the candidate universe" has two answers. They are candidates and
    not members: the cohort was selected pathogen-first, so every one still has to be screened
    against the field/aerial criteria. Their `undermind_decision` stays empty to say so.
    """
    if not COHORT.exists():
        print(f"warning: {COHORT.name} missing; cohort papers will be absent", file=sys.stderr)
        return []
    with open(COHORT) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_papers(refs: dict, links: list[dict], cohort: list[dict],
                 bp_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """(papers, links) keyed on DOI, merging the Undermind corpus and the Kraken cohort.

    Keyed on DOI, not on the Undermind reference. The reference only exists for papers the
    search returned, so keying on it silently drops every cohort paper, and a paper found by
    both routes would appear twice. A paper with no DOI keeps its reference as the key, which
    is the only identifier it has.
    """
    bps = {b["BioProject"]: b for b in bp_rows}
    papers: dict[str, dict] = {}
    extra_links: list[dict] = []

    def row(key: str) -> dict:
        return papers.setdefault(key, dict(
            paper_key=key, doi="", paper_ref="", title="", sources=set(),
            bioprojects=set(), undermind_decision=set(),
            pdf="", supplement="", library_prep="", sampling_design="", status=""))

    bp_by_ref = collections.defaultdict(set)
    for l in links:
        bp_by_ref[l["paper_ref"]].add(l["BioProject"])

    for ref, meta in refs.items():
        key = _norm_doi(meta.get("doi")) or ref
        r = row(key)
        r["doi"] = r["doi"] or _norm_doi(meta.get("doi"))
        r["paper_ref"] = ref
        r["title"] = r["title"] or meta.get("citation", "")
        r["bioprojects"] |= bp_by_ref.get(ref, set())
        for acc in bp_by_ref.get(ref, set()):
            src = bps.get(acc, {}).get("accession_source", "")
            if src:
                r["sources"].add(src if src != "both" else "accession_list")
                if src == "both":
                    r["sources"].add("undermind_search")
            d = bps.get(acc, {}).get("undermind_decision", "")
            if d:
                r["undermind_decision"].add(d)
        if not r["sources"]:
            r["sources"].add("undermind_report")

    # Papers found by Europe PMC full-text search enter as rows too, or they have no worklist
    # entry, no directory, and no place in studies/by-project/ — which showed ONE symlink for
    # a project cited by eleven papers. They are candidates like any other: a mention is not
    # evidence of scope, and `undermind_decision` stays empty to say nobody has screened them.
    for l in links:
        if l.get("relation") != "mentions":
            continue
        key = _norm_doi(l["paper_ref"])
        if not key:
            continue
        r = row(key)
        r["doi"] = r["doi"] or key
        r["sources"].add("europepmc_mention")
        r["bioprojects"].add(l["BioProject"])
        if not r["title"]:
            ev = l.get("evidence", "")
            r["title"] = ev.split("; ", 1)[1] if "; " in ev else ev

    for c in cohort:
        key = _norm_doi(c.get("doi"))
        if not key:
            continue
        r = row(key)
        r["doi"] = r["doi"] or key
        r["title"] = r["title"] or c.get("title", "")
        r["sources"].add("kraken_cohort")
        acc = c.get("BioProject", "")
        if acc:
            r["bioprojects"].add(acc)
            # The cohort link came from NCBI's literature resolution for that BioProject, not
            # from anyone reading the paper. Name it for what it is rather than calling it
            # `generated`, which is a claim about who produced the data.
            extra_links.append(dict(paper_ref=key, BioProject=acc, relation="sra_linked",
                                    section="kraken_cohort",
                                    evidence=f"resolved by 02_literature/01_search for {acc}"))

    out = []
    for r in papers.values():
        out.append(dict(r, sources=";".join(sorted(r["sources"])),
                        n_bioprojects=len(r["bioprojects"]),
                        bioprojects=";".join(sorted(r["bioprojects"])),
                        undermind_decision=";".join(sorted(r["undermind_decision"]))))
    out.sort(key=lambda r: (-r["n_bioprojects"], r["paper_key"]))
    return out, extra_links


def _write(path: Path, recs: list[dict], fields: list[str], key: str | None = None) -> None:
    """Write a table, merging over any existing rows when `key` is given.

    Merging is not a nicety. A scoped run (`--accessions PRJNA746402`) computes rows for one
    project, and a plain overwrite then replaces a 327-project table with a 1-project table.
    That happened on 2026-10-05 and silently destroyed the registry; the tables had to be
    rebuilt. A scoped run must only ever update the rows it actually recomputed.
    """
    merged = recs
    if key and path.exists():
        with open(path) as fh:
            old = {r[key]: r for r in csv.DictReader(fh, delimiter="\t") if r.get(key)}
        n_before = len(old)
        for r in recs:
            old[r[key]] = r
        merged = list(old.values())
        if len(merged) > len(recs):
            print(f"  merged {len(recs)} recomputed row(s) into {n_before} existing")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(merged)
    print(f"  wrote {path.relative_to(ROOT)}  ({len(merged)} rows)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", nargs="*", default=[], metavar="PRJ")
    ap.add_argument("--from-undermind", action="store_true",
                    help="scrape accessions and provenance from the Undermind .md files")
    ap.add_argument("--refresh", action="store_true", help="ignore the runinfo cache")
    ap.add_argument("--no-cohort", action="store_true",
                    help="exclude the Kraken cohort studies; Undermind candidates only")
    ap.add_argument("--full-biosamples", action="store_true",
                    help="fetch every BioSample, not a sample of BS_SAMPLE. Use once a study "
                         "has survived screening and its per-run metadata is actually needed")
    ap.add_argument("--full-only", nargs="*", metavar="PRJ",
                    help="restrict --full-biosamples to these BioProjects")
    ap.add_argument("--no-biosamples", action="store_true",
                    help="skip the BioSample fetch, so no sra-complete/sra-partial split")
    args = ap.parse_args()

    mds = [m for m in sorted(HERE.glob("*.md")) if m.name != "README.md"]
    prov, refs, links = {}, {}, []
    accs = list(dict.fromkeys(args.accessions))
    if args.from_undermind:
        prov = undermind_provenance(mds)
        refs = undermind_references(mds)
        links = paper_links(mds) + resolved_links() + mention_links() + found_accessions()
        accs = list(dict.fromkeys(accs + [l["BioProject"] for l in found_accessions()]))
        accs = list(dict.fromkeys(accs + sorted(prov)
                                  + [l["BioProject"] for l in resolved_links()]))
    cohort = [] if args.no_cohort else cohort_studies()
    if cohort:
        # Cohort accessions join the SAME funnel. Their `screened`/`gate_passed` are already
        # known locally, so the only new cost is runinfo, and it is cached.
        accs = list(dict.fromkeys(accs + [c["BioProject"] for c in cohort if c.get("BioProject")]))
    if not accs:
        sys.exit("no accessions; pass --accessions or --from-undermind")
    if not os.environ.get("NCBI_API_KEY"):
        print("warning: no NCBI_API_KEY, rate limit is 2.5 req/s. "
              "`source ~/.profile` first.", file=sys.stderr)

    print("reading local STAT counts ...", file=sys.stderr)
    screened, passed = local_counts()

    BS_CACHE.mkdir(parents=True, exist_ok=True)
    bioprojects, runs = [], []
    runs_by_bp: dict[str, list[dict]] = {}
    for i, acc in enumerate(accs, 1):
        try:
            rows = runinfo(acc, args.refresh)
        except Exception as exc:
            print(f"  [{i}/{len(accs)}] {acc}: {type(exc).__name__}: {exc}", file=sys.stderr)
            time.sleep(0.5)
            continue

        attrs = {}
        if rows and not args.no_biosamples:
            bs_all = sorted({r["BioSample"] for r in rows if r.get("BioSample")})
            cache = BS_CACHE / f"{acc}.json"
            want_all = args.full_biosamples and (not args.full_only
                                                 or acc in set(args.full_only))
            if not args.refresh:
                attrs, complete, _ = load_attr_cache(cache, bs_all)
                # NCBI efetch answers some SAME*/SAMD* queries with records keyed under a
                # DIFFERENT accession, so the cache looks populated but every lookup misses and
                # the project reads as having no metadata. Nine projects were in that state,
                # including PRJEB39201 at 550 runs. Validate by intersection, not by existence.
                if attrs and not (set(attrs) & set(bs_all)):
                    print(f"    stale cache for {acc} (no key overlap); refetching",
                          file=sys.stderr)
                    attrs, complete = {}, False
                    cache.unlink()
                if attrs and want_all and not complete:
                    attrs = {}      # a sampled cache cannot answer a full request
            if not attrs and len(bs_all) > 1:
                take = bs_all if want_all else bs_all[:BS_SAMPLE]
                attrs = biosample_attrs(take, project=acc)
                save_attr_cache(cache, attrs, bs_all)

        state, geo_pct, date_pct = triage_state(rows, attrs)
        s_n, p_n = screened.get(acc, 0), passed.get(acc, 0)
        pv = prov.get(acc, {})
        bioprojects.append(dict(
            BioProject=acc, triage=state,
            sra_runs=len(rows),
            sra_rnaseq=sum(1 for r in rows if r.get("LibraryStrategy") == "RNA-Seq"),
            sra_biosamples=len({r["BioSample"] for r in rows if r.get("BioSample")}),
            screened=s_n, gate_passed=p_n, runs_unexamined=max(0, len(rows) - p_n),
            gate_kept_pct=round(100 * p_n / s_n, 2) if s_n else "",
            geo_pct=geo_pct, date_pct=date_pct,
            organism=collections.Counter(
                r.get("ScientificName", "") for r in rows).most_common(1)[0][0] if rows else "",
            centre=collections.Counter(
                r.get("CenterName", "") for r in rows).most_common(1)[0][0] if rows else "",
            in_cohort="yes" if p_n else "",
            accession_source=pv.get("source", ""),
            undermind_decision=pv.get("undermind_decision", ""),
            undermind_decision_raw=pv.get("undermind_decision_raw", ""),
        ))

        for r in rows:
            bs_acc = r.get("BioSample", "")
            a = attrs.get(bs_acc, {})
            # Three states, not two. "" with source "" means the BioSample was never fetched,
            # so we do not know; "absent" means it was fetched and carries no such field. The
            # sampling cap makes the first case common (490 of PRJEB39201's 550 runs), and
            # collapsing them understates coverage while looking like data.
            fetched = bs_acc in attrs

            def fld(key):
                if not fetched:
                    return "", "not_fetched"
                v = a.get(key, "")
                return v, ("biosample" if v else "absent")

            tis, tis_src = fld("tissue")
            loc, loc_src = fld("geo_loc_name")
            dat, dat_src = fld("collection_date")
            runs.append(dict(
                Run=r["Run"], BioProject=acc, BioSample=r.get("BioSample", ""),
                Experiment=r.get("Experiment", ""), LibraryName=r.get("LibraryName", ""),
                SampleName=r.get("SampleName", ""),
                LibraryStrategy=r.get("LibraryStrategy", ""),
                LibrarySelection=r.get("LibrarySelection", ""), spots=r.get("spots", ""),
                # Resolved fields. Each is filled from BioSample where SRA has it, and from a
                # per-study adapter over the paper's supplement where it does not. The _source
                # column records which, so a rate can always be recomputed on one source alone.
                tissue=tis, tissue_source=tis_src,
                location=loc, location_source=loc_src,
                collection_date=dat, date_source=dat_src,
                setting="", setting_source="",
                # Whether the SAMPLING conditioned on disease. Surveillance cohorts are
                # collected BECAUSE a plant is symptomatic, which makes them the best
                # co-infection test set available and an unusable prevalence denominator.
                # Neither `setting` nor `tissue` distinguishes the two, so it is its own
                # field: a field study can be either.
                sampling_selection="", sampling_selection_source="",
            ))

        runs_by_bp[acc] = rows
        if i % 10 == 0 or i == len(accs):
            print(f"  [{i}/{len(accs)}] {acc}", file=sys.stderr)
        time.sleep(0.12)

    bioprojects.sort(key=lambda r: -r["runs_unexamined"])

    papers, cohort_links = build_papers(refs, links, cohort, bioprojects)
    links = links + cohort_links

    # Attach the run scope to every claim. Unscoped claims say so explicitly rather than
    # leaving the column blank, because "all 433" and "we never checked" must not look alike.
    scopes = run_scopes()
    for l in links:
        sc = scopes.get((l["paper_ref"], l["BioProject"]))
        rows_for = runs_by_bp.get(l["BioProject"], [])
        if sc:
            hits = expand_spec(sc["spec"], rows_for)
            l["run_scope"] = "subset"
            l["run_spec"] = sc["spec"]
            l["n_runs_claimed"] = len(hits)
            l["scope_evidence"] = sc.get("evidence", "")[:300]
            if not hits:
                print(f"  WARNING: scope {sc['spec']!r} for {l['paper_ref']} / "
                      f"{l['BioProject']} matched 0 runs", file=sys.stderr)
        else:
            l["run_scope"] = "all"
            l["run_spec"] = ""
            l["n_runs_claimed"] = len(rows_for)
            l["scope_evidence"] = ""

    # A BioProject reached by the cohort but never named by Undermind has no accession_source.
    # Label it, or "where did this come from" has a blank answer for a third of the table.
    named = set(prov)
    for b in bioprojects:
        if not b["accession_source"]:
            b["accession_source"] = "kraken_cohort" if b["BioProject"] not in named else ""

    # Which paper actually PRODUCED each dataset, as distinct from which papers mention it.
    # Attribution is a counting problem, not tidiness: Ada21 names ten BioProjects and
    # generated one, having reprocessed 486 older datasets, so letting it own all ten
    # double-counts four datasets against Bue17, Hub15, Rad19 and Car26, who made them.
    #
    # Precedence, strongest evidence first. Blank is a real answer and means a human decides;
    # defaulting to the paper that merely reprocessed the data is the error being avoided.
    by_bp = collections.defaultdict(list)
    for l in links:
        by_bp[l["BioProject"]].append(l)
    for b in bioprojects:
        claims = by_bp.get(b["BioProject"], [])
        gen = [l for l in claims if l["relation"] == "generated"]
        sra = [l for l in claims if l["relation"] == "sra_linked"]
        others = [l for l in claims if l["relation"] not in ("sra_linked",)]
        if len(gen) == 1:
            b["primary_paper"], b["primary_source"] = gen[0]["paper_ref"], "relation=generated"
        elif sra:
            # the DOI 02_literature/01_search resolved from NCBI for this BioProject
            b["primary_paper"], b["primary_source"] = sra[0]["paper_ref"], "sra_resolved_doi"
        elif len(others) == 1:
            b["primary_paper"], b["primary_source"] = others[0]["paper_ref"], "sole_claimant"
        else:
            b["primary_paper"], b["primary_source"] = "", "unresolved" if claims else "no_claim"
        b["n_claims"] = len(claims)
        # Two sources naming different primaries is exactly the case meta_search.py recorded,
        # where a bioRxiv survey mining other people's RNA-seq was matched as the primary paper
        # for five unrelated BioProjects. Flag it rather than picking a winner silently.
        cands = {l["paper_ref"] for l in gen} | {l["paper_ref"] for l in sra}
        b["primary_conflict"] = "yes" if len(cands) > 1 else ""

    print()
    _write(HERE / "data/bioprojects.tsv", bioprojects, list(bioprojects[0]), key="BioProject")
    if not args.from_undermind:
        print("  scoped run: papers.tsv and paper_bioproject.tsv left untouched "
              "(they are derived from the full corpus, not from one accession)")
        return
    _write(HERE / "data/papers.tsv", papers,
           ["paper_key", "doi", "paper_ref", "sources", "n_bioprojects", "bioprojects",
            "undermind_decision", "title", "pdf", "supplement", "library_prep",
            "sampling_design", "status"])
    _write(HERE / "data/paper_bioproject.tsv", links,
           ["paper_ref", "BioProject", "relation", "run_scope", "run_spec", "n_runs_claimed",
            "section", "evidence", "scope_evidence"])
    _write(HERE / "data/runs.tsv", runs, list(runs[0]) if runs else ["Run"], key="Run")

    print(f"\n{len(bioprojects)} BioProjects, {len(papers)} papers, {len(links)} claims, "
          f"{len(runs):,} runs")
    for state, n in collections.Counter(r["triage"] for r in bioprojects).most_common():
        rr = sum(r["sra_runs"] for r in bioprojects if r["triage"] == state)
        print(f"  {state:<22}{n:>4} projects{rr:>9,} runs")
    if links:
        print("  relation:", dict(collections.Counter(l["relation"] for l in links)))
    print("  accession source:",
          dict(collections.Counter(r["accession_source"] or "-" for r in bioprojects)))
    print("  primary attribution:",
          dict(collections.Counter(r["primary_source"] for r in bioprojects).most_common()))
    conflicts = [r["BioProject"] for r in bioprojects if r["primary_conflict"]]
    if conflicts:
        print(f"  primary CONFLICT ({len(conflicts)}): {', '.join(conflicts[:8])}")
    print("  paper sources:",
          dict(collections.Counter(r["sources"] for r in papers).most_common(6)))
    nf = sum(1 for r in runs if r["location_source"] == "not_fetched")
    if nf:
        print(f"  runs whose BioSample was never fetched: {nf:,} of {len(runs):,} "
              f"({100*nf/len(runs):.0f}%). Not missing metadata, UNFETCHED metadata: "
              f"rerun with --full-biosamples for the studies you keep.")
    tot_r = sum(r["sra_runs"] for r in bioprojects)
    tot_p = sum(r["gate_passed"] for r in bioprojects)
    print(f"\n  SRA runs {tot_r:,}; gate-passed {tot_p:,}; never examined {tot_r-tot_p:,}")


if __name__ == "__main__":
    main()
