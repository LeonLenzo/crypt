#!/usr/bin/env python3
"""Bucket BioProjects by setting evidence, using only text already on disk.

Written 2026-10-08 for the STAT-frame tail: 69 cereal projects carrying 552 already-
classified runs, nearly none with a known paper. Hunting a paper for each is the expensive
step, so this runs first and does the cheap half.

What it reads: the cached BioProject XML (Name/Title/Description - the Description is
SUBMITTER free text and routinely states the methods outright), plus the SRA study titles
and the BioSample attribute values from the runinfo/biosample caches. No network.

THE ASYMMETRY. Archive metadata may rule a project OUT of the cohort but not IN
(leon 2026-10-06): a deposit still needs an available manuscript to enter. So:

  controlled  -> actionable. A `setting=controlled` rule can be written from this text.
  field       -> NOT actionable. Marks the project as worth a paper hunt, nothing more.
  mixed       -> both in one sentence. Needs the paper; this is the field-phenotype/
                 controlled-sequencing trap, which metadata cannot resolve.
  none        -> no signal either way. Needs the paper.

THE LINKED-SENTENCE RULE. A growth-environment term only counts when a sampling or
sequencing term appears in the SAME sentence. "Plants were evaluated in field trials;
leaf tissue for RNA-seq was taken from greenhouse-grown seedlings" has both environments,
and only the linked one describes what was sequenced. Scoring the document as a bag of
words called that project field. Eight of twenty hand-read projects turned on this.

Guards that exist because they fired: `bright field` / `dark field` / `field of view` are
microscopy, and `in the field of` is rhetoric. All are excluded by lookaround, not by
hoping they are rare.

Run:  python 01a_Literature/retrieve/setting_screen.py PRJEB45007 PRJNA327013
      python 01a_Literature/retrieve/setting_screen.py --file pending.txt --quotes
"""

from __future__ import annotations

import argparse, collections, csv, json, re, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import BIOPROJECT_XML, BIOSAMPLE, RUNINFO, SILVER, STUDIES, TEXT_CACHE

# --- vocabulary -------------------------------------------------------------------------

CONTROLLED = re.compile(r"""
    greenhouse | glass.?house | growth\s+(chamber|cabinet|room) | climate\s+chamber
  | phytotron | controlled\s+(environment|condition|climate) | growth\s+conditions?
  | incubator | photoperiod | \b\d{1,2}\s*(h|hr|hours?)\s*(light|photoperiod|day)
  | hydroponic | potting | \bpots?\b | \bpotted\b | in\s+vitro | axenic
  | petri | culture\s+(medium|media|plate|flask) | \bPDA\b | \bagar\b
  | inoculat | \bspray.?inoculat | detached\s+leaf | seedling\s+assay
  | artificially\s+infect | controlled\s+inoculation
""", re.I | re.X)

# `field` needs the most care; `plot` and `nursery` are weak on their own.
FIELD = re.compile(r"""
    (?<!bright\s)(?<!dark\s)(?<!bright-)(?<!dark-)
    (?<!in\sthe\s)\bfield\b(?!\s+of\s+(view|research|study|science))
  | field.?(grown|collected|sampled|trial|plot|site|experiment|condition)
  | natural(ly)?\s+(infect|occurr|incidence)
  | \bpaddock | \bfarm(er|s|ers)?\b | commercial\s+(field|crop|farm|plantation)
  | experimental\s+station | research\s+station | \brain.?fed\b
  | \bsurvey(ed|s)?\b | \borchard\b | \bvineyard\b
""", re.I | re.X)

# The link requirement: what makes a growth term describe the SEQUENCED material.
SAMPLING = re.compile(r"""
    sampl | collect | harvest | \bRNA\b | sequenc | librar | transcriptom
  | tissue | \bleaf\b | \bleaves\b | flash.?frozen | snap.?frozen
  | liquid\s+nitrogen | extract
""", re.I | re.X)

# For a PAPER body the link term has to be sequencing-grade. `tissue` and `leaf` are
# everywhere in a results section and in every figure legend, so the loose SAMPLING rule
# scored "rice leaf sheath tissue inoculated with conidia (Scale bar, 10 um.)" as evidence
# about the sequenced material. Requiring RNA extraction or sequencing in the sentence
# keeps the screen on the methods, where the answer actually is.
SEQUENCING = re.compile(r"""
    RNA.?seq | RNA\s+sequenc | transcriptome\s+sequenc | total\s+RNA
  | RNA\s+was\s+(extracted|isolated|prepared) | RNA\s+extraction | RNA\s+isolation
  | librar(y|ies)\s+(was|were|prepar|construct) | cDNA\s+librar
  | sequenc(ed|ing)\s+(on|using|was|at|by) | Illumina | NovaSeq | HiSeq | BGI
  | (deposit|submitt)ed | accession\s+(number|code)
""", re.I | re.X)

SENT = re.compile(r"(?<=[.;!?])\s+|\n+")


def sentences(text: str) -> list[str]:
    return [s.strip() for s in SENT.split(text or "") if s.strip()]


def project_text(bp: str) -> list[tuple[str, str]]:
    """(locus, text) for every local source that might state a setting."""
    out = []
    x = BIOPROJECT_XML / f"{bp}.xml"
    if x.exists():
        blob = x.read_text(errors="replace")
        # Only the ProjectDescr block: Submission/Organization carries institute names
        # ("... Research Station") that trip FIELD with no bearing on the experiment.
        m = re.search(r"<ProjectDescr>(.*?)</ProjectDescr>", blob, re.S)
        for tag in ("Name", "Title", "Description"):
            for v in re.findall(rf"<{tag}>(.*?)</{tag}>", m.group(1) if m else "", re.S):
                v = re.sub(r"<[^>]+>", " ", v)
                v = re.sub(r"\s+", " ", v).strip()
                if v:
                    out.append((f"bioproject/{tag}", v))
    ri = RUNINFO / f"{bp}.json"
    if ri.exists():
        try:
            rows = json.loads(ri.read_text())
        except Exception:
            rows = []
        rows = rows if isinstance(rows, list) else rows.get("rows", [])
        seen = set()
        for r in rows[:400]:
            for col in ("Study_Pubmed_id", "StudyTitle", "ExperimentTitle",
                        "LibraryName", "SampleName", "Title"):
                v = str(r.get(col) or "").strip()
                if v and v not in seen and len(v) > 12:
                    seen.add(v)
                    out.append((f"runinfo/{col}", v))
    ba = BIOSAMPLE / f"{bp}.json"
    if ba.exists():
        try:
            d = json.loads(ba.read_text())
        except Exception:
            d = {}
        # Two cache shapes are in the store, and reading the second one as the first is
        # not a no-op. 203 files are a bare {SAMN: {attr: value}} map; 164 are the newer
        # envelope {sampled, total, complete, attrs: {SAMN: {...}}}. Iterating the
        # envelope's values yields ints and the whole attrs dict, so every per-sample
        # record got stringified into ONE pseudo-value. That value carries all of a
        # sample's text at once, which let a controlled term from one sample pair with a
        # field term from another and produced a spurious `mixed`. It also matched
        # regexes by accident, so the bug inflated hits rather than losing them.
        samples = d.get("attrs") if isinstance(d.get("attrs"), dict) else d
        seen = set()
        for attrs in (samples.values() if isinstance(samples, dict) else []):
            if not isinstance(attrs, dict):
                continue
            for k, v in attrs.items():
                v = str(v).strip()
                if v and (k, v) not in seen:
                    seen.add((k, v))
                    out.append((f"biosample/{k}", v))
    return out


def paper_text(bp: str) -> list[tuple[str, str]]:
    """(locus, text) for every cached full text linked to this BioProject.

    Added 2026-10-08. The archive text settles the controlled projects and leaves the rest
    silent, and 30 of the 42 silent ones turned out to have a paper once Europe PMC was
    searched on the accession. Reading 30 papers by hand is the whole cost of the tail, so
    the same linked-sentence rule is applied to the paper body.

    A hit here is actionable in BOTH directions, unlike the archive: a paper IS the
    available manuscript that no-manuscript-no-entry asks for. It is still a screen and
    not a verdict - the matched sentence has to be read before anything is curated, because
    a Methods section describing a companion greenhouse experiment reads identically to one
    describing the sequenced material.
    """
    ap = SILVER / "accession_papers.tsv"
    if not ap.exists() or not TEXT_CACHE.exists():
        return []
    dois = []
    with ap.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("accession") == bp and (r.get("doi") or "").strip():
                dois.append(r["doi"].strip().lower())
    if not dois:
        return []
    want, out, got = set(dois), [], {}
    with TEXT_CACHE.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            key, sep, blob = line.partition("\t")
            k = key.strip().lower()
            if k not in want or not blob.strip():
                continue
            try:
                d = json.loads(blob)
            except Exception:
                continue
            t = d if isinstance(d, str) else next(
                (d[x] for x in ("fulltext", "full_text", "text", "body", "xml")
                 if isinstance(d.get(x), str) and d[x].strip()), "")
            if t:
                got[k] = t
    # The cache stores at most 60,000 characters and 72% of its texts sit at that cap, which
    # in a Wiley or Nature article removes the Methods outright. Prefer a full copy fetched
    # into studies/ wherever one exists, and fall back to the capped text otherwise.
    for k in dois:
        t = got.get(k, "")
        if not t or len(t) >= 59_980:
            hits = sorted(f for f in (STUDIES / ("doi_" + k.replace("/", "_"))).glob("*")
                          if f.is_file() and f.name != "adapter.yaml")
            if hits:
                t = hits[0].read_text(encoding="utf-8", errors="replace")
        if t:
            out.append((f"paper/{k}", re.sub(r"<[^>]+>", " ", t)))
    return out


def classify(bp: str, papers: bool = False) -> dict:
    ctrl, fld, mixed = [], [], []
    srcs = project_text(bp) + (paper_text(bp) if papers else [])
    for locus, text in srcs:
        for s in sentences(text):
            # A BioSample attribute is a value, not prose: the field IS the statement,
            # so the link requirement cannot apply to it.
            if locus.startswith("biosample/"):
                # An attribute is a value, not prose: the field IS the statement.
                linked = True
            elif locus.startswith("paper/"):
                linked = bool(SEQUENCING.search(s))
            else:
                linked = bool(SAMPLING.search(s))
            if not linked:
                continue
            c, f = bool(CONTROLLED.search(s)), bool(FIELD.search(s))
            if c and f:
                mixed.append((locus, s))
            elif c:
                ctrl.append((locus, s))
            elif f:
                fld.append((locus, s))
    if mixed or (ctrl and fld):
        call = "mixed"
    elif ctrl:
        call = "controlled"
    elif fld:
        call = "field"
    else:
        call = "none"
    return dict(bp=bp, call=call, controlled=ctrl, field=fld, mixed=mixed,
                n_sources=len(srcs), n_papers=sum(1 for l, _ in srcs
                                                  if l.startswith("paper/")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", nargs="*")
    ap.add_argument("--file", help="one accession per line")
    ap.add_argument("--quotes", action="store_true", help="print the matched sentences")
    ap.add_argument("--only", help="print only this call (controlled/field/mixed/none)")
    ap.add_argument("--papers", action="store_true",
                    help="also screen the cached full text of linked papers")
    a = ap.parse_args()

    accs = list(a.accessions)
    if a.file:
        accs += [l.strip() for l in Path(a.file).read_text().split() if l.strip()]
    if not accs:
        sys.exit("usage: setting_screen.py <BioProject>... | --file list.txt")

    results = [classify(bp, papers=a.papers) for bp in accs]
    tally = collections.Counter(r["call"] for r in results)
    print("  ".join(f"{k}:{tally[k]}" for k in ("controlled", "field", "mixed", "none")))
    print()
    for r in results:
        if a.only and r["call"] != a.only:
            continue
        print(f"{r['bp']:<14}{r['call']:<11}{r['n_sources']:>4} sources"
              f"{'  ' + str(r['n_papers']) + ' paper(s)' if r.get('n_papers') else ''}")
        if a.quotes:
            for tag, hits in (("CTRL", r["controlled"]), ("FIELD", r["field"]),
                              ("MIXED", r["mixed"])):
                for locus, s in hits[:3]:
                    print(f"    {tag:<6}{locus:<26}{s[:170]}")


if __name__ == "__main__":
    main()
