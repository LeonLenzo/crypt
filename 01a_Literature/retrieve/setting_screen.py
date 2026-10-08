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

from _layout import BIOPROJECT_XML, BIOSAMPLE, RUNINFO

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


def classify(bp: str) -> dict:
    ctrl, fld, mixed = [], [], []
    for locus, text in project_text(bp):
        for s in sentences(text):
            # A BioSample attribute is a value, not prose: the field IS the statement,
            # so the link requirement cannot apply to it.
            linked = locus.startswith("biosample/") or bool(SAMPLING.search(s))
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
                n_sources=len(project_text(bp)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", nargs="*")
    ap.add_argument("--file", help="one accession per line")
    ap.add_argument("--quotes", action="store_true", help="print the matched sentences")
    ap.add_argument("--only", help="print only this call (controlled/field/mixed/none)")
    a = ap.parse_args()

    accs = list(a.accessions)
    if a.file:
        accs += [l.strip() for l in Path(a.file).read_text().split() if l.strip()]
    if not accs:
        sys.exit("usage: setting_screen.py <BioProject>... | --file list.txt")

    results = [classify(bp) for bp in accs]
    tally = collections.Counter(r["call"] for r in results)
    print("  ".join(f"{k}:{tally[k]}" for k in ("controlled", "field", "mixed", "none")))
    print()
    for r in results:
        if a.only and r["call"] != a.only:
            continue
        print(f"{r['bp']:<14}{r['call']:<11}{r['n_sources']:>4} sources")
        if a.quotes:
            for tag, hits in (("CTRL", r["controlled"]), ("FIELD", r["field"]),
                              ("MIXED", r["mixed"])):
                for locus, s in hits[:3]:
                    print(f"    {tag:<6}{locus:<26}{s[:170]}")


if __name__ == "__main__":
    main()
