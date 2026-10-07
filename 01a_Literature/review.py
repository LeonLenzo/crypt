#!/usr/bin/env python3
"""review.py — step through candidate studies and record a human decision on each.

The decision is yours. Everything around it is not: ranking the queue, assembling the
evidence, proposing an answer with its reasoning, writing the record, and never asking the
same question twice. This is the stepper that does the rest.

**The record, not the interface, is the deliverable.** Decisions append to
`data/decisions.tsv` with a reason, the evidence behind them, and a date, on the same terms
as `_exclusions.py`: an entry without evidence is an opinion. The file is append-only, so a
changed mind shows up as a change rather than silently replacing what you thought before.
That file is what the methods section is written from.

**Three levels, and only one is per-sample.**

    paper       is this field bulk RNA-seq at all?            one call
    bioproject  which of this paper's projects are field?     one call per project
    sample      which runs within a project?                  a RULE, never a row

The third never comes here. PRJNA383416 splits 740 field against 1,220 controlled, and you
resolve that by writing `LMAD, LMAN, Kern -> field` once in an adapter, not by judging 1,960
runs. The adapter IS the sample-level decision record.

**Nothing is ever dropped from the data.** A `drop` decision flags runs, it does not remove
them, so the denominator stays computable and any rate can be reported with and without an
exclusion. That is the only way to show an exclusion is not producing the result.

Queue order is deliberate. `--need none` first: 66 papers that need no supplement, pure
yes/no, each immediately adding runs to the cohort. They calibrate your criteria cheaply
before you spend an afternoon on a supplement. Deciding the hard ones first is how criteria
drift.

Usage:

    python 01a_Literature/review.py --show            # the next card, no prompt
    python 01a_Literature/review.py --need none       # work the quick queue
    python 01a_Literature/review.py --key Ada21       # jump to one study
    python 01a_Literature/review.py --status          # what is decided, what is left
"""
import argparse, collections, csv, datetime, json, re, sys, textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _paths import ROOT          # repo-wide
from _layout import SEEDS      # this module's layer dirs
from _exclusions import EXCLUDED

DATA      = HERE / "data"
STUDIES   = HERE / "studies"
DECISIONS = SEEDS / "decisions.tsv"
FOUND     = SEEDS / "found.tsv"

DECISION_FIELDS = ["ts", "level", "key", "decision", "reason", "evidence", "decided_by"]
FOUND_FIELDS    = ["ts", "paper_key", "kind", "value", "note", "found_by"]

# Facts are not judgments, and keeping them apart is the point. A DOI, an accession or a
# supplement path is something the pipeline can act on: triage.py re-reads found.tsv and the
# study changes state. "keep" and "drop" are terminal human calls that no amount of fetching
# will produce. Collapsing the two is how `accession-hunt` came to mean "I went looking and
# found a PDF", which is a fact wearing a verdict's clothes.
FOUND_KINDS = {
    "doi":        "a DOI for a paper we had none for",
    "bioproject": "a BioProject accession named in the paper or its data statement",
    "other_acc":  "a GEO/SRA/ArrayExpress accession, to be resolved to a BioProject",
    "supplement": "a supplementary file now sitting in this study's directory",
    "source_doc": "another document that holds the metadata (a thesis, a data paper)",
    "none":       "looked for it and it is genuinely not there",
}

# What a decision can be. `split` means the study is partly in scope and an adapter has to
# carry the rule; `accession-hunt` means we cannot judge it until an accession is found.
CHOICES = {
    "k": ("keep",           "in scope; its runs join the cohort"),
    "d": ("drop",           "out of scope; runs are FLAGGED, never deleted"),
    "s": ("split",          "partly in scope; needs an adapter rule to say which runs"),
    "a": ("accession-hunt", "cannot judge until a data accession is found"),
    "f": ("defer",          "come back to it"),
    "b": ("blocked",        "cannot resolve: no full text, no data, dead accession"),
}


def read(name: str) -> list[dict]:
    p = DATA / name
    if not p.exists():
        sys.exit(f"{p} missing; run triage.py and scaffold.py first")
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def decisions() -> dict:
    """Latest decision per key, with retractions removed.

    `retracted` is not a verdict, it is the withdrawal of one, so a retracted study returns to
    the queue as though it had never been decided. The row stays in the file: the history of
    what was claimed and then withdrawn is worth more than a tidy table.
    """
    if not DECISIONS.exists():
        return {}
    out = {}
    with open(DECISIONS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r["level"], r["key"])
            if r["decision"] == "retracted":
                out.pop(k, None)
            else:
                out[k] = r
    return out


def record(level: str, key: str, decision: str, reason: str, evidence: str,
           who: str = "leon") -> None:
    """Append a decision. Verdicts are Leon's; Claude may only append `retracted`.

    Enforced rather than documented, because on 2026-10-05 Claude recorded a `drop` on
    10.1038/srep37302 after assembling the evidence for it. Assembling evidence and reaching a
    verdict are different jobs, and the second one is the whole reason this file exists.
    """
    if who.startswith("claude") and decision != "retracted":
        raise PermissionError(
            f"refusing to record '{decision}' as {who}: verdicts are the human's. "
            "Present the evidence and let Leon decide.")
    new = not DECISIONS.exists()
    DECISIONS.parent.mkdir(parents=True, exist_ok=True)
    with open(DECISIONS, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=DECISION_FIELDS, delimiter="\t")
        if new:
            w.writeheader()
        w.writerow(dict(ts=datetime.datetime.now().isoformat(timespec="seconds"),
                        level=level, key=key, decision=decision,
                        reason=reason, evidence=evidence, decided_by=who))


def registered_exclusions(projects: list[dict]) -> list[tuple[str, dict]]:
    """Projects already in _exclusions.py, with the recorded reason.

    The registry exists so a decision made once stays made. A stepper that does not read it
    would hand back PRJNA1217477, judged unusable on 2026-10-02 with its evidence written
    down, as though it were an open question.
    """
    return [(p["BioProject"], EXCLUDED[p["BioProject"]])
            for p in projects if p["BioProject"] in EXCLUDED]


def note_found(key: str, kind: str, value: str, note: str, who: str = "leon") -> None:
    new = not FOUND.exists()
    FOUND.parent.mkdir(parents=True, exist_ok=True)
    with open(FOUND, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FOUND_FIELDS, delimiter="\t")
        if new:
            w.writeheader()
        w.writerow(dict(ts=datetime.datetime.now().isoformat(timespec="seconds"),
                        paper_key=key, kind=kind, value=value, note=note, found_by=who))


def missing_fields(p: dict) -> list[str]:
    """Which per-run fields this project is short of, by name.

    `geo 100% date 0%` makes the reader do the inference. Naming the field is what tells you
    what to go and look for: PRJNA746402 does not need an accession, it needs collection dates.
    """
    out = []
    if p["triage"] in ("one-biosample", "few-biosamples"):
        return ["tissue", "location", "collection_date"]
    for field, col in (("location", "geo_pct"), ("collection_date", "date_pct")):
        v = p.get(col, "")
        if v == "" or (v.replace(".", "").isdigit() and float(v) < 90):
            out.append(field)
    return out


def suggest(w: dict, projects: list[dict]) -> tuple[str, str]:
    """A proposed decision and the evidence for it. Confirming is fast; deriving 229 times
    is not. Never applied automatically: this fills the prompt, it does not answer it."""
    dec = (w.get("undermind_decision") or "").strip()
    ex = registered_exclusions(projects)
    if ex and len(ex) == len(projects):
        return "drop", f"all projects registered in _exclusions.py: {ex[0][1]['kind']}"
    if ex:
        return "split", (f"{len(ex)} of {len(projects)} projects registered in "
                         f"_exclusions.py; the rest need a separate call")
    if w["need"] == "no-accession":
        return "accession-hunt", "no BioProject linked; check the data-availability statement"
    if dec == "does_not_meet":
        return "drop", "Undermind screened it as not meeting the field/aerial criteria"
    if dec == "duplicate_dataset":
        return "drop", "Undermind flagged it as a reanalysis of an already-counted dataset"
    states = {p["triage"] for p in projects}
    if dec == "meets" and states <= {"sra-complete"}:
        return "keep", "Undermind says meets, and every project has per-run geo and date"
    if states & {"one-biosample", "few-biosamples"}:
        return "split", "pooled or single BioSample; per-run scope needs an adapter rule"
    if states & {"sra-partial"}:
        gaps = sorted({f for p in projects for f in missing_fields(p)})
        return "keep", ("in scope but per-run metadata is incomplete; missing "
                        + ", ".join(gaps) + " — a supplement would fill it")
    if states <= {"sra-complete"}:
        return "keep", "every project has per-run geography and a collection date"
    return "", ""


_ACC_ANY = re.compile(r"\b(PRJ[DEN][A-Z]\d+|[SED]RP\d{5,}|GSE\d{3,}|E-[A-Z]{4}-\d+)\b", re.I)


def _ask(prompt: str) -> str:
    try:
        return input(f"  {prompt} > ").strip()
    except (EOFError, KeyboardInterrupt):
        raise KeyboardInterrupt


def interview(w: dict, projects: list[dict]) -> tuple[str, str, str] | None:
    """Walk the open questions for one study; return (decision, reason, evidence).

    Only what is actually unresolved gets asked. A study with a DOI, an accession and complete
    per-run metadata is one keystroke; a study missing all three is four. The flat six-way menu
    asked the same question of every study and offered no way to record an answer that was a
    FACT rather than a verdict, which is why "I found the paper" had to be filed as
    `accession-hunt`.

    Returns None if the study is skipped.
    """
    key = w["paper_key"]

    # 1. identity. Only when we have no DOI: without one the study cannot be looked up again.
    if not w["doi"]:
        a = _ask("no DOI on record. paste one, or [n] none findable / [s] skip")
        if a.lower() == "s":
            return None
        if a.lower() == "n":
            note_found(key, "none", "doi", "no DOI findable for this paper")
        elif a:
            note_found(key, "doi", a, "added during review")
            print("     recorded. re-run triage.py to pick it up.")

    # 2. data. Only when nothing is linked yet.
    if not projects:
        a = _ask("no accession linked. paste a BioProject/SRP/GSE/E-MTAB, "
                 "or [n] none in the paper / [s] skip")
        if a.lower() == "s":
            return None
        if a.lower() == "n":
            note_found(key, "none", "bioproject", "no accession stated in the paper")
            return ("blocked", "no data accession stated in the paper", "data-availability statement")
        if a:
            m = _ACC_ANY.search(a)
            if not m:
                print("     that does not look like an accession; not recorded.")
            else:
                acc = m.group(1).upper()
                kind = "bioproject" if acc.startswith("PRJ") else "other_acc"
                note_found(key, kind, acc, "found during review")
                print(f"     recorded {acc}. re-run triage.py "
                      + ("(and resolve_accessions.py) " if kind == "other_acc" else "")
                      + "to pick it up.")

    # 3. scope. The judgment, always asked.
    a = _ask("in scope (field, aerial, bulk RNA-seq)? "
             "[y] keep  [n] drop  [p] partly  [u] unsure, defer  [s] skip").lower()
    if a == "s":
        return None
    if a == "u":
        return ("defer", _ask("what would settle it?"), "")
    if a == "n":
        return ("drop", _ask("why out of scope?"), _ask("evidence (section, sentence)"))
    decision = "split" if a == "p" else "keep"
    reason = _ask("reason")

    # 4. metadata. Only when something is actually missing.
    gaps = sorted({f for p in projects for f in missing_fields(p)})
    if gaps:
        a = _ask(f"needs {', '.join(gaps)}. supplement filename in this directory, "
                 f"or [n] none available / [enter] leave open")
        if a.lower() == "n":
            note_found(key, "none", "supplement", f"no source found for {', '.join(gaps)}")
        elif a:
            d = STUDIES / w["dir"]
            kind = "supplement" if (d / a).exists() else "source_doc"
            if not (d / a).exists():
                print(f"     no file named {a} in {d}; recorded as a pointer, not a file.")
            note_found(key, kind, a, f"supplies {', '.join(gaps)}")

    return (decision, reason, _ask("evidence (section, table, sentence)"))


def card(w: dict, projects: list[dict], runs_by_bp: dict, prior: dict | None) -> str:
    L = []
    add = L.append
    add("=" * 78)
    add(f"  {w.get('title_short') or w['paper_key']}")
    add("=" * 78)
    if w.get("paper_ref"):
        add(f"  ref {w['paper_ref']}")
    add(f"  doi  {('https://doi.org/' + w['doi']) if w['doi'] else '(none recorded)'}")
    if prior:
        add(f"  ALREADY DECIDED {prior['decision']} on {prior['ts'][:10]}: {prior['reason']}")
        add("")
    for bp, e in registered_exclusions(projects):
        add(f"  !! {bp} IS A REGISTERED EXCLUSION ({e['kind']}, decided {e['decided']})")
        for line in textwrap.wrap(e["reason"], 70):
            add(f"     {line}")
        add("")
    # The single most important two lines on the card. The accession WAS shown before, buried
    # in a block of run counts three lines down, while the menu offered `accession-hunt`. That
    # combination sent someone hunting for an accession we already had. State plainly what is
    # in hand and what is actually missing, before anything else.
    have, need = [], []
    (have if w["doi"] else need).append("paper DOI")
    if projects:
        have.append("accession " + ", ".join(
            f"{p['BioProject']} ({p['sra_runs']} runs)" for p in projects))
    else:
        need.append("a data accession")
    gaps = sorted({f for p in projects for f in missing_fields(p)})
    for f in ("tissue", "location", "collection_date"):
        if projects and f not in gaps:
            have.append(f)
    need += gaps
    add(f"  HAVE   {'; '.join(have) if have else 'nothing'}")
    add(f"  NEED   {', '.join(need) if need else 'nothing — this study is complete'}")
    add("")
    add(f"  need            {w['need']}")
    add(f"  Undermind       {w.get('undermind_decision') or '(not screened)'}")
    add(f"  runs affected   {int(w['runs_affected'] or 0):,}"
        f"     never examined {int(w['runs_unexamined'] or 0):,}")
    add("")
    for p in projects:
        add(f"  {p['BioProject']:<14} {p['triage']:<14} {p['organism'][:30]}")
        add(f"  {'':<14} {p['sra_runs']:>5} runs, {p['sra_biosamples']:>4} BioSamples, "
            f"{p['gate_passed']:>4} passed the old gate")
        gaps = missing_fields(p)
        add(f"  {'':<14} geo {p['geo_pct'] or '?'}%  date {p['date_pct'] or '?'}%"
            f"   -> needs: {', '.join(gaps) if gaps else 'nothing, complete'}")
        add(f"  {'':<14} primary: {p['primary_paper'] or '(unresolved)'} "
            f"[{p['primary_source']}]")
        if p.get("primary_conflict"):
            add(f"  {'':<14} ** PRIMARY CONFLICT: generating paper and resolved DOI disagree")
        rr = runs_by_bp.get(p["BioProject"], [])
        if p["triage"] in ("one-biosample", "few-biosamples") and rr:
            keys = [c for c in ("LibraryName", "SampleName", "Experiment")
                    if len({r[c] for r in rr if r.get(c)}) == len(rr)]
            add(f"  {'':<14} join on: {', '.join(keys) or 'Run/Experiment only'}")
        add("")
    raw = (w.get("undermind_decision_raw") or "").strip()
    if raw:
        add("  Undermind's reasoning:")
        for line in textwrap.wrap(raw, 70):
            add(f"    {line}")
        add("")
    # Create the drop zone on sight rather than only in scaffold.py. A paper arriving in a
    # later Undermind batch would otherwise have nowhere to put its PDF until someone
    # remembered to re-scaffold, and "there was no directory" is a silly reason to lose a
    # supplement someone just downloaded.
    d = STUDIES / w["dir"]
    d.mkdir(parents=True, exist_ok=True)
    supp = [f.name for f in d.iterdir()
            if f.is_file() and ":Zone.Identifier" not in f.name
            and f.suffix.lower() in (".xls", ".xlsx", ".csv", ".tsv")]
    pdf = [f.name for f in d.iterdir()
           if f.is_file() and f.suffix.lower() == ".pdf"]
    add(f"  dir             studies/{w['dir']}")
    add(f"  pdf             {', '.join(pdf) or '(none dropped in yet)'}")
    add(f"  supplement      {', '.join(supp) or '(none dropped in yet)'}")
    return "\n".join(L)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--need", help="restrict the queue to one need class")
    ap.add_argument("--key", help="jump to one paper_key, ref or directory name")
    ap.add_argument("--show", action="store_true", help="print the next card and exit")
    ap.add_argument("--status", action="store_true", help="progress summary")
    ap.add_argument("--all", action="store_true", help="include already-decided studies")
    args = ap.parse_args()

    work = read("worklist.tsv")
    # The title lives in papers.tsv and the card never showed it, so identifying a study meant
    # looking the DOI up by hand. That is how a search for Caicedo et al. 2021 came back with
    # a different author's doctoral thesis.
    titles = {r["paper_key"]: r.get("title", "") for r in read("papers.tsv")}
    for w in work:
        t = titles.get(w["paper_key"], "")
        m = re.search(r"[\u201c\"]([^\u201d\"]{10,150})[\u201d\"]", t)
        w["title_short"] = (m.group(1) if m else t)[:110]
        w["title_full"] = t
    bps = {b["BioProject"]: b for b in read("bioprojects.tsv")}
    runs_by_bp = collections.defaultdict(list)
    for r in read("runs.tsv"):
        runs_by_bp[r["BioProject"]].append(r)
    done = decisions()

    if args.status:
        n = collections.Counter(v["decision"] for (lvl, _), v in done.items() if lvl == "paper")
        print(f"decided {sum(n.values())} of {len(work)} papers")
        for k, v in n.most_common():
            print(f"  {k:<16}{v:>4}")
        left = collections.Counter(w["need"] for w in work
                                   if ("paper", w["paper_key"]) not in done)
        print("\nremaining by need:")
        for k, v in left.most_common():
            print(f"  {k:<16}{v:>4}")
        return

    queue = [w for w in work
             if (args.all or ("paper", w["paper_key"]) not in done)
             and (not args.need or w["need"] == args.need)]
    if args.key:
        queue = [w for w in work if args.key in (w["paper_key"], w.get("paper_ref"), w["dir"])]
        if not queue:
            sys.exit(f"no study matching {args.key!r}")
    if not queue:
        print("queue empty: everything matching is decided. --status for the summary.")
        return

    for w in queue:
        projects = [bps[a] for a in w["bioprojects"].split(";") if a in bps]
        prior = done.get(("paper", w["paper_key"]))
        print(card(w, projects, runs_by_bp, prior))
        prop, why = suggest(w, projects)
        if prop:
            print(f"  suggested: {prop}  ({why})")
        if args.show:
            return
        print()
        try:
            got = interview(w, projects)
        except KeyboardInterrupt:
            print("\n  stopped; nothing further recorded.")
            return
        if got is None:
            print("  skipped.\n")
            continue
        decision, reason, ev = got
        record("paper", w["paper_key"], decision, reason, ev)
        print(f"  recorded: {decision}\n")


if __name__ == "__main__":
    main()
