#!/usr/bin/env python3
"""apply_curation.py — write curated per-run values into runs.tsv, with the chain intact.

Every value this writes can be traced back to a sentence in a named paper. That is the whole
point: a cohort table saying 248 maize runs were collected at Baihe in 2016 is worthless if
nobody can say how it knows.

The chain, and the file that holds each link:

    data/runs.tsv         location = "China: Shanghai, Baihe Experimental Station"
                          location_source = "luo2020-site"       <- a rule id
    data/curation.tsv     rule luo2020-site: for PRJNA306542 where tissue=Leaf and
                          dev_stage=heading, set location to that value, evidence ev-luo2020-site
    data/provenance.tsv   ev-luo2020-site: Luo 2020, 10.1111/eva.13054, Methods 2.2,
                          quote "...drought resistance screening facility at Baihe
                          Experimental Station, Shanghai, in 2016 (May to October)"
    data/decisions.tsv    why this study is in scope at all, who decided, when

So `location_source` is no longer a bare token like "biosample". It names the rule that put
the value there, and the rule names the evidence, and the evidence quotes the paper.

**Rules are predicates on BioSample fields, not accession ranges.** That is what the data
actually needs. PRJNA306542 holds six experiments from six papers, and they are separated by
`tissue` and `dev_stage`, not by contiguous accession blocks: `Leaf`/`heading` is Luo's 248,
`leaf`/`Panicle initiation Period` is Ma's 51, `mesocotyl` is a lab light-exposure study.
`run_scope` in the link table cannot express any of that.

**Nothing is overwritten silently.** A rule that would replace an existing non-curated value
reports the conflict and skips, unless the rule sets `override: yes`. The case this guards
against is real: `geo_loc_name` for Ma's 51 runs reads "Brazil, China, Cote d'Ivoire", which
is where the CULTIVARS came from, and the paper says the plants were in Shanghai. Both are
"locations" and only one is a collection site.

Usage:

    python 01a_Literature/apply_curation.py --dry-run
    python 01a_Literature/apply_curation.py
    python 01a_Literature/apply_curation.py --explain SRR9129858
"""
import argparse, collections, csv, json, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _paths import ROOT

DATA       = HERE / "data"
RUNS       = DATA / "runs.tsv"
CURATION   = DATA / "curation.tsv"
PROVENANCE = DATA / "provenance.tsv"

# Fields a rule may set, each paired with its _source column in runs.tsv.
SETTABLE = {"tissue": "tissue_source", "location": "location_source",
            "collection_date": "date_source", "setting": "setting_source"}


def read(p: Path) -> list[dict]:
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def attrs_for(bp: str) -> dict:
    f = DATA / "biosample_attrs" / f"{bp}.json"
    if not f.exists():
        return {}
    o = json.loads(f.read_text())
    return o.get("attrs", o) if isinstance(o, dict) else {}


def matches(rule: dict, run: dict, a: dict) -> bool:
    """Does this run satisfy the rule's predicate?

    Clauses are joined by `;` and all must hold. Two operators:

        field=value     exact equality
        field~regex     regular expression search, case-sensitive

    Each is evaluated against the run's BioSample attributes first and its runinfo columns
    second. An empty predicate matches every run in the BioProject, which is what a
    whole-project rule looks like.

    The regex form exists for PRJNA383416, where 1,960 runs share ONE BioSample and the only
    per-run handle is a token inside `LibraryName`: KAKRNA6_10346261_GRoot_B79_CGTCGC. Those
    tokens separate a growth chamber from a greenhouse from the Musgrave field, so without
    substring matching the largest field block we have found cannot be expressed at all.

    Token boundaries matter and the delimiter is `_`, so write `(^|_)LMAD8(_|$)` rather than
    `LMAD8`: a bare `LMAD` would also swallow `LMAD26`, silently merging two collection dates
    two and a half weeks apart.
    """
    if run["BioProject"] != rule["BioProject"]:
        return False
    pred = (rule.get("predicate") or "").strip()
    if not pred:
        return True
    for clause in (c.strip() for c in pred.split(";") if c.strip()):
        if "~" in clause and ("=" not in clause or clause.index("~") < clause.index("=")):
            k, rx = (x.strip() for x in clause.split("~", 1))
            if not re.search(rx, a.get(k, run.get(k, "")) or ""):
                return False
        elif "=" in clause:
            k, v = (x.strip() for x in clause.split("=", 1))
            if (a.get(k, run.get(k, "")) or "") != v:
                return False
        else:
            return False
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--explain", metavar="RUN", help="trace every curated value on one run")
    args = ap.parse_args()

    runs = read(RUNS)
    rules = read(CURATION)
    prov = {p["evidence_id"]: p for p in read(PROVENANCE)}
    if not rules:
        sys.exit(f"{CURATION} is empty; nothing to apply")

    by_run = {r["Run"]: r for r in runs}
    attr_cache = {bp: attrs_for(bp) for bp in {r["BioProject"] for r in runs}}

    if args.explain:
        r = by_run.get(args.explain)
        if not r:
            sys.exit(f"{args.explain} not in runs.tsv")
        print(f"{r['Run']}   {r['BioProject']}   BioSample {r['BioSample']}")
        for field, srccol in SETTABLE.items():
            val, src = r.get(field, ""), r.get(srccol, "")
            print(f"\n  {field} = {val!r}\n    source: {src or '(none)'}")
            rule = next((x for x in rules if x["rule_id"] == src), None)
            if rule:
                print(f"    rule {rule['rule_id']}: where {rule['predicate'] or '(all runs)'}")
                ev = prov.get(rule.get("evidence_id", ""))
                if ev:
                    print(f"    evidence {ev['evidence_id']}  {ev['doi']}  {ev['locus']}")
                    print(f"      \"{ev['quote'][:220]}\"")
                    print(f"    entered by {ev['entered_by']} on {ev['ts'][:10]}")
        return

    applied = collections.Counter()
    conflicts = []
    for rule in rules:
        a_all = attr_cache.get(rule["BioProject"], {})
        field = rule["field"]
        if field not in SETTABLE:
            sys.exit(f"rule {rule['rule_id']}: field {field!r} is not settable")
        if rule.get("evidence_id") and rule["evidence_id"] not in prov:
            sys.exit(f"rule {rule['rule_id']}: evidence_id {rule['evidence_id']!r} "
                     f"is not in provenance.tsv — refusing to write an unevidenced value")
        srccol = SETTABLE[field]
        for r in runs:
            a = a_all.get(r["BioSample"], {})
            if not matches(rule, r, a):
                continue
            existing, esrc = r.get(field, ""), r.get(srccol, "")
            if existing and esrc not in ("", "not_fetched", "absent", rule["rule_id"]):
                conflicts.append((rule["rule_id"], r["Run"], field, existing, rule["value"]))
                if (rule.get("override") or "").lower() not in ("yes", "true", "1"):
                    continue
            r[field] = rule["value"]
            r[srccol] = rule["rule_id"]
            applied[rule["rule_id"]] += 1

    print(f"{len(rules)} rules over {len(runs):,} runs\n")
    for rid, n in applied.most_common():
        rule = next(x for x in rules if x["rule_id"] == rid)
        print(f"  {rid:<26}{n:>6} runs   {rule['field']} = {rule['value'][:46]}")
    if conflicts:
        print(f"\n  {len(conflicts)} conflicts with an existing curated value "
              f"(skipped unless override=yes):")
        for rid, run, f, old, new in conflicts[:6]:
            print(f"    {rid} {run} {f}: {old!r} -> {new!r}")
    if args.dry_run:
        print("\ndry run; runs.tsv not written")
        return

    with open(RUNS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(runs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(runs)
    print(f"\nwrote {RUNS.relative_to(ROOT)}")
    print("trace any value with:  python 01a_Literature/apply_curation.py --explain <RUN>")


if __name__ == "__main__":
    main()
