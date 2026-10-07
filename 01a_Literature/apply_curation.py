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
from _paths import ROOT          # repo-wide
from _layout import CURATION, DATA, JOINS, PROVENANCE, STUDIES

RUNS       = DATA / "runs.tsv"

# A join must match most of the project's samples or it is not describing them. PRJNA1217477's
# supplement matched 0 of 450 because it tabulated inoculum isolates rather than the sequenced
# plants; merging it would have stamped Californian vineyard coordinates onto Arabidopsis.
MIN_JOIN_RATE = 0.60

# How much a value is worth, by where it came from. A writer at a higher tier replaces a
# lower one FREELY; replacing an equal or higher tier needs an explicit `override: yes`.
#
# This is structural on purpose. Before it, every rule needed override=yes to beat even an
# archive value, so precedence lived in whether the author remembered the column rather than
# in the code — and a hand rule read out of a paper could be silently blocked by a bulk
# import that had got there first. Reading the paper is the whole point of the manual pass;
# it must win.
TIER = {"": 0, "absent": 0, "not_fetched": 0,
        "biosample": 1,          # the archive record
        "_import": 2,            # import_provenance.py bulk batches
        "_join": 3,              # a rate-checked join to a paper's supplement
        "_rule": 4}              # a hand-written rule quoting a paper's methods


def tier_of(src: str, rule_ids: set, join_ids: set) -> int:
    if src in rule_ids:
        return TIER["_rule"]
    if src in join_ids:
        return TIER["_join"]
    if src.startswith("import-"):
        return TIER["_import"]
    return TIER.get(src, 1)


def may_write(new_src: str, held_src: str, rule: dict,
              rule_ids: set, join_ids: set) -> bool:
    """May a writer replace what is already there?"""
    if held_src in ("", "absent", "not_fetched") or held_src == new_src:
        return True
    hi, lo = tier_of(new_src, rule_ids, join_ids), tier_of(held_src, rule_ids, join_ids)
    if hi > lo:
        return True
    return (rule.get("override") or "").lower() in ("yes", "true", "1")

# Fields a rule may set, each paired with its _source column in runs.tsv.
SETTABLE = {"tissue": "tissue_source", "location": "location_source",
            "collection_date": "date_source", "setting": "setting_source",
            "sampling_selection": "sampling_selection_source",
            "library_representative": "library_representative_source"}


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


def _sra_key(spec: str, run: dict, a: dict) -> str:
    """The join key on the SRA side.

    Either a bare column name (`SampleName`), or `column~regex` where capture group 1 is the
    key. The regex form is needed because submitters bury the key in free text: PRJDB7234's
    DDBJ descriptions end "... Sample ID: 20001", and that integer is what the paper's
    supplementary table is keyed on. Nothing in the structured fields carries it.
    """
    if "~" in spec:
        col, rx = (x.strip() for x in spec.split("~", 1))
        m = re.search(rx, a.get(col, run.get(col, "")) or "")
        return m.group(1) if m and m.groups() else ""
    return (a.get(spec.strip(), run.get(spec.strip(), "")) or "").strip()


def load_join_table(rule: dict) -> dict:
    """{key -> row} from a supplement, using the rule's file, sheet and header row."""
    path = STUDIES / rule["dir"] / rule["file"]
    if not path.exists():
        print(f"  {rule['rule_id']}: no such file {path}", file=sys.stderr)
        return {}
    try:
        import pandas as pd
    except ImportError:
        sys.exit("pandas needed for join rules")
    hdr = int(rule.get("header_row") or 0)
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path, sheet_name=rule["sheet"] or 0, header=hdr)
    else:
        df = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",", header=hdr)
    kc = rule["key_col"]
    if kc not in df.columns:
        print(f"  {rule['rule_id']}: key column {kc!r} not in {list(df.columns)[:8]}",
              file=sys.stderr)
        return {}
    # The supplement's key often needs reshaping to meet SRA's. Harris's table is keyed
    # `A1Y1_001_L` while SRA's LibraryName is `1L`: same sample, different spelling of the
    # same vine number and tissue letter. `key_regex` names a pattern whose capture groups are
    # concatenated to form the key, so `_0*(\d+)_([LR])$` turns one into the other and the
    # leading-zero strip happens in the pattern rather than in code.
    krx = re.compile(rule["key_regex"]) if rule.get("key_regex") else None
    out, unmatched = {}, 0
    for rec in df.dropna(subset=[kc]).to_dict("records"):
        k = rec[kc]
        k = str(int(k)) if isinstance(k, float) and k.is_integer() else str(k).strip()
        if krx:
            m = krx.search(k)
            if not m:
                unmatched += 1
                continue
            k = "".join(g for g in m.groups() if g is not None)
        out[k] = rec
    if unmatched:
        print(f"  {rule['rule_id']}: key_regex did not match {unmatched} of "
              f"{unmatched + len(out)} table rows", file=sys.stderr)
    return out


def apply_joins(runs: list[dict], attr_cache: dict, prov: dict, dry: bool,
                rule_ids: set, join_ids: set) -> tuple[collections.Counter, list]:
    """Apply every join rule: match SRA runs to supplement rows, then set fields.

    Reports the match rate and REFUSES below MIN_JOIN_RATE rather than merging a table that is
    not about these samples. A join that matches nothing is the single most dangerous thing in
    this module, because a merge still "succeeds" and writes confident wrong values.
    """
    applied, problems = collections.Counter(), []
    rules = read(JOINS)
    for rule in rules:
        if rule.get("evidence_id") and rule["evidence_id"] not in prov:
            sys.exit(f"join {rule['rule_id']}: evidence_id {rule['evidence_id']!r} not in "
                     f"provenance.tsv")
        field = rule["field"]
        if field not in SETTABLE:
            sys.exit(f"join {rule['rule_id']}: field {field!r} is not settable")
        table = load_join_table(rule)
        if not table:
            problems.append((rule["rule_id"], "table unreadable or key column missing"))
            continue
        a_all = attr_cache.get(rule["BioProject"], {})
        mine = [r for r in runs if r["BioProject"] == rule["BioProject"]]
        hits = 0
        for r in mine:
            k = _sra_key(rule["sra_key"], r, a_all.get(r["BioSample"], {}))
            if k and k in table:
                hits += 1
        rate = hits / len(mine) if mine else 0
        if rate < MIN_JOIN_RATE:
            problems.append((rule["rule_id"],
                             f"REFUSED: joined {hits}/{len(mine)} ({rate:.0%}), "
                             f"below {MIN_JOIN_RATE:.0%}. Is this table about these samples?"))
            continue
        srccol = SETTABLE[field]
        for r in mine:
            k = _sra_key(rule["sra_key"], r, a_all.get(r["BioSample"], {}))
            rec = table.get(k)
            if not rec:
                continue
            val = rec.get(rule["value_col"])
            if val is None or str(val).strip() in ("", "nan"):
                continue
            val = str(int(val)) if isinstance(val, float) and val.is_integer() else str(val).strip()
            if rule.get("value_prefix"):
                val = rule["value_prefix"] + val
            held = r.get(srccol, "")
            if r.get(field) and not may_write(rule["rule_id"], held, rule,
                                              rule_ids, join_ids):
                continue
            if not dry:
                r[field], r[srccol] = val, rule["rule_id"]
            applied[rule["rule_id"]] += 1
        print(f"  {rule['rule_id']:<26}joined {hits}/{len(mine)} ({rate:.0%})  "
              f"{field} <- {rule['value_col']}")
    return applied, problems


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


def check_sources(prov: list[dict]) -> list[str]:
    """Warn where an evidence row's paper has no copy on disk.

    An evidence row is a quote attributed to a paper. If nobody can open that paper locally,
    the quote cannot be checked against it, which makes the row the weakest kind of link in
    the chain: plausible, attributed, and unverifiable. This does not block a write - the
    quote may well be right, and some sources are paywalled - but it is worth saying out loud
    every run rather than discovering it months later.

    Added 2026-10-06, after an ad-hoc version of this check found that seven of the fourteen
    papers curated that day had been read from Europe PMC into a temporary directory and never
    filed. Directories are matched on the DOI slug AND on the paper_key, because some
    directories predate the DOI naming convention (see studies/Cai21b/NAMING.txt), and the
    paper_ref is resolved through papers.tsv because some provenance rows carry the DOI in
    paper_key rather than the short ref.
    """
    refs = {}
    for r in read(DATA / "papers.tsv"):
        d, ref = (r.get("doi") or "").strip(), (r.get("paper_ref") or "").strip()
        if d and ref:
            refs[d] = ref
    seen, warn = {}, []
    for r in prov:
        doi, key = (r.get("doi") or "").strip(), (r.get("paper_key") or "").strip()
        if not doi or doi == "(multiple)":
            continue
        if doi in seen:
            continue
        cands = [STUDIES / ("doi_" + doi.replace("/", "_"))]
        for alt in (key, refs.get(doi)):
            if alt and alt != doi:
                cands.append(STUDIES / alt)
        ok = any(d.is_dir() and any(f.is_file() and f.name != "adapter.yaml"
                                    for f in d.rglob("*")) for d in cands)
        seen[doi] = ok
        if not ok:
            warn.append(f"{doi}  ({key or 'no paper_key'}), {sum(1 for x in prov if x.get('doi') == doi)} evidence rows")
    return warn


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--explain", metavar="RUN", help="trace every curated value on one run")
    args = ap.parse_args()

    runs = read(RUNS)
    rules = read(CURATION)
    prov_rows = read(PROVENANCE)
    prov = {p["evidence_id"]: p for p in prov_rows}
    if not rules:
        sys.exit(f"{CURATION} is empty; nothing to apply")

    unsourced = check_sources(prov_rows)
    if unsourced:
        print(f"WARNING: {len(unsourced)} evidenced paper(s) have no copy under studies/, so "
              f"their quotes cannot be checked:", file=sys.stderr)
        for w in unsourced:
            print(f"  {w}", file=sys.stderr)
        print(file=sys.stderr)

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
            rule = next((x for x in rules + read(JOINS) if x["rule_id"] == src), None)
            # import_provenance.py writes an evidence_id straight into the _source column:
            # a bulk import is not a predicate over runs, so it has no rule to name. Resolve
            # it directly, or 2,407 imported values read as unevidenced.
            if rule is None and src in prov:
                ev = prov[src]
                print(f"    bulk import {src}")
                print(f"    evidence {ev['evidence_id']}  {ev['locus']}")
                print(f"      \"{ev['quote'][:260]}\"")
                print(f"    entered by {ev['entered_by']} on {ev['ts'][:10]}")
                continue
            if rule:
                if rule.get("sra_key"):
                    print(f"    join {rule['rule_id']}: {rule['file']} sheet {rule['sheet']}, "
                          f"{rule['key_col']} <-> {rule['sra_key']}, value {rule['value_col']}")
                else:
                    print(f"    rule {rule['rule_id']}: where {rule['predicate'] or '(all runs)'}")
                ev = prov.get(rule.get("evidence_id", ""))
                if ev:
                    print(f"    evidence {ev['evidence_id']}  {ev['doi']}  {ev['locus']}")
                    print(f"      \"{ev['quote'][:220]}\"")
                    print(f"    entered by {ev['entered_by']} on {ev['ts'][:10]}")
        return

    applied = collections.Counter()
    conflicts = []
    replaced = collections.Counter()
    rule_ids = {x["rule_id"] for x in rules}
    join_ids = {x["rule_id"] for x in read(JOINS)}
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
            if existing and not may_write(rule["rule_id"], esrc, rule, rule_ids, join_ids):
                conflicts.append((rule["rule_id"], r["Run"], field, existing, rule["value"]))
                continue
            if existing and esrc != rule["rule_id"]:
                replaced[(rule["rule_id"], esrc)] += 1
            r[field] = rule["value"]
            r[srccol] = rule["rule_id"]
            applied[rule["rule_id"]] += 1

    japplied, jproblems = apply_joins(runs, attr_cache, prov, args.dry_run,
                                      rule_ids, join_ids)
    for rid, n in japplied.items():
        applied[rid] += n
    for rid, msg in jproblems:
        print(f"  !! {rid}: {msg}", file=sys.stderr)

    print(f"{len(rules)} predicate rules + {len(read(JOINS))} join rules "
          f"over {len(runs):,} runs\n")
    all_rules = rules + read(JOINS)
    for rid, n in applied.most_common():
        rule = next(x for x in all_rules if x["rule_id"] == rid)
        # join rules have no literal `value`; they report the supplement column instead
        shown = rule.get("value") or f"<- {rule.get('value_col','?')} (join)"
        print(f"  {rid:<26}{n:>6} runs   {rule['field']} = {shown[:46]}")
    if replaced:
        print("\n  replaced a lower-tier value (precedence working as intended):")
        for (new, old), n in replaced.most_common(8):
            print(f"    {n:>6} x  {new} over {old}")
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
