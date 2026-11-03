#!/usr/bin/env python3
"""The shared tail of a per-study resolver, and the ONE definition of a join key.

A resolver exists when a study hides its per-sample metadata somewhere the archive does not
carry it: a supplement, a GEO record, a design code. Each one's body is necessarily bespoke -
the refusal guards are its whole value, and they must stay specific to the study (Kashima
checks every supplement day against the archive's UTC day, Sato checks every library appears
exactly twice). What is NOT bespoke is the tail: write the CSV, then prove the table will
actually join to our runs.

Two things pushed this into a shared module on 2026-10-06.

**Three of the five resolvers never checked their join rate.** Only range28 and epicon did,
and each had its own copy of the arithmetic. A resolver that writes a beautiful table keyed on
something absent from `runs.tsv` reports success, and the failure only surfaces later as a
refusal inside `apply_curation.py` - or worse, as a 61% join nobody looked at. `finish()`
makes the check unconditional.

**The key resolution had drifted out of sync with the engine.** `apply_curation.py` owned
`_sra_key`, and the resolvers approximated it. The two must agree exactly, or a resolver's
reported rate is not the rate the engine will achieve. Both now import from here, so there is
one definition. The same argument applies to `attrs_for`: when the BioSample cache moved from
`data/` to `bronze/`, `apply_curation.py` kept a stale path and two joins covering 1,657 runs
silently fell to 0%. One loader means a path change breaks everything at once, loudly, rather
than one caller quietly.
"""

from __future__ import annotations

import csv, json, re, sys
from pathlib import Path

from _layout import BIOSAMPLE, ROOT, RUNS

# A join that matches nothing is the single most dangerous thing in this module, because the
# merge still "succeeds" and writes confident wrong values. Shared with apply_curation.py so
# the resolver's prediction and the engine's execution use the same bar.
MIN_JOIN_RATE = 0.60


def attrs_for(bp: str) -> dict:
    """The cached BioSample attributes for one BioProject, or {} if never fetched."""
    f = BIOSAMPLE / f"{bp}.json"
    if not f.exists():
        return {}
    o = json.loads(f.read_text())
    return o.get("attrs", o) if isinstance(o, dict) else {}


def sra_key(spec: str, run: dict, a: dict) -> str:
    """The join key on the SRA side.

    Either a bare column name (`SampleName`), or `column~regex` where capture group 1 is the
    key. The regex form is needed because submitters bury the key in free text: PRJDB7234's
    DDBJ descriptions end "... Sample ID: 20001", and that integer is what the paper's
    supplementary table is keyed on. Nothing in the structured fields carries it.

    BioSample attributes are consulted BEFORE the runinfo columns, which is why a design code
    living in `source_material_id` is a zero-regex join key.
    """
    if "~" in spec:
        col, rx = (x.strip() for x in spec.split("~", 1))
        m = re.search(rx, a.get(col, run.get(col, "")) or "")
        return m.group(1) if m and m.groups() else ""
    return (a.get(spec.strip(), run.get(spec.strip(), "")) or "").strip()


def finish(rows: list[dict], out: Path, *, key_col: str, sra_key_spec: str | None = None,
           bioprojects: list[str] | None = None, min_rate: float = MIN_JOIN_RATE,
           sort_key=None) -> int:
    """Write a resolver's join table, then prove it joins. Returns the matched run count.

    `key_col`   the column in THIS table that `joins.tsv` will key on.
    `sra_key_spec`  how the same key is found on the SRA side; defaults to `key_col`.
    `bioprojects`   which projects to measure the rate against. None measures every run,
                    which is right for a table spanning projects.

    Refuses rather than returns on: no rows, a non-unique key, or a rate below `min_rate`.
    """
    if not rows:
        sys.exit("REFUSED: no rows to write")
    spec = sra_key_spec or key_col
    if key_col not in rows[0]:
        sys.exit(f"REFUSED: key column {key_col!r} is not in the table; "
                 f"columns are {list(rows[0])}")

    keys = [str(r[key_col]).strip() for r in rows]
    if "" in keys:
        sys.exit(f"REFUSED: {keys.count('')} row(s) have an empty {key_col}")
    if len(set(keys)) != len(keys):
        dup = sorted({k for k in keys if keys.count(k) > 1})[:5]
        sys.exit(f"REFUSED: {key_col} is the join key and is not unique, e.g. {dup}")

    if sort_key:
        rows.sort(key=sort_key)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    with open(RUNS) as fh:
        every = list(csv.DictReader(fh, delimiter="\t"))
    cache = {bp: attrs_for(bp) for bp in {r["BioProject"] for r in every}}
    kset = set(keys)

    def matches(r: dict) -> bool:
        return sra_key(spec, r, cache.get(r["BioProject"], {}).get(r["BioSample"], {})) in kset

    if bioprojects is None:
        # A table may span several BioProjects (ada21 covers four). Scope the rate to the
        # projects this table actually reaches, so "did it cover its own projects" is the
        # question asked, rather than diluting it against all of runs.tsv.
        bioprojects = sorted({r["BioProject"] for r in every if matches(r)})
        if not bioprojects:
            sys.exit(f"REFUSED: not one run in runs.tsv matches this table on {spec!r}")
    ours = [r for r in every if r["BioProject"] in set(bioprojects)]
    if not ours:
        sys.exit(f"REFUSED: runs.tsv has no rows for {bioprojects}; the table cannot join")
    hit = sum(1 for r in ours if matches(r))
    rate = hit / len(ours)

    try:
        where = out.relative_to(ROOT)
    except ValueError:
        where = out
    print(f"\nwrote {where}: {len(rows)} rows", file=sys.stderr)
    print(f"  joins {hit}/{len(ours)} ({rate:.0%}) of our runs on {spec} "
          f"across {len(bioprojects)} project(s)", file=sys.stderr)
    if rate < min_rate:
        sys.exit(f"REFUSED: {rate:.0%} is below {min_rate:.0%}. The table is keyed on "
                 f"something {spec!r} does not carry, or it is not about these samples.")
    return hit
