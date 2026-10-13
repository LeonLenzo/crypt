#!/usr/bin/env python3
"""build_seeds.py — regenerate the GENERATED half of 03_kraken/seeds/ from holobase.

seeds/ is the single source of truth for what is in and what is out of the CDS database.
Half of it is AUTHORED and must never be overwritten by a script (sources.tsv, hosts.tsv,
exclusions.tsv, residents.tsv); half is a reproducible extraction from holobase. This script
writes only the second half. See README.md for which is which and why.

THE SCOPE RULE, and it is the thing most likely to be broken by a one-line change:

    Plant-associated, ANY host. NEVER scope by the cohort's own hosts.

USDA Fungus-Host Distributions is host-ASSOCIATION data. It is used here to ask "is this
fungus ever recorded on any plant", which is a scope gate. It must never be used to ask "is
it recorded on wheat/maize/rice/sorghum", which is a plausibility filter, because the study
exists to find UNREPORTED host associations and such a filter makes them undetectable by
construction. CLAUDE.md already rejects host-plausibility filtering of CALLS for this reason
(Leon's objection). Moving the same filter earlier, into reference selection, is worse: a
filter on calls can be relaxed afterwards, but a missing reference is invisible and
unrecoverable without a rebuild, because the reads go silently to the nearest taxon that IS
present. That failure mode is documented in 01_search/README.md.

Two traps in holobase, both of which returned confidently wrong answers before being caught:

  * `taxon.taxon_id` is an INTERNAL SURROGATE, not an NCBI taxid — the schema says so, and
    the real key is `taxon.ncbi_taxid`. Joining db_v3's taxids against taxon_id returns a
    confident 0% overlap instead of an error. Always join on ncbi_taxid.
  * `assembly.is_annotated` is 1 on all 36,709 rows, so holobase CANNOT tell you whether an
    assembly has a gene model. The database is CDS-only, so annotation is the binding
    constraint, and the only way to size a build is a real --select-only run. Availability
    counts overestimated the resident set by 3x (508 assemblies against the selector's 161).

Also: `taxon.kingdom` spells plants `Viridiplantae`, not `Plants`. Filtering on 'Plants'
silently returns zero.

Two passes, deliberately separate files, because they are two different kinds of truth:

    fungi_usda.tsv       a projection of holobase: WHICH fungi are plant-associated
    fungi_annotated.tsv  NCBI truth: which of them have a GENE MODEL, and so a CDS

holobase cannot answer the second question at all (`is_annotated` is 1 on all 36,709 rows),
and annotation is the only thing that matters for a CDS-only database: we do not care about
annotation for its own sake, we care because it is the only way to get CDS. So the second pass
asks NCBI directly, per genus rather than per taxon (584 calls instead of 1,488), and keeps
its answer in its own file so re-running one pass does not invalidate the other.

Run from the repository root:

    python 03_kraken/seeds/build_seeds.py            # -> fungi_usda.tsv  (holobase)
    python 03_kraken/seeds/build_seeds.py --check    # verify, write nothing
    python 03_kraken/seeds/build_seeds.py --annotate # -> fungi_annotated.tsv  (NCBI, ~10 min)
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def _root() -> Path:
    here = Path(__file__).resolve()
    for d in (here, *here.parents):
        if d.is_dir() and (d / "CLAUDE.md").exists():
            return d
    raise SystemExit("repo root not found")


ROOT = _root()
sys.path.insert(0, str(ROOT))
from _paths import DB_INSPECT, REFDATA  # noqa: E402

HERE = Path(__file__).resolve().parent
HOLOBASE = REFDATA / "holobase.db"
USDA_TSV = HERE / "fungi_usda.tsv"
ANNOT_TSV = HERE / "fungi_annotated.tsv"
ANNOT_CACHE = HERE / ".annotation_cache.jsonl"   # gitignored; a rerun costs nothing

# Assembly levels a reference may be built from. Contig-only taxa are excluded here because
# 01_search/kraken_search.py's own gate is scaffold_plus(), so listing them as seeds would
# promise taxa the selector then silently drops. If that gate is ever relaxed (there is an
# open decision about it for six resident genera whose only annotated assembly is a contig,
# including Papiliotrema, rice's most abundant genus), relax it here too and in step.
SCAFFOLD_PLUS = ("Scaffold", "Chromosome", "Complete Genome")

# Ranks worth seeding. formae speciales are kept deliberately: the rusts' f.sp. split is a
# real signal in 05_filter, and `_expand_taxids()` must pass intermediate_nodes=True for the
# same reason. Dropping them would lose 38 taxa including wheat rust lineages.
RANKS = ("species", "forma specialis", "varietas", "forma")


def db_v3_fungal_species() -> set[str]:
    """NCBI taxids of db_v3's fungal species nodes, from its own inspect table."""
    out, stack = set(), []
    with open(DB_INSPECT) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            raw = p[5]
            depth = (len(raw) - len(raw.lstrip(" "))) // 2
            while stack and stack[-1][0] >= depth:
                stack.pop()
            anc = {t for _, t in stack}
            stack.append((depth, p[4].strip()))
            if "4751" in anc and p[3].strip() == "S":
                out.add(p[4].strip())
    return out


# FUNGI ONLY. Leon's call 2026-10-13: this is a fungi-only paper. Oomycetes are out of scope,
# which also removes a holobase problem rather than inheriting it — holobase holds NO oomycete
# assemblies at all (its assembly sources are the genbank/refseq FUNGI and PLANT summaries and
# nothing else), so the 738 plant-associated oomycete taxa USDA records would have had to be
# resolved against NCBI separately. Dropping them costs nothing measurable: ZERO of db_v3's 29
# oomycete species pass both 05_filter floors on any cereal host. See exclusions.tsv.
#
# The kingdom filter on the FUNGUS side is not optional. Without it the list picks up plant
# parasites: USDA legitimately records Orobanche cernua var. cumana (broomrape) and Capsicum
# annuum as "recorded on" a host, and Capsicum was being selected as an addable reference with
# 6 annotated assemblies before the filter was added.
_SELECT = """
SELECT tf.ncbi_taxid                          AS ncbi_taxid,
       tf.scientific_name                     AS scientific_name,
       tf.rank                                AS rank,
       tf.kingdom                             AS kingdom,
       COUNT(DISTINCT th.ncbi_taxid)          AS n_plant_hosts,
       (SELECT COUNT(*) FROM assembly a
          WHERE a.taxon_id = tf.taxon_id
            AND a.assembly_level IN {levels}) AS n_assemblies
  FROM host_association h
  JOIN taxon tf ON tf.taxon_id = h.fungus_id
  JOIN taxon th ON th.taxon_id = h.host_id
 WHERE h.source = 'usda_fungus_host'
   AND tf.ncbi_taxid IS NOT NULL
   AND th.kingdom = 'Viridiplantae'          -- NOT 'Plants'; that returns zero
   AND tf.kingdom = '{kingdom}'
   AND tf.rank IN {ranks}
 GROUP BY tf.ncbi_taxid
{having}
 ORDER BY n_assemblies DESC, n_plant_hosts DESC
"""

QUERY_FUNGI = _SELECT.format(levels=SCAFFOLD_PLUS, ranks=RANKS, kingdom="Fungi",
                             having="HAVING n_assemblies > 0")


def _binomial(name: str) -> str:
    m = re.match(r"^([A-Z][a-z]+ [a-z][a-z-]+)", name or "")
    return m.group(1) if m else (name or "")


def fungal_genus_taxids() -> dict:
    """genus name -> NCBI taxid of the FUNGAL genus of that name.

    Needed for two kinds of name the datasets CLI refuses outright, both of which otherwise
    look like a genus with no assemblies:

      homonyms        Pichia, Porina, Rhizophagus, Rhizosphaera, Sphaerellopsis, Sugiyamaella,
                      Venturia, Volutella all match more than one taxid ("exact match for more
                      than one taxid"). Cryptococcus and Rhizophagus are the known-costly pair.
      bracketed names [Candida], [Emmonsia] — NCBI brackets a genus whose placement is
                      provisional, and the brackets are part of the scientific_name, so the
                      literal string is "not exact".

    holobase's own tree carries `kingdom`, so the fungal one can be picked without guessing.
    """
    con = sqlite3.connect(f"file:{HOLOBASE}?mode=ro", uri=True)
    rows = con.execute(
        "SELECT scientific_name, ncbi_taxid FROM taxon "
        "WHERE rank='genus' AND kingdom='Fungi' AND ncbi_taxid IS NOT NULL").fetchall()
    con.close()
    return {n: t for n, t in rows}


def _genus_assemblies(genus: str) -> list:
    """Every current assembly in a genus, from the datasets CLI. [] on failure, with a warning.

    Queried per GENUS, not per taxon: 584 calls against 1,488, and one call returns the
    species-level detail needed to match our seed list anyway.
    """
    r = subprocess.run(
        ["datasets", "summary", "genome", "taxon", genus, "--as-json-lines",
         "--limit", "all", "--assembly-version", "latest"],
        capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        # A taxonomic HOMONYM exits non-zero ("exact match for more than one taxid") and would
        # otherwise look like a genus with no assemblies. Cryptococcus and Rhizophagus both do
        # this; see exclusions.tsv.
        # The CLI writes a "New version of client available" notice to stderr too, so take
        # the first line that actually looks like an error rather than err[0].
        err = [ln for ln in (r.stderr or "").strip().splitlines()
               if ln.lstrip().startswith("Error")] or (r.stderr or "").strip().splitlines()
        print(f"\n  WARNING {genus}: {err[0] if err else 'exit ' + str(r.returncode)}",
              flush=True)
        return []
    out = []
    for line in r.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        out.extend(d["reports"]) if "reports" in d else out.append(d)
    return out


def annotate(workers: int = 6):
    """Ask NCBI which seed taxa have a gene model, and so a usable CDS."""
    seeds = list(csv.DictReader(open(USDA_TSV), delimiter="\t"))
    if not seeds:
        raise SystemExit(f"{USDA_TSV} is empty; run without --annotate first")
    by_taxid = {r["ncbi_taxid"]: r for r in seeds}
    by_binom = defaultdict(list)
    for r in seeds:
        by_binom[_binomial(r["scientific_name"])].append(r)
    genera = sorted({r["scientific_name"].split()[0] for r in seeds})

    done = {}
    if ANNOT_CACHE.exists():
        for line in open(ANNOT_CACHE):
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            done[d["genus"]] = d["rows"]
        print(f"cache: {len(done)}/{len(genera)} genera already queried")

    todo = [g for g in genera if g not in done]
    print(f"querying {len(todo)} genera at {workers} workers …")

    genus_taxid = fungal_genus_taxids()

    def work(genus):
        # ACCESSIONS, not counts: an assembly must be attributed to at most one seed, and the
        # download size is a set union rather than a sum.
        hits = defaultdict(lambda: {"sp": set(), "ann": set()})
        assemblies = _genus_assemblies(genus)
        if not assemblies and genus in genus_taxid:
            # Homonym or bracketed name: retry by taxid, which is unambiguous.
            assemblies = _genus_assemblies(str(genus_taxid[genus]))
            if assemblies:
                print(f"\n  recovered {genus} via taxid {genus_taxid[genus]} "
                      f"({len(assemblies)} assemblies)", flush=True)
        for a in assemblies:
            info = a.get("assembly_info") or {}
            org = a.get("organism") or {}
            acc = a.get("accession")
            if not acc or (info.get("assembly_level") or "") not in SCAFFOLD_PLUS:
                continue
            annotated = bool(a.get("annotation_info"))
            # Match a seed by its own taxid first: 38 seeds are formae speciales or varietas
            # whose names are not binomials, so a binomial collapse would merge them into the
            # parent species.
            key = str(org.get("tax_id")) if str(org.get("tax_id")) in by_taxid else None
            if key:
                targets = [by_taxid[key]]
            else:
                cand = by_binom.get(_binomial(org.get("organism_name") or ""), [])
                # If the binomial maps to SEVERAL seeds, they are formae speciales of one
                # species and this assembly sits at the species node, so it cannot be
                # attributed to any one of them. Crediting all of them is what inflated
                # Fusarium oxysporum: every one of its 10 f.sp. seeds was handed the same
                # ~129 assemblies, turning ~130 real assemblies into ~1,290.
                targets = cand if len(cand) == 1 else []
            for t in targets:
                h = hits[t["ncbi_taxid"]]
                h["sp"].add(acc)
                if annotated:
                    h["ann"].add(acc)
        return genus, {k: {"sp": sorted(v["sp"]), "ann": sorted(v["ann"])}
                       for k, v in hits.items()}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work, g): g for g in todo}
        with open(ANNOT_CACHE, "a") as cache:
            for n, fut in enumerate(as_completed(futs), 1):
                genus, rows = fut.result()
                done[genus] = rows
                cache.write(json.dumps({"genus": genus, "rows": rows}) + "\n")
                cache.flush()
                if n % 20 == 0 or n == len(todo):
                    print(f"  {n}/{len(todo)} genera …", end="\r", flush=True)
    print()

    merged = defaultdict(lambda: {"sp": set(), "ann": set()})
    for rows in done.values():
        for tid, v in rows.items():
            # Tolerate the pre-fix cache shape ([sp, ann] counts) by ignoring it: a stale
            # cache entry would silently reintroduce the fan-out bug.
            if not isinstance(v, dict):
                continue
            merged[tid]["sp"].update(v["sp"])
            merged[tid]["ann"].update(v["ann"])

    out = []
    for r in seeds:
        m = merged.get(r["ncbi_taxid"], {"sp": set(), "ann": set()})
        sp, ann = len(m["sp"]), len(m["ann"])
        out.append({
            "ncbi_taxid": r["ncbi_taxid"], "scientific_name": r["scientific_name"],
            "rank": r["rank"], "in_db_v3": r["in_db_v3"],
            "n_scaffold_plus": sp, "n_annotated": ann,
            # The ONLY thing that matters for a CDS-only database.
            "route": ("present" if r["in_db_v3"] == "yes"
                      else "add" if ann > 0
                      else "miniprot" if sp > 0
                      else "none"),
        })
    out.sort(key=lambda r: (r["route"], -r["n_annotated"], r["scientific_name"]))
    with open(ANNOT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "ncbi_taxid", "scientific_name", "rank", "in_db_v3",
            "n_scaffold_plus", "n_annotated", "route"])
        w.writeheader()
        w.writerows(out)

    print(f"\n{'route':10s} {'taxa':>6s} {'annotated assemblies':>21s}")
    for route in ("present", "add", "miniprot", "none"):
        sel = [r for r in out if r["route"] == route]
        if sel:
            print(f"{route:10s} {len(sel):6d} {sum(r['n_annotated'] for r in sel):21d}")
    add = [r for r in out if r["route"] == "add"]
    uniq = set()
    for r in add:
        uniq |= merged.get(r["ncbi_taxid"], {"ann": set()})["ann"]
    print(f"\nADDABLE NOW: {len(add)} taxa, {len(uniq)} UNIQUE annotated assemblies "
          f"(pre-floor; the per-taxon column sums higher where seeds share an assembly)")
    print(f"wrote {ANNOT_TSV.relative_to(ROOT)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    ap.add_argument("--annotate", action="store_true",
                    help="second pass: ask NCBI which seeds have a gene model (needs network)")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()

    if a.annotate:
        annotate(a.workers)
        return

    if not HOLOBASE.exists():
        raise SystemExit(f"holobase not found at {HOLOBASE}")
    con = sqlite3.connect(f"file:{HOLOBASE}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    rows = [dict(r) for r in con.execute(QUERY_FUNGI)]
    con.close()

    v3 = db_v3_fungal_species()
    for r in rows:
        r["in_db_v3"] = "yes" if str(r["ncbi_taxid"]) in v3 else "no"

    new = [r for r in rows if r["in_db_v3"] == "no"]
    print(f"USDA plant-associated fungi with a scaffold+ assembly : {len(rows):,}")
    print(f"  already in db_v3                                   : {len(rows) - len(new):,}")
    print(f"  NEW to db_v3                                       : {len(new):,}")
    print(f"  assemblies they hold (pre-floor, pre-annotation)    : "
          f"{sum(r['n_assemblies'] for r in new):,}")
    print("\nNOTE: n_assemblies ignores annotation, which holobase cannot report "
          "(is_annotated is 1 on\nevery row). The database is CDS-only, so the real "
          "addition needs a --select-only run.")

    if a.check:
        print("\n--check: nothing written")
        return
    with open(USDA_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "ncbi_taxid", "scientific_name", "rank", "kingdom", "n_plant_hosts",
            "n_assemblies", "in_db_v3"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {USDA_TSV.relative_to(ROOT)}  ({len(rows):,} rows)")


if __name__ == "__main__":
    main()
