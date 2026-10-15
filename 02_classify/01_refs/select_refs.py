#!/usr/bin/env python3
"""select_refs.py — db_v4's reference set: every annotated assembly the survey found.

    python 02_classify/01_refs/select_refs.py --select-only     # write the candidate table
    python 02_classify/01_refs/select_refs.py --download        # ... then fetch the CDS

Step 1c of 02_classify, after `routes.py` (who is a member, and why) and `survey.py` (what
exists for them at NCBI). It queries nothing: it filters `data/silver/assemblies.tsv`, which
is the survey's own record. That is deliberate. NCBI's catalogue moves, so a selector that
re-queried could disagree with the survey that justified it; filtering a file instead makes the
selection reproducible and the disagreement impossible.

## The rule, in full

**Take every assembly that carries a gene model. No floor, no cap, no quality gate.**

That is the whole policy, and it replaces three mechanisms that each needed defending:

  * **The depth floor and its no-cap exception are gone.** db_v3 took `min(floor, available)`
    for most species and unlimited depth for the 38 genera STAT had detected, which made STAT a
    gate on database COMPOSITION, and it failed permissive when `01_stat` moved away: the
    loader read a missing file as "all genera detected", lifting the cap everywhere. Taking
    everything is simpler AND removes the gate. It is affordable for the reason the floor was
    never the real constraint on breadth: the median species has ONE assembly, so the floor
    only ever clipped a handful of deep ones, and because pangenomes are open with Heaps gamma
    ~0.16 the extra assemblies on a deep species add sublinear k-mer content.
  * **The assembly-level gate is gone.** db_v3 required `scaffold_plus AND has_annotation`, and
    for six resident genera those never co-occur, which excluded *Papiliotrema* (rice's most
    abundant genus at 14.5%) and *Dioszegia* (a wheat core yeast). A reference exists to give a
    read somewhere nearer to park, and that needs real k-mers rather than a good assembly.
    Reference incompleteness largely cancels in the coverage ratio anyway, because a fragmented
    assembly shrinks numerator and denominator together.
  * **BUSCO is not a gate and never was.** At the 50% threshold 2 of 1,018 assemblies fail, and
    a per-taxid fallback meant 0 of 309 taxa were ever lost to it. It stays as a SCORE, because
    it is the tie-breaker deciding which single assembly represents a taxon, which `03_align`'s
    kallisto index depends on. Do not filter on duplication either: the worst offenders are
    *P. graminis* at 86.6% duplicated and *V. longisporum* at 71.3%, which is dikaryotic and
    allodiploid biology, and such a filter would preferentially remove the rusts whose absence
    caused the *Phakopsora* artefact. Kraken2 stores each distinct minimizer once, so
    duplicated sequence adds no weight regardless.

**Annotation is the one hard requirement, and it is not ours.** `kraken_build.py` skips any
assembly with no gene model, so an unannotated assembly contributes exactly nothing to a
CDS-only index. Those species are not rejected, they are deferred: `survey.py` writes them to
`miniprot_worklist.tsv`, and gene prediction is where the remaining breadth lives.

## What depth means now, and the caveat that survives

Taking everything does not EQUALISE depth, it maximises it. Depth then tracks research
attention rather than biology, so the asymmetry that drives misassignment between congeners is
still there, just set by NCBI's sampling instead of by our floor. The reporting rule therefore
stands unchanged: **species-level counts for shallow taxa are not reliable and should be
aggregated to genus**, and misassignment risk is highest where a shallow species sits beside a
deeply-sampled congener.

## The route column travels with the assembly

Every row carries the routes that made its species a member. `05_detect` must read it, because
a taxon present ONLY as a `congener` is in the database to absorb reads, and reporting a
detection on it the way a target is reported would be an artefact of our own database design.
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _layout import SILVER, CDS_PATHOGEN, REF_CANDIDATES, LOGS   # noqa: E402
from _util import _Tee, make_log_dir, link_latest                # noqa: E402
import _assembly as A                                            # noqa: E402

ASSEMBLIES_TSV = SILVER / "assemblies.tsv"
OUT_COLS = A.CANDIDATE_COLS + ["routes"]


def load_assemblies(path: Path) -> list:
    if not path.exists():
        raise SystemExit(
            f"{path} not found.\n"
            f"Run the two steps that produce it first:\n"
            f"  python 02_classify/01_refs/routes.py --write\n"
            f"  python 02_classify/01_refs/survey.py")
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"{path} has no rows")
    for col in ("accession", "fasta_type", "taxid", "routes"):
        if col not in rows[0]:
            raise SystemExit(f"{path}: missing `{col}`. Found {list(rows[0])}")
    return rows


def select(rows: list, require_scaffold: bool = False) -> tuple:
    """Keep every annotated assembly. Returns (kept, {reason: count}) for what was not kept."""
    kept, dropped = [], defaultdict(int)
    seen = set()
    for r in rows:
        acc = r["accession"]
        if not acc or acc in seen:
            dropped["duplicate_accession"] += 1
            continue
        # fasta_type is `cds` exactly when the assembly carries a gene model, which is the only
        # condition kraken_build.py can act on. Tested POSITIVELY: an absent or unexpected value
        # must fall through to dropped, never into the build. The same discipline as the NaN
        # score floors in 05_detect, where `if x < FLOOR: continue` admitted every NaN row.
        if r.get("fasta_type") != "cds":
            dropped["no_gene_model"] += 1
            continue
        if require_scaffold and r.get("assembly_level") in ("Contig", ""):
            dropped["below_scaffold_strict_gate"] += 1
            continue
        seen.add(acc)
        kept.append(r)
    return kept, dict(dropped)


def report(kept: list, dropped: dict, rows: list) -> None:
    by_species, by_route, by_level = defaultdict(int), defaultdict(int), defaultdict(int)
    for r in kept:
        by_species[r["taxid"]] += 1
        by_level[r.get("assembly_level") or "?"] += 1
        for route in (r.get("routes") or "via_clade").split(";"):
            by_route[route] += 1

    print(f"\n-- db_v4 reference set -----------------------------------------")
    print(f"  assemblies taken    {len(kept):>8,}   of {len(rows):,} surveyed")
    print(f"  species             {len(by_species):>8,}")
    print(f"\n  not taken:")
    for reason, n in sorted(dropped.items(), key=lambda kv: -kv[1]):
        print(f"    {reason:<28} {n:>8,}")
    print(f"\n  by assembly level:")
    for lvl, n in sorted(by_level.items(), key=lambda kv: -kv[1]):
        print(f"    {lvl:<28} {n:>8,}")
    print(f"\n  by membership route (a species may be in by several):")
    for route, n in sorted(by_route.items(), key=lambda kv: -kv[1]):
        print(f"    {route:<28} {n:>8,}")

    depth = defaultdict(int)
    for n in by_species.values():
        key = ("1" if n == 1 else "2-3" if n < 4 else "4-9" if n < 10
               else "10-49" if n < 50 else "50+")
        depth[key] += 1
    print(f"\n  depth per species (no floor, no cap):")
    for key in ("1", "2-3", "4-9", "10-49", "50+"):
        if depth[key]:
            print(f"    {key:>6} assemblies: {depth[key]:>6,} species")
    print(f"\n  db_v3 for scale: 2,063 assemblies / 877 fungal species / 4.80 B minimizers")
    print(f"  Depth tracks research attention, not biology. Aggregate shallow species to")
    print(f"  genus when reporting; taking everything does not change that caveat.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--strict-gate", action="store_true",
                    help="also require scaffold-plus, as db_v3 did. For comparison only; it "
                         "costs Papiliotrema and Dioszegia, whose only annotated assembly is "
                         "contig-level")
    ap.add_argument("--in", dest="infile", default=None,
                    help=f"survey output to filter (default {ASSEMBLIES_TSV})")
    ap.add_argument("--out", default=None, help=f"candidate table (default {REF_CANDIDATES})")
    ap.add_argument("--workers", type=int, default=8,
                    help="download concurrency. 8, not 16: 16 measured worse for NCBI")
    ap.add_argument("--select-only", action="store_true",
                    help="write the table, download nothing")
    ap.add_argument("--download", action="store_true",
                    help="after selecting, download CDS for every candidate. Resumable: "
                         "skips accessions already on disk")
    args = ap.parse_args()
    if not (args.select_only or args.download):
        ap.error("pass --select-only or --download")

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "select_refs.log")
    link_latest(LOGS, log_dir / "select_refs.log")
    sys.stdout = log
    try:
        src = Path(args.infile) if args.infile else ASSEMBLIES_TSV
        rows = load_assemblies(src)
        print(f"  read {len(rows):,} surveyed assemblies from {src}")
        kept, dropped = select(rows, args.strict_gate)
        report(kept, dropped, rows)

        out = Path(args.out) if args.out else REF_CANDIDATES
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(kept)
        print(f"\n  -> {out}")

        if args.download:
            # Cross-reference Acacia before paying NCBI for bytes we already own. The CDS tree
            # there is keyed by ACCESSION under db/search (5,641 objects / 40.1 GB, cds/ only),
            # so anything already backed up needs no refetch. download_cds is resumable and
            # skips what is on disk, but on a purged scratch the Acacia copy has to come down
            # first. The --endpoint-url is MANDATORY: without it the request goes to Amazon,
            # which returns InvalidAccessKeyId and reads exactly like expired Pawsey keys.
            print(f"\n  NOTE check Acacia first:")
            print(f"    module load aws-cli/2.13.0")
            print(f"    aws s3 ls s3://pawsey1168-kraken-db/db/search/cds/ \\")
            print(f"        --profile pawsey1168 --endpoint-url https://projects.pawsey.org.au")
            A.run_download(kept, CDS_PATHOGEN, args.workers)
    finally:
        log.close()


if __name__ == "__main__":
    main()
