#!/usr/bin/env python3
"""fetch_cds.py — put a CDS FASTA on disk for every accession in ref_candidates.tsv.

Step 1 of 3 in 02_build (fetch -> saturate -> build). The only step that moves bulk bytes,
and the only one that may talk to NCBI.

**Three sources, tried in cost order, and the order is the whole point.** Of the 4,286
accessions db_v4 wants, 1,753 were already downloaded for db_v3 and are sitting in Acacia's
project bucket. Pulling those from Acacia is free and internal; re-pulling them from NCBI is
neither, and NCBI is also the one source that can rate-limit, truncate or silently hand back
a different assembly than the one requested. So:

    1. DISK    already extracted under bronze/cds/pathogen/<accession>/ -> nothing to do
    2. ACACIA  in s3://<bucket>/db/search/cds/pathogen/<accession>/ -> aws s3 cp
    3. NCBI    everything left -> `datasets download genome accession`

## Presence is a non-empty file, never a directory

This is the one rule the step is built around, because the 03_kraken tree violated it. On
2026-10-16 Setonix scratch held **1,997 accession directories but only 1,162 non-empty CDS
FASTAs**: the other 835 were directory skeletons left behind when /scratch's 21-day purge
took the files inside them. Every `if dir.exists()` check counts those 835 as held, finds no
sequence at build time, and the build quietly proceeds 835 assemblies short — which is
exactly how db_v3 came to be built from 2,081 of 2,115 selected assemblies with the shortfall
recorded only as a line in a log nobody re-read.

So `have_cds()` tests for a non-empty `.fna` and nothing else, and the counts printed at the
end must add up to the wanted total or the step exits non-zero.

## It fetches everything, and there is no per-species ceiling

Every accession in ref_candidates.tsv is fetched, including the 294 *Aspergillus fumigatus*
genomes. A per-taxon assembly cap is **already a rejected decision** (CLAUDE.md, "Decisions
worth not relitigating"): it cannot buy breadth and it costs depth in exactly the well-sampled
taxa where depth was measured to matter. Re-measured 2026-10-16 from
`03_kraken/utilities/saturation/`, which agrees:

  * A flat cap loses real k-mer space, unevenly. Capping at one assembly per species retains
    **49%** of the profiled pangenome k-mer space; at five, 73%. The spread is the problem, not
    the mean: `pangenome_ratio` runs from 1.02 (*Fulvia fulva*) to 14.34 (*Rhizoctonia solani*),
    so at a cap of one, *R. solani* keeps 7% of its k-mer space and *Verticillium nonalfalfae*
    keeps 87%. No single number is right for both.
  * The measured per-species rule — cap at `n_for_95pct`, the N reaching 95% of that species'
    own pangenome — keeps 610 of 857 assemblies, 71%. The pangenomes are open (Heaps gamma
    0.16), so there is no cheap ceiling to find.
  * Depth is nearly free in accuracy terms. `depth_bias.tsv` measures a deep species stealing a
    shallow congener's reads across 3,143 pairs: median misassignment 0.065% at
    `--minimum-hit-groups 1` and 0.001% at 4, against a production setting of 3. It grows
    monotonically with depth but only about 2x from depth 1 to 12.

So depth costs index size and bandwidth, not correctness, and both are affordable: ~50 GB to
fetch and an index that still fits one 250 GB node. Saturation stays a QC RECORD, the same
status BUSCO has and for the same reason — a screen that can silently drop a taxon is a
liability, and this one would drop exactly the well-sampled taxa we most rely on.

Run on Setonix, from the repo root. Not on a login node: this moves ~50 GB.

    sbatch 02_classify/slurm/fetch_cds.slurm
    python 02_classify/02_build/fetch_cds.py --inventory        # report only, moves nothing
"""

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path

MODULE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import (ACACIA_BUCKET, ACACIA_CDS_PREFIX, ACACIA_ENDPOINT,  # noqa: E402
                     CDS_MANIFEST, CDS_PATHOGEN, LOGS, REF_CANDIDATES, SILVER)
from _util import _Tee, atomic_write, link_latest, make_log_dir  # noqa: E402

# A CDS FASTA under this size is not a fungal gene set; the smallest real one measured across
# the 1,019 backed-up assemblies is 1.6 MB. 100 kB is a floor for "something arrived but it is
# not what we asked for" (an error page, a truncated transfer, a genome with 40 genes).
MIN_FASTA_BYTES = 100_000

# Names written by downstream steps into the same directory, never the download itself.
DERIVED_SUFFIXES = ("_clean.fna", "_combined.fna", ".tagged.fna")


def cds_files(acc_dir: Path) -> list:
    """Every non-derived .fna under an accession directory, largest first."""
    if not acc_dir.is_dir():
        return []
    out = [f for f in acc_dir.glob("**/*.fna")
           if not f.name.endswith(DERIVED_SUFFIXES) and f.is_file()]
    return sorted(out, key=lambda f: -f.stat().st_size)


def have_cds(acc_dir: Path) -> Path | None:
    """The accession's CDS FASTA, or None. Tests for a NON-EMPTY FILE, not a directory.

    See the module docstring: 835 of 1,997 directories on scratch were empty skeletons left
    by the purge, and a directory check calls all of them present.
    """
    for f in cds_files(acc_dir):
        if f.stat().st_size >= MIN_FASTA_BYTES:
            return f
    return None


def load_wanted(path: Path) -> list:
    """The accessions db_v4 wants, from ref_candidates.tsv.

    `fasta_type` is tested POSITIVELY for "cds". The negative form (`!= "genome"`) silently
    admits a blank, which is how an unannotated assembly reached a CDS-only database before.
    """
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("fasta_type") or "").strip() == "cds":
                rows.append(r)
    if not rows:
        sys.exit(f"no cds rows in {path} — wrong table, or 01_refs has not run")
    seen, uniq = set(), []
    for r in rows:
        if r["accession"] not in seen:
            seen.add(r["accession"])
            uniq.append(r)
    if len(uniq) != len(rows):
        print(f"  note: {len(rows) - len(uniq):,} duplicate accession rows collapsed")
    return uniq


def acacia_inventory(bucket: str, prefix: str, profile: str) -> dict:
    """{accession: bytes} for every non-empty CDS FASTA in the bucket.

    One recursive listing, not one HEAD per accession: 4,286 round trips would dominate the
    runtime of a step whose actual work is a bulk copy.
    """
    cmd = ["aws", "s3", "ls", "--recursive", "--endpoint-url", ACACIA_ENDPOINT,
           "--profile", profile, f"s3://{bucket}/{prefix}/"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        err = (r.stderr or "").strip()
        hint = ""
        if "NoneType" in err:
            hint = ("\n  This is aws-cli 2.13 crashing while formatting a Ceph denial with an "
                    "empty <Message>. It means the profile cannot read this bucket, not a "
                    "network fault. Check the profile owns the PROJECT allocation.")
        elif "Unable to locate credentials" in err:
            hint = ("\n  On Setonix the CLI is a Singularity wrapper that overrides --home: it "
                    "reads $MYSOFTWARE/setonix/<ver>/.shpc_home/.aws, not ~/.aws.")
        sys.exit(f"cannot list s3://{bucket}/{prefix}: {err}{hint}")

    out = {}
    for line in r.stdout.splitlines():
        f = line.split()
        if len(f) < 4:
            continue
        size, key = int(f[2]), f[-1]
        name = key.rsplit("/", 1)[-1]
        if not name.endswith(".fna") or name.endswith(DERIVED_SUFFIXES):
            continue
        if size < MIN_FASTA_BYTES:
            continue
        parts = key.split("/")
        if len(parts) < 2:
            continue
        acc = parts[-2]
        if size > out.get(acc, 0):
            out[acc] = size
    return out


def pull_from_acacia(acc: str, dest: Path, bucket: str, profile: str) -> Path | None:
    """Copy one accession's CDS tree out of Acacia. Returns the FASTA, or None."""
    dest.mkdir(parents=True, exist_ok=True)
    cmd = ["aws", "s3", "cp", "--recursive", "--endpoint-url", ACACIA_ENDPOINT,
           "--profile", profile, "--only-show-errors",
           f"s3://{bucket}/{ACACIA_CDS_PREFIX}/{acc}/", str(dest)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        print(f"  {acc}  acacia pull FAILED: {(r.stderr or '').strip()[:160]}", flush=True)
        return None
    return have_cds(dest)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--candidates", default=None,
                    help=f"input table (default {REF_CANDIDATES})")
    ap.add_argument("--cds-dir", default=None,
                    help=f"destination tree (default {CDS_PATHOGEN})")
    ap.add_argument("--bucket", default=ACACIA_BUCKET,
                    help=f"Acacia bucket (default {ACACIA_BUCKET}; the PROJECT one)")
    ap.add_argument("--profile", default="pawsey1168",
                    help="aws profile; Setonix uses pawsey1168, which owns the project bucket")
    ap.add_argument("--workers", type=int, default=8,
                    help="parallel NCBI downloads. 8 measured faster than 16: beyond ~8 the "
                         "NCBI throttles and retries cost more than the concurrency gains")
    ap.add_argument("--inventory", action="store_true",
                    help="report what each source holds and exit; moves no bytes")
    ap.add_argument("--no-acacia", action="store_true",
                    help="skip the Acacia pass (for a deliberate clean re-fetch from NCBI)")
    ap.add_argument("--no-ncbi", action="store_true",
                    help="reconcile and pull from Acacia only; leaves the NCBI gap unfetched")
    args = ap.parse_args()

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "fetch_cds.log")
    link_latest(LOGS, log_dir)
    try:
        cand_path = Path(args.candidates) if args.candidates else REF_CANDIDATES
        cds_dir = Path(args.cds_dir) if args.cds_dir else CDS_PATHOGEN
        wanted = load_wanted(cand_path)
        print(f"wanted: {len(wanted):,} cds accessions, "
              f"{len({r['organism_name'] for r in wanted}):,} species  ({cand_path})")

        # ── 1. disk ─────────────────────────────────────────────────────────────────────
        on_disk = {}
        for r in wanted:
            f = have_cds(cds_dir / r["accession"])
            if f:
                on_disk[r["accession"]] = f
        print(f"on disk already:        {len(on_disk):,}   ({cds_dir})")

        # Count the skeletons explicitly rather than letting them hide inside "missing".
        skeletons = [r["accession"] for r in wanted
                     if r["accession"] not in on_disk
                     and (cds_dir / r["accession"]).is_dir()]
        if skeletons:
            print(f"  of which EMPTY directories (purged, not held): {len(skeletons):,}")

        # ── 2. acacia ───────────────────────────────────────────────────────────────────
        in_acacia = {}
        if not args.no_acacia:
            inv = acacia_inventory(args.bucket, ACACIA_CDS_PREFIX, args.profile)
            print(f"acacia bucket holds:    {len(inv):,}   (s3://{args.bucket}/{ACACIA_CDS_PREFIX})")
            in_acacia = {r["accession"]: inv[r["accession"]] for r in wanted
                         if r["accession"] in inv and r["accession"] not in on_disk}
            print(f"  pullable from acacia: {len(in_acacia):,}  "
                  f"({sum(in_acacia.values())/1e9:.1f} GB)")
            orphans = len(inv) - len([a for a in inv if a in {r['accession'] for r in wanted}])
            print(f"  in acacia, no longer wanted: {orphans:,}  (kept, never deleted here)")

        gap = [r for r in wanted
               if r["accession"] not in on_disk and r["accession"] not in in_acacia]
        print(f"must come from NCBI:    {len(gap):,}")

        if args.inventory:
            print("\n--inventory: nothing moved.")
            return

        # ── pull ────────────────────────────────────────────────────────────────────────
        source = {a: "disk" for a in on_disk}
        if in_acacia:
            print(f"\nPulling {len(in_acacia):,} from Acacia →  {cds_dir}")
            ok = 0
            for i, acc in enumerate(sorted(in_acacia), 1):
                f = pull_from_acacia(acc, cds_dir / acc, args.bucket, args.profile)
                if f:
                    on_disk[acc] = f
                    source[acc] = "acacia"
                    ok += 1
                else:
                    gap.append(next(r for r in wanted if r["accession"] == acc))
                if i % 200 == 0 or i == len(in_acacia):
                    print(f"  [{i:,}/{len(in_acacia):,}]  {ok:,} ok", flush=True)
            print(f"  acacia: {ok:,} of {len(in_acacia):,}")

        if gap and not args.no_ncbi:
            from _assembly import run_download  # noqa: E402  imported late: NCBI-only path
            print(f"\nDownloading {len(gap):,} from NCBI ({args.workers} workers)")
            run_download(gap, cds_dir, args.workers)
            for r in gap:
                f = have_cds(cds_dir / r["accession"])
                if f:
                    on_disk[r["accession"]] = f
                    source[r["accession"]] = "ncbi"
        elif gap:
            print(f"\n--no-ncbi: {len(gap):,} accessions left unfetched.")

        # ── manifest ────────────────────────────────────────────────────────────────────
        SILVER.mkdir(parents=True, exist_ok=True)
        cols = ["accession", "taxid", "organism_name", "routes", "source", "bytes", "path"]
        missing = []
        with atomic_write(CDS_MANIFEST) as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(cols)
            for r in wanted:
                acc = r["accession"]
                f = on_disk.get(acc)
                if not f:
                    missing.append(r)
                w.writerow([acc, r["taxid"], r["organism_name"], r.get("routes", ""),
                            source.get(acc, "MISSING"),
                            f.stat().st_size if f else 0,
                            str(f.relative_to(cds_dir)) if f else ""])
        print(f"\nmanifest → {CDS_MANIFEST}")

        tally = Counter(source.get(r["accession"], "MISSING") for r in wanted)
        print("\n── summary ──────────────────────────────────────────────")
        for k in ("disk", "acacia", "ncbi", "MISSING"):
            if tally.get(k):
                print(f"  {k:8s} {tally[k]:6,d}")
        held = sum(v for k, v in tally.items() if k != "MISSING")
        print(f"  {'total':8s} {held:6,d} of {len(wanted):,} wanted  "
              f"({sum(f.stat().st_size for f in on_disk.values())/1e9:.1f} GB)")

        if missing:
            path = SILVER / "cds_missing.tsv"
            with atomic_write(path) as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["accession", "taxid", "organism_name", "routes"])
                for r in missing:
                    w.writerow([r["accession"], r["taxid"], r["organism_name"],
                                r.get("routes", "")])
            print(f"\n{len(missing):,} accessions have NO CDS on disk → {path}")
            print("  The build must not proceed as if these were present. Re-run to retry; "
                  "a persistent failure is a real NCBI withdrawal and belongs in "
                  "seeds/exclusions.tsv with the evidence.")
            sys.exit(1)
    finally:
        log.close()


if __name__ == "__main__":
    main()
