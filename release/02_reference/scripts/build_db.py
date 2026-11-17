#!/usr/bin/env python3
"""Build one Kraken2 index from the CDS the fetch steps put on disk.

    sbatch release/02_reference/scripts/slurm/build_db.slurm        # NAME=db1p
    python release/02_reference/scripts/build_db.py --name db1p \
        --tier db0_baseline,db1_congener --dry-run

Reads  data/curation/build_plan.tsv        which assemblies, which tier, which taxid
       data/raw/cds_fetch_manifest.tsv     what the fungal fetch actually landed
       data/raw/host_cds_manifest.tsv      the resolved Ensembl CDS filenames
       data/raw/raw_taxdmp.zip             the PINNED taxonomy, for db/taxonomy/
Writes data/db/<name>/                     the index (gitignored; Setonix scratch)
       data/outputs/db_inspect_<name>.tsv  minimizers per taxon, TRACKED

Fetches nothing. Chooses nothing. The plan decides what goes in and this puts it there, which
is the same split the module already runs on everywhere else: a screen decides, a fetch acts.
Keeping them apart is not tidiness. `build_plan.tsv` has to be IDENTICAL across the arms of
the ladder for db1'-db0 and db2-db1' to isolate one effect each, and a plan regenerated as a
side effect of a six-hour build can differ between arms with nothing saying so. It has gone
stale in git once already.

## THE PLAN IS THE WHOLE INPUT, AND IT CARRIES THE RIGHT TAXID

Every row gives `species_taxid` and where its CDS comes from, so nothing here consults the
research tree:

    cds_source=ncbi            -> data/raw/cds/<accession>/cds_from_genomic.fna
    cds_source=ensembl_plants  -> data/raw/host_cds/<cds_ref>/<file from the host manifest>

`species_taxid` is the COHORT's taxid, which matters for exactly one row and would be a silent
defect: Ensembl files barley under 112509 (*H. vulgare* subsp. *vulgare*) while the cohort
records 4513, and a host call has to resolve at the rank the chapter reports. The plan already
holds 4513. A build driven from an Ensembl-derived table would tag it 112509 and the host call
would land one node below where every downstream table looks for it.

**The host tier is in EVERY database and is never subtracted.** db0 and db1 were built with no
plant CDS at all (`Fungi 100%` in their inspect tables), which is a regression: these are plant
RNA-seq libraries, so without a sink a plant read is absorbed by the nearest fungal CDS and
reported as a fungal hit. Hosts therefore go in regardless of `--tier`, so that a tier
subtraction measures the tier and not the sink.

## FOUR THINGS THIS KEEPS FROM THE 03_kraken BUILD, EACH PAID FOR

1. **A failure is FATAL, never a warning.** db_v3 printed `Skipped: 32 / Missing: 2` and built
   anyway, and nothing downstream ever learned the index was 34 assemblies short; the shortfall
   survived only as a line in a log. An index that is quietly short is worse than no index,
   because every coverage denominator from it is wrong in an unrecorded direction.
2. **The download is never mutated.** The old tagger wrote a `.tagged.fna` then renamed it over
   the original, so a rebuild could not verify what NCBI sent. Here the tagged copy is derived:
   written to a staging directory, handed to `--add-to-library` (which copies it into
   `<db>/library/` anyway), then deleted.
3. **Re-tagging keys on the tag's VALUE, not its presence.** The old tagger read the first 30
   bytes, saw `>kraken:taxid|`, and returned unchanged, so a taxid written for a previous
   database survived into the next one. Most such tags are right, which is what makes it
   dangerous. Here the existing tag is parsed, compared to the plan, replaced on any
   difference, and the count is reported.
4. **k=35, l=31, stated rather than inherited.** The detection thresholds are calibrated
   against these, so a change invalidates every coverage number in the chapter. `--no-masking`
   is deliberate too: Kraken2's LCA already resolves a shared k-mer to the common ancestor, so
   low-complexity masking removes real sequence for nothing.

## The taxonomy comes from the module's pinned taxdump

`raw/raw_taxdmp.zip` is the archive every screen in this module already reads, md5-pinned in
`raw_sources_manifest.tsv`. The 03_kraken build downloaded `taxdump.tar.gz` fresh at build
time instead, which means the index and the tables describing it could disagree about the
tree. Only `nodes.dmp` and `names.dmp` are extracted, which is all kraken2-build reads.

## db2

db2 = db1' + the 135 gene-predicted CDS from `research/predict`. Those are NOT in the plan
(plan_build.py only emits CDS that exists to download), so they enter through `--extra-cds`,
a three-column TSV of accession, taxid and path. Predicted CDS is a COMPETITOR in the
database and never a positive call, which is a reporting rule and not something this script
can enforce.
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
PLAN = HERE / "data/curation/build_plan.tsv"
CDS = HERE / "data/raw/cds"
HOST_CDS = HERE / "data/raw/host_cds"
HOST_MANIFEST = HERE / "data/raw/host_cds_manifest.tsv"
CDS_MANIFEST = HERE / "data/raw/cds_fetch_manifest.tsv"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
DB_ROOT = HERE / "data/db"
OUTPUTS = HERE / "data/outputs"

TAG_PREFIX = "kraken:taxid|"
KMER_LEN = 35
MINIMIZER_LEN = 31
HOST_TIER = "host"


@contextmanager
def atomic_write(path):
    """Write via a temp file in the same directory, rename on success.

    A killed job leaves the previous file intact rather than a half-written one, which matters
    for db_inspect: a truncated coverage denominator is read as a real one.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name + ".")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            yield fh
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def rd(path):
    with open(path) as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                   delimiter="\t"))


def host_files():
    """ensembl_species -> the resolved CDS file, from the host manifest.

    The filename cannot be constructed: Ensembl's species table and the FASTA disagree on the
    assembly token (Landmark is PGSBv2.0_Landmark in the table and PGSBv2.1 in the file), which
    is why fetch_host_cds.py records what it actually downloaded.
    """
    if not HOST_MANIFEST.exists():
        return {}
    out = {}
    for r in rd(HOST_MANIFEST):
        if r.get("file"):
            out[r["ensembl_species"]] = HOST_CDS / r["ensembl_species"] / r["file"]
    return out


def plan_rows(tiers, extra):
    """The accessions to build, as (accession, taxid, source_path, label).

    Hosts are added whatever `tiers` says; see the docstring.
    """
    if not PLAN.exists():
        sys.exit(f"no plan at {PLAN} — run plan_build.py first")
    hosts = host_files()
    rows, missing = [], []
    for r in rd(PLAN):
        tier = r["db_tier"]
        if tier != HOST_TIER and tier not in tiers:
            continue
        taxid = (r.get("species_taxid") or "").strip()
        if not taxid.isdigit():
            sys.exit(f"{r['accession']}: no usable species_taxid in the plan")
        if r.get("cds_source") == "ensembl_plants":
            src = hosts.get(r.get("cds_ref", ""))
            if src is None:
                missing.append((r["accession"], tier,
                                f"no host manifest row for {r.get('cds_ref')!r}"))
                continue
        else:
            src = CDS / r["accession"] / "cds_from_genomic.fna"
        if not src.is_file() or not src.stat().st_size:
            missing.append((r["accession"], tier, f"no CDS at {src}"))
            continue
        rows.append((r["accession"], int(taxid), src, r.get("organism_name", "")))
    for acc, taxid, path in extra:
        p = Path(path)
        if not p.is_file() or not p.stat().st_size:
            missing.append((acc, "extra", f"no CDS at {p}"))
            continue
        rows.append((acc, int(taxid), p, "(predicted)"))
    return rows, missing


def read_extra(path):
    if not path:
        return []
    out = []
    for r in rd(path):
        out.append((r["accession"], r["taxid"], r["path"]))
    return out


def opener(path):
    return gzip.open(path, "rt", errors="replace") if path.suffix == ".gz" \
        else open(path, "r", errors="replace")


def existing_tag(fasta):
    """The taxid already in this FASTA's first header, or None if untagged."""
    with opener(fasta) as fh:
        first = fh.readline()
    if not first.startswith(">" + TAG_PREFIX):
        return None
    tid = first[len(">" + TAG_PREFIX):].split("|", 1)[0]
    return int(tid) if tid.isdigit() else None


def tag_fasta(src, taxid, dest):
    """Copy src to dest with every header tagged `>kraken:taxid|TAXID|...`.

    Tagging lets Kraken2 assign sequence to a taxon directly and skips the 15 GB
    accession2taxid lookup. Any tag already present is REPLACED, never preserved. Host CDS
    arrives gzipped and is decompressed here, because `--add-to-library` wants plain FASTA.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    pre = ">" + TAG_PREFIX
    with opener(src) as fin, open(dest, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                if header.startswith(TAG_PREFIX):
                    parts = header.split("|", 2)
                    header = parts[2] if len(parts) == 3 else header
                fout.write(f"{pre}{taxid}|{header}\n")
            else:
                fout.write(line)


def add_to_library(fna, db_dir):
    """kraken2-build --add-to-library. Returns "" on success or the error text."""
    r = subprocess.run(["kraken2-build", "--add-to-library", str(fna),
                        "--db", str(db_dir), "--no-masking"],
                       capture_output=True, text=True)
    return "" if r.returncode == 0 else (r.stderr or r.stdout or "unknown").strip()[:200]


def install_taxonomy(db_dir):
    """nodes.dmp + names.dmp out of the module's pinned taxdump archive."""
    tax = db_dir / "taxonomy"
    if (tax / "nodes.dmp").exists() and (tax / "names.dmp").exists():
        print("taxonomy already installed")
        return
    if not TAXDMP.exists():
        sys.exit(f"no taxdump at {TAXDMP} — run fetch_species.py first")
    tax.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(TAXDMP) as z:
        for name in ("nodes.dmp", "names.dmp"):
            if name not in z.namelist():
                sys.exit(f"{name} is not in {TAXDMP}")
            with z.open(name) as src, open(tax / name, "wb") as out:
                shutil.copyfileobj(src, out)
    print(f"taxonomy installed from {TAXDMP.name} (the pinned archive, not a fresh download)")


def build(db_dir, threads):
    print(f"\nkraken2-build --build  (k={KMER_LEN}, l={MINIMIZER_LEN}, threads={threads})",
          flush=True)
    r = subprocess.run(["kraken2-build", "--build", "--db", str(db_dir),
                        "--kmer-len", str(KMER_LEN), "--minimizer-len", str(MINIMIZER_LEN),
                        "--threads", str(threads)], text=True)
    if r.returncode != 0:
        sys.exit(f"kraken2-build --build failed (exit {r.returncode})")


def inspect(db_dir, out, threads):
    """kraken2-inspect -> the minimizers-per-taxon table detection divides by.

    Not optional bookkeeping: coverage is measured against the distinct minimizers a taxon
    holds in THIS database, so a stale inspect table computes every coverage number against a
    previous database's denominators. Tracked for that reason, and because it has to exist on
    both machines.
    """
    print(f"\nkraken2-inspect -> {out}", flush=True)
    r = subprocess.run(["kraken2-inspect", "--db", str(db_dir), "--threads", str(threads)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"kraken2-inspect failed: {(r.stderr or '').strip()[:300]}")
    with atomic_write(out) as fh:
        fh.write(r.stdout)
    total, kingdoms = 0, []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) > 5 and f[3] == "" and f[5].strip() == "root":
            total = int(f[1])
        if len(f) > 5 and f[3] == "K":
            kingdoms.append((f[5].strip(), f[0].strip()))
    print(f"  root minimizers: {total:,}")
    for name, pct in kingdoms:
        print(f"    {name:20s} {pct:>7s}%")
    if len(kingdoms) == 1:
        print("  ONE KINGDOM ONLY. If that is Fungi, the host sink is missing again: see "
              "the host tier in build_plan.tsv and fetch_host_cds.py.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--name", required=True,
                    help="database name, e.g. db1p; the index goes to data/db/<name>/")
    ap.add_argument("--tier", default="db0_baseline,db1_congener",
                    help="comma-separated ladder tiers; the host tier is always included")
    ap.add_argument("--db-dir", type=Path, default=None,
                    help="override the index location (default data/db/<name>/)")
    ap.add_argument("--extra-cds", type=Path,
                    help="TSV of accession, taxid, path to add beyond the plan (db2's "
                         "predicted CDS)")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve and tag-check every input, add nothing, build nothing")
    ap.add_argument("--build-only", action="store_true",
                    help="skip the library pass: taxonomy, --build and inspect on an existing "
                         "library, for resuming a build that died in --build")
    args = ap.parse_args()

    tiers = {t.strip() for t in args.tier.split(",") if t.strip()}
    db_dir = args.db_dir or DB_ROOT / args.name
    db_dir.mkdir(parents=True, exist_ok=True)
    stage = db_dir / ".staging"
    inspect_out = OUTPUTS / f"db_inspect_{args.name}.tsv"

    if not args.build_only:
        rows, missing = plan_rows(tiers, read_extra(args.extra_cds))
        if missing:
            path = db_dir / "build_failures.tsv"
            with atomic_write(path) as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["accession", "db_tier", "error"])
                w.writerows(missing)
            for acc, tier, err in missing[:8]:
                print(f"  MISSING {acc} ({tier}): {err}")
            sys.exit(f"\n{len(missing):,} inputs have no CDS on disk -> {path}\n"
                     "NOT building: an index that is quietly short makes every coverage "
                     "denominator wrong in an unrecorded direction. Run the fetch step for "
                     "the tier and resubmit.")

        by_tier = Counter()
        for r in rd(PLAN):
            if r["db_tier"] == HOST_TIER or r["db_tier"] in tiers:
                by_tier[r["db_tier"]] += 1
        print(f"{args.name}: {len(rows):,} CDS files, "
              f"{len({t for _, t, _, _ in rows}):,} distinct taxids")
        for tier in sorted(by_tier):
            print(f"  {tier:16s} {by_tier[tier]:>6,}")
        print(f"  index -> {db_dir}")

        tags = Counter()
        failures, added = [], 0
        for i, (acc, taxid, src, name) in enumerate(rows, 1):
            old = existing_tag(src)
            tags["untagged" if old is None
                 else "tag matches" if old == taxid else "RE-TAGGED"] += 1
            if old is not None and old != taxid:
                print(f"  re-tag {acc}: {old} -> {taxid}  ({name})")
            if args.dry_run:
                continue
            dest = stage / acc / "tagged.fna"
            tag_fasta(src, taxid, dest)
            err = add_to_library(dest, db_dir)
            shutil.rmtree(dest.parent, ignore_errors=True)
            if err:
                failures.append((acc, "", err))
            else:
                added += 1
            if i % 250 == 0 or i == len(rows):
                print(f"  [{i:,}/{len(rows):,}]  added {added:,}", flush=True)

        print("\n-- header tags --")
        for k in ("untagged", "tag matches", "RE-TAGGED"):
            if tags.get(k):
                print(f"  {k:12s} {tags[k]:6,d}")
        if tags.get("RE-TAGGED"):
            print("  RE-TAGGED rows carried a taxid from a previous database. The plan wins; "
                  "this is the bug the 03_kraken tagger hid by checking only for a prefix.")
        shutil.rmtree(stage, ignore_errors=True)

        if failures:
            path = db_dir / "build_failures.tsv"
            with atomic_write(path) as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["accession", "db_tier", "error"])
                w.writerows(failures)
            sys.exit(f"\n{len(failures):,} of {len(rows):,} did not reach the library -> "
                     f"{path}\nNOT building. db_v3 warned here and built anyway, which is why "
                     "nothing downstream knew it was 34 assemblies short.")
        if args.dry_run:
            print("\n--dry-run: nothing added, nothing built.")
            return
        print(f"\nlibrary complete: {added:,} of {len(rows):,}")

    install_taxonomy(db_dir)
    build(db_dir, args.threads)
    inspect(db_dir, inspect_out, args.threads)

    hash_k2d = db_dir / "hash.k2d"
    if hash_k2d.exists():
        print(f"\ndb ready: {db_dir}")
        print(f"  hash.k2d {hash_k2d.stat().st_size / 1e9:.1f} GB "
              "(this is the RAM a classify job must hold)")
    print(f"  -> {inspect_out}")


if __name__ == "__main__":
    main()
