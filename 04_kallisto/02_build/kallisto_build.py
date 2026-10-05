#!/usr/bin/env python3
"""kallisto_build.py — one kallisto index per host, from tagged CDS.

Each index holds three groups of reference, for three different reasons:

  host       so host reads are absorbed rather than forced onto a pathogen
  any-level  every species Kraken2 saw in that host AT ALL, not only the flagged ones. A
             species that is genuinely present but missing from the index has its reads
             pushed onto whatever relative IS present, inflating exactly what this module
             is meant to test. One assembly each makes this cheap insurance.
  neighbour  the measured >=1% k-mer containment partners of the flagged species, from
             02_build/data/neighbours.tsv. Rank-agnostic on purpose: relatedness does not
             predict shared 31-mers, containment does. See kallisto_neighbours.py.

Scope comes from 01_select: one index per host in the stratified set, because that is what
03_quant will need. The selection guarantees every flagged host is represented, so this is the
same 19 hosts as `flagged_hosts()` and no host can be silently dropped.

Run on Setonix (kallisto via spack; CDS in 03_kraken/01_search/data/cds/):
    sbatch 04_kallisto/slurm/kallisto_build_index.slurm
Or directly, from the repo root:
    python 04_kallisto/01_build/kallisto_build.py --host "Triticum aestivum" --threads 32
    python 04_kallisto/01_build/kallisto_build.py --all --dry-run
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (CDS_DIR, HOST_CDS, HOST_V3, IDX_DIR, NEIGHBOURS, SCORE_MIN, COV_MIN,
                     STRATA, busco_ranked, classified, flagged_hosts, num, slug,
                     step_log)


def hosts_to_build():
    """Hosts needing an index: those in the stratified set, falling back to every flagged host
    if 01_select has not run yet."""
    if STRATA.exists():
        return sorted({r["host"] for r in csv.DictReader(open(STRATA), delimiter="\t")})
    return flagged_hosts()

# Species whose one-assembly k-mer coverage is too poor for a single representative. Both are
# known species complexes: saturation_summary puts R. solani at pangenome_ratio 14.3 (one
# assembly = 7% of species k-mer space) and F. oxysporum at 5.25 (19%). MULTI_N is a
# compromise, not a fix: their pangenomes are open, and n_for_95pct is 48 and 56 respectively.
MULTI_ASSEMBLY = {"Rhizoctonia solani", "Fusarium oxysporum"}
MULTI_N = 8
NEIGHBOUR_MIN = 0.01   # containment floor for putting a cross-talk partner in the index


def cds_path(accession):
    """CDS fasta for a pathogen accession, as laid down by the db-build download."""
    return CDS_DIR / accession / "ncbi_dataset/data" / accession / "cds_from_genomic.fna"


def host_cds(host):
    """CDS fasta for a host.

    Two shapes, because the db build drew hosts from two places: NCBI-sourced hosts nest like
    the pathogens (<accession>/ncbi_dataset/data/<accession>/cds_from_genomic.fna), Ensembl
    ones sit as <species_slug>/<species_slug>.cds.fna. Glob rather than assume either.
    """
    for r in csv.DictReader(open(HOST_V3), delimiter="\t"):
        if r["organism_name"] == host and r.get("accession"):
            d = HOST_CDS / r["accession"]
            if d.is_dir():
                for f in sorted(d.rglob("*.fna")):
                    return f
    s = slug(host)
    if HOST_CDS.is_dir():
        for d in sorted(HOST_CDS.iterdir()):
            if d.name == s or d.name.startswith(s):
                for f in sorted(d.rglob("*.fna")):
                    return f
    return None


def reference_set(host):
    """Returns (species -> [CDS fasta], flagged species, species with no CDS on disk)."""
    rows = classified(host)
    anylevel = {r["species"] for r in rows}
    flag = {r["species"] for r in rows
            if num(r, "score") >= SCORE_MIN and num(r, "coverage") >= COV_MIN}

    neigh = set()
    if NEIGHBOURS.exists():
        for r in csv.DictReader(open(NEIGHBOURS), delimiter="\t"):
            if r["species"] in flag and num(r, "containment_a_in_b") >= NEIGHBOUR_MIN:
                neigh.add(r["neighbour"])

    ranked = busco_ranked()
    out, missing = {}, []
    for sp in sorted(anylevel | neigh):
        accs = ranked.get(sp)
        if not accs:
            missing.append(sp)
            continue
        paths = [p for p in (cds_path(a) for a in
                             accs[:MULTI_N if sp in MULTI_ASSEMBLY else 1]) if p.exists()]
        (out.setdefault(sp, paths) if paths else missing.append(sp))
    return out, sorted(flag), missing


def retag(fasta, species, accession, fh):
    """Stream a CDS fasta into fh, rewriting headers to {species_slug}|{accession}|{id}.

    kallisto reports per-target abundance, so the species has to live in the target name for
    abundances to aggregate up, and the original transcript id has to survive after it because
    gene breadth counts distinct genes. Incoming headers carry the db build's `kraken:taxid|N|`
    prefix, which is stripped so the tag is unambiguous.
    """
    tag, n = slug(species), 0
    with open(fasta) as src:
        for line in src:
            if line.startswith(">"):
                rest = line[1:].strip()
                if rest.startswith("kraken:taxid|"):
                    parts = rest.split("|", 2)
                    rest = parts[2] if len(parts) > 2 else parts[-1]
                fh.write(f">{tag}|{accession}|{rest.split()[0]}\n")
                n += 1
            else:
                fh.write(line)
    return n


def index_is_valid(idx):
    """True if idx exists and kallisto can read it.

    Resume must test integrity, not existence: the 2026-09-14 truncated-gzip incident came
    from a resume check that only asked whether a file was there. `kallisto inspect` parses the
    header, so a half-written index fails it.
    """
    if not idx.exists() or idx.stat().st_size == 0:
        return False
    r = subprocess.run(["kallisto", "inspect", str(idx)],
                       capture_output=True, text=True)
    return r.returncode == 0


def build_index(host, threads=8, dry_run=False, force=False):
    idx_existing = IDX_DIR / f"{slug(host)}.idx"
    if not force and not dry_run and index_is_valid(idx_existing):
        print(f"{host}: index already built and valid "
              f"({idx_existing.stat().st_size / 1e9:.2f} GB), skipping")
        return idx_existing
    refs, flag, missing = reference_set(host)
    hcds = host_cds(host)
    if hcds is None:
        print(f"  !! no host CDS for {host}; its reads cannot be absorbed", file=sys.stderr)

    IDX_DIR.mkdir(parents=True, exist_ok=True)
    fa  = IDX_DIR / f"{slug(host)}.cds.fna"
    idx = IDX_DIR / f"{slug(host)}.idx"

    nbytes = sum(p.stat().st_size for ps in refs.values() for p in ps) + \
             (hcds.stat().st_size if hcds else 0)
    print(f"{host}: {len(refs) + (1 if hcds else 0)} species ({len(flag)} flagged), "
          f"{sum(len(v) for v in refs.values())} assemblies, {nbytes / 1e9:.2f} GB of CDS")
    if missing:
        print(f"  {len(missing)} species skipped, no CDS on disk: {', '.join(missing[:6])}"
              f"{' ...' if len(missing) > 6 else ''}")
    if dry_run:
        return idx

    ntx = 0
    with open(fa, "w") as fh:
        if hcds:
            ntx += retag(hcds, host, "host", fh)
        for sp, paths in sorted(refs.items()):
            for p in paths:
                # .../cds/pathogen/<acc>/ncbi_dataset/data/<acc>/cds_from_genomic.fna, so the
                # accession is the PARENT directory. It matters for MULTI_ASSEMBLY species,
                # where several assemblies are concatenated and target ids must stay unique.
                ntx += retag(p, sp, p.parent.name, fh)
    print(f"  wrote {fa.name}: {ntx} transcripts, {fa.stat().st_size / 1e9:.2f} GB")

    cmd = ["kallisto", "index", "-i", str(idx), "-t", str(threads), str(fa)]
    print(f"  $ {' '.join(cmd)}", flush=True)
    if subprocess.run(cmd).returncode != 0:
        print(f"  !! kallisto index failed for {host}", file=sys.stderr)
        return None
    print(f"  built {idx.name}: {idx.stat().st_size / 1e9:.2f} GB")
    return idx


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", help="one host species")
    ap.add_argument("--all", action="store_true",
                    help="every host in the stratified set (see 01_select)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true",
                    help="report index composition without writing or indexing")
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if a valid index is already present")
    args = ap.parse_args()

    hosts = [args.host] if args.host else hosts_to_build() if args.all else None
    if not hosts:
        ap.error("pass --host or --all")
    log = step_log(__file__, "kallisto_build")
    try:
        bad = sum(1 for h in hosts
                  if build_index(h, args.threads, args.dry_run, args.force) is None
                  and not args.dry_run)
    finally:
        log.close()
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
