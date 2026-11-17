#!/usr/bin/env python3
"""Download the host CDS from Ensembl Plants, one FASTA per selected genome.

    python release/02_reference/scripts/fetch_host_cds.py
    python release/02_reference/scripts/fetch_host_cds.py --dry-run
    python release/02_reference/scripts/fetch_host_cds.py --have acacia.txt

Reads  data/inputs/scope_rule.tsv              source=ensembl_plants, the pinned release
       data/curation/screen_hosts.tsv          which hosts and which Ensembl genomes
       --have FILE (optional)                  filenames already available, one per line
Writes data/raw/host_cds/<ensembl_species>/<file>.cds.all.fa.gz
       data/raw/host_cds_manifest.tsv           species -> file, bytes, md5, TRACKED

CDS ONLY, never host genomes. These are RNA-seq libraries, so the transcriptome is what reads
come from: wheat CDS is ~0.2 GB against ~15 GB genomic, and the full host genome set would be
~275 GB for no gain.

## THE FILENAME HAS TO BE DISCOVERED, NOT CONSTRUCTED

`species_EnsemblPlants.txt` carries an `assembly` column that does NOT match the assembly
token in the FASTA name. Landmark is `PGSBv2.0_Landmark` in the table and `PGSBv2.1` in the
file; Weebil is `EI_WeebilV1` against `WeebilV1`. Constructing
`<Species>.<assembly>.cds.all.fa.gz` 404s on 2 of the 19 wheat cultivars, and silently,
because a 404 on one cultivar of nineteen looks like nothing. So each `cds/` directory is
listed and the `*.cds.all.fa.gz` entry taken from it, and the resolved name is recorded in the
manifest so a reviewer can see exactly which file was used.

`.cds.all.fa.gz` and not `.cds.abinitio.fa.gz`: `all` is the gene-build CDS, `abinitio` is
unreviewed ab-initio predictions Ensembl ships alongside and does not include in its gene set.

## Rate limiting

One request at a time with a pause between. EBI shares a 2 req/s limiter across EuropePMC and
ENA (CLAUDE.md), and this project has already collapsed one EBI job by running 48 parallel
streams against it (1,012 connection refusals). 21 files is not worth risking that for.

Idempotent: a file already present with a non-zero size is skipped, so a killed run resumes.
"""

import argparse
import csv
import hashlib
import re
import shutil
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parent.parent
RULE = HERE / "data/inputs/scope_rule.tsv"
HOSTS = HERE / "data/curation/screen_hosts.tsv"
DEST = HERE / "data/raw/host_cds"
MANIFEST = HERE / "data/raw/host_cds_manifest.tsv"
UA = {"User-Agent": "crypt/02_reference (leon.lenzo@curtin.edu.au)"}
PAUSE = 0.6
FASTA_RE = re.compile(r"[A-Za-z0-9_.\-]+\.cds\.all\.fa\.gz")


def read_rule():
    out = {}
    with open(RULE) as fh:
        for row in csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                  delimiter="\t"):
            out.setdefault(row["source"], {})[row["field"]] = row["value"]
    return out


def rd(p):
    with open(p) as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                   delimiter="\t"))


def get(url, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


def resolve(base, species):
    """The real CDS filename for a species, by listing its cds/ directory."""
    listing = get(base.replace("{species}", species))
    names = sorted(set(FASTA_RE.findall(listing)))
    if not names:
        return None
    # One per species in practice; take the longest name if ever more, which is the most
    # specific rather than an arbitrary pick.
    return max(names, key=len)


def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=900) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh, length=1 << 20)
    tmp.replace(dest)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve filenames and sizes, download nothing")
    ap.add_argument("--have", type=Path,
                    help="file of CDS filenames already available; those are skipped")
    args = ap.parse_args()

    ens = read_rule()["ensembl_plants"]
    base = ens["cds_url"]
    have = set()
    if args.have and args.have.exists():
        have = {ln.strip() for ln in args.have.read_text().splitlines() if ln.strip()}

    jobs = []
    for h in rd(HOSTS):
        if h["decision"] != "keep" or not h["ensembl_species"]:
            continue
        for sp in h["ensembl_species"].split(";"):
            if sp:
                jobs.append((h["host"], h["role"], h["ncbi_taxid"], sp))

    print(f"Ensembl Plants {ens['release']}: {len(jobs)} host CDS files to resolve\n")
    rows, n_new, n_skip, n_fail, total = [], 0, 0, 0, 0
    for host, role, tid, sp in jobs:
        try:
            fname = resolve(base, sp)
        except Exception as e:
            print(f"  FAILED listing {sp}: {type(e).__name__}: {e}")
            n_fail += 1
            rows.append({"ensembl_species": sp, "host": host, "role": role,
                         "ncbi_taxid": tid, "file": "", "bytes": "", "md5": "",
                         "status": f"listing failed: {e}"})
            continue
        time.sleep(PAUSE)
        if not fname:
            print(f"  NO CDS for {sp}")
            n_fail += 1
            rows.append({"ensembl_species": sp, "host": host, "role": role,
                         "ncbi_taxid": tid, "file": "", "bytes": "", "md5": "",
                         "status": "no .cds.all.fa.gz in the directory"})
            continue
        dest = DEST / sp / fname
        url = base.replace("{species}", sp) + fname
        if fname in have:
            print(f"  have   {fname}")
            n_skip += 1
            rows.append({"ensembl_species": sp, "host": host, "role": role,
                         "ncbi_taxid": tid, "file": fname, "bytes": "", "md5": "",
                         "status": "supplied via --have"})
            continue
        if dest.exists() and dest.stat().st_size:
            n_skip += 1
            b = dest.stat().st_size
            total += b
            print(f"  skip   {fname}  {b / 1e6:.0f} MB (present)")
            rows.append({"ensembl_species": sp, "host": host, "role": role,
                         "ncbi_taxid": tid, "file": fname, "bytes": b,
                         "md5": md5_of(dest), "status": "present"})
            continue
        if args.dry_run:
            print(f"  would fetch {fname}")
            rows.append({"ensembl_species": sp, "host": host, "role": role,
                         "ncbi_taxid": tid, "file": fname, "bytes": "", "md5": "",
                         "status": "dry-run"})
            continue
        try:
            download(url, dest)
            b = dest.stat().st_size
            total += b
            n_new += 1
            print(f"  got    {fname}  {b / 1e6:.0f} MB")
            rows.append({"ensembl_species": sp, "host": host, "role": role,
                         "ncbi_taxid": tid, "file": fname, "bytes": b,
                         "md5": md5_of(dest), "status": "fetched"})
        except Exception as e:
            print(f"  FAILED {fname}: {type(e).__name__}: {e}")
            n_fail += 1
            rows.append({"ensembl_species": sp, "host": host, "role": role,
                         "ncbi_taxid": tid, "file": fname, "bytes": "", "md5": "",
                         "status": f"failed: {e}"})
        time.sleep(PAUSE)

    cols = ["ensembl_species", "host", "role", "ncbi_taxid", "file", "bytes", "md5", "status"]
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\n  {n_new} fetched, {n_skip} already present, {n_fail} failed; "
          f"{total / 1e9:.2f} GB gz on disk")
    print(f"  -> {MANIFEST.relative_to(REPO)}")
    if n_fail:
        raise SystemExit(f"{n_fail} host CDS files are missing; the sink would be incomplete")
    print("Next: plan_build.py, then the Kraken2 build")


if __name__ == "__main__":
    main()
