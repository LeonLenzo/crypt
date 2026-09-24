#!/usr/bin/env python3
"""
kraken_db_hosts.py — select host plant CDS for the Kraken2 database.

db_v2 excluded hosts. kraken_db_build.py still says so ("Host sequences are
intentionally excluded"). That was wrong for the co-infection question: a host read
with no correct home in the database is not left unclassified, it is absorbed by the
nearest fungal CDS and reported as a fungal call. With 1,349 wheat runs and no
Triticum in db_v2, that is a live source of false positives and a likely contributor
to the Phakopsora-in-every-wheat-sample artefact. See the memory notes.

The host is the most abundant organism in a plant RNA-seq library, so it is the single
most valuable LCA competitor available. This script gives it a home.

Design, and how it differs from the pathogen sweep in kraken_db_search.py:

  Curate first.   The resolved host list in run_list.tsv is submitter metadata and
                  carries errors: Bos grunniens (yak) is not a plant host, and several
                  entries are genus-level taxids with no species assembly. --audit
                  reports these; they are dropped or flagged, never silently fetched.

  A pangenome floor. The dominant hosts have real pan-transcriptomes on Ensembl
                  (wheat ~19 cultivars, barley and rice similar), and each cultivar
                  adds accessory transcripts that catch host reads a single reference
                  misses. --big-floor cultivars for the dominant hosts (wheat, rice,
                  barley), 1 reference for the rest. Kraken2 stores each minimizer
                  once, so shared transcripts across cultivars add no weight; only the
                  cultivar-specific sequence costs anything, and there is no
                  discrimination downside since host identity comes from metadata.

  Two sources.    Ensembl Plants first, NCBI second. Ensembl carries curated plant
                  annotation and the cultivar pan-transcriptomes NCBI lacks; it is also
                  where NCBI plant annotation is patchiest. CDS is fetched from the
                  EnsemblGenomes FTP directly (its release numbering and species-dir
                  names differ from the ensembl.org REST list, and its annotation does
                  not appear in NCBI's record for the same GCA — routing an Ensembl GCA
                  through `datasets --include cds` returns nothing). NCBI `datasets`
                  covers the long tail Ensembl does not carry.

Output mirrors the pathogen table so the download and build steps treat hosts
identically: data/host_candidates.tsv with the same columns, source=host_ensembl or
host_ncbi, then CDS into data/cds/host/{accession}/.

Usage:
    python kraken/db/kraken_db_hosts.py --audit          # report the host list, fetch nothing
    python kraken/db/kraken_db_hosts.py                  # select, write host_candidates.tsv
    python kraken/db/kraken_db_hosts.py --download       # select + fetch CDS
"""

import argparse
import csv
import gzip
import re
import shutil
import sys
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Reuse the pathogen script's download + taxonomy helpers rather than fork them.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kraken_db_search as S  # noqa: E402

OUT_DIR   = Path("kraken/output/db/search")
DATA_DIR  = OUT_DIR / "data"
RUN_LIST  = Path("kraken/output/run/select/data/run_list.tsv")
HOST_TSV  = DATA_DIR / "host_candidates.tsv"
HOST_CDS  = DATA_DIR / "cds" / "host"
ENSEMBL_FTP  = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants"

# Hosts that get more than one assembly: the few that dominate the cohort. Wheat, rice,
# barley are >50%, >5%, >3% of runs respectively; the rest have a long thin tail where
# one reference is plenty. Keyed by NCBI taxid.
BIG_HOSTS = {4565, 4530, 4513}   # Triticum aestivum, Oryza sativa, Hordeum vulgare

# Resolved host taxids that are not plant hosts — metadata errors in run_list.tsv that
# resolved faithfully from wrong submitter input. Their reads stay in the cohort but
# get no host CDS. Bos grunniens (yak): all 9 runs are PRJNA994854, mislabelled at
# submission. Add others here as --audit surfaces them.
NOT_A_HOST = {30521}


def load_host_runs(run_list: Path) -> dict:
    """{ncbi taxid: run count} from run_list.tsv, unresolved rows dropped."""
    counts = defaultdict(int)
    with open(run_list) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = (r.get("host_taxid") or "").strip()
            if t.isdigit():
                counts[int(t)] += 1
    return dict(counts)


def _ftp_listing(url: str) -> list:
    """Return the directory/file names in an Ensembl FTP HTML index."""
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=120) as resp:
        html = resp.read().decode("utf-8", "replace")
    # Apache autoindex: href="name/" for dirs, href="name" for files.
    return re.findall(r'href="([^"?/][^"]*)"', html)


def ensembl_current_release() -> int:
    """Highest release-N under the EnsemblGenomes plants FTP.

    EnsemblGenomes numbers releases separately from ensembl.org (the REST species
    list's `release` field is the ensembl.org number and does NOT match the FTP path),
    so the release must be read from the FTP itself.
    """
    rels = [int(m.group(1)) for name in _ftp_listing(f"{ENSEMBL_FTP}/")
            if (m := re.match(r"release-(\d+)/?$", name))]
    if not rels:
        raise RuntimeError(f"no release-N dirs under {ENSEMBL_FTP}")
    return max(rels)


def ensembl_species_index(release: int) -> list:
    """Sorted list of species directory names under a release's fasta/ tree.

    These are production names: the reference is `genus_species`, and pangenome
    cultivars are separate dirs like `triticum_aestivum_jagger`. Wheat alone has ~19.
    """
    base = f"{ENSEMBL_FTP}/release-{release}/fasta/"
    return sorted(n.rstrip("/") for n in _ftp_listing(base) if n.endswith("/"))


def ensembl_dirs_for(latin_name: str, species_index: list, floor: int) -> list:
    """Up to `floor` Ensembl species dirs for a host, reference first then cultivars.

    Matches `genus_species` exactly (the reference assembly) and as a prefix of
    `genus_species_<cultivar>` (pangenome members), so the dominant hosts pull their
    whole pan-transcriptome up to the floor while the tail takes one reference.
    """
    prefix = "_".join(latin_name.lower().split()[:2])
    exact = [d for d in species_index if d == prefix]
    cultivars = sorted(d for d in species_index if d.startswith(prefix + "_"))
    return (exact + cultivars)[:floor]


def download_ensembl_cds(species_dir: str, release: int, dest: Path) -> list:
    """Download the CDS FASTA for one Ensembl species dir. Resumable: skips if a
    .fna is already present. Returns [.fna path] or [] on failure."""
    dest.mkdir(parents=True, exist_ok=True)
    existing = [f for f in dest.glob("*.fna")]
    if existing:
        return existing
    cds_url = f"{ENSEMBL_FTP}/release-{release}/fasta/{species_dir}/cds/"
    try:
        names = _ftp_listing(cds_url)
    except Exception as e:
        print(f"  [{species_dir}] listing FAILED ({e})", flush=True)
        return []
    gz = next((n for n in names if n.endswith(".cds.all.fa.gz")), None)
    if not gz:
        print(f"  [{species_dir}] no cds.all.fa.gz", flush=True)
        return []
    out = dest / (species_dir + ".cds.fna")
    try:
        with urllib.request.urlopen(cds_url + gz, timeout=1800) as resp, \
                gzip.GzipFile(fileobj=resp) as fin, open(out, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    except Exception as e:
        print(f"  [{species_dir}] download FAILED ({e})", flush=True)
        out.unlink(missing_ok=True)
        return []
    return [out]


def ncbi_best_annotated(taxid: int, floor: int) -> list:
    """Up to `floor` annotated assemblies for taxid from NCBI, best first.

    Best = assembly level then protein-coding gene count, the same quality order the
    pathogen selector uses. Returns [(accession, gene_count), ...].
    """
    level = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}
    scored = []
    for a in S.datasets_query(taxid):
        ann = a.get("annotation_info")
        if not ann:
            continue
        genes = (ann.get("stats", {}).get("gene_counts", {}) or {}).get("protein_coding") or 0
        lvl = level.get(a.get("assembly_info", {}).get("assembly_level", ""), 0)
        scored.append(((lvl, genes), a.get("accession"), genes))
    scored.sort(reverse=True)
    return [(acc, genes) for _, acc, genes in scored[:floor]]


def select_hosts(run_list: Path, big_floor: int) -> tuple:
    """Return (candidate rows, audit dict). One row per chosen assembly, columns
    matching CANDIDATE_COLS so the download and build steps are source-agnostic."""
    to_species, sci_name = S.holobase_species_map()
    runs = load_host_runs(run_list)
    release = ensembl_current_release()
    species_index = ensembl_species_index(release)
    print(f"Ensembl Plants release-{release}: {len(species_index)} species dirs", flush=True)

    rows, audit = [], {"not_host": [], "no_assembly": []}
    for taxid, n_runs in sorted(runs.items(), key=lambda kv: -kv[1]):
        name = sci_name.get(taxid, str(taxid))
        if taxid in NOT_A_HOST:
            audit["not_host"].append((taxid, name, n_runs))
            continue

        floor = big_floor if taxid in BIG_HOSTS else 1

        # Ensembl first — it carries curated plant annotation and the pangenome
        # cultivar transcriptomes NCBI lacks (wheat alone has ~19). NCBI covers the
        # tail Ensembl does not annotate. accession is the Ensembl species dir or the
        # GCA; the download step routes on `source`.
        picks, source = [], None
        ens_dirs = ensembl_dirs_for(name, species_index, floor)
        if ens_dirs:
            source = "host_ensembl"
            picks = [(d, "ensembl:" + str(release), 0) for d in ens_dirs]
        else:
            got = ncbi_best_annotated(taxid, floor)
            if got:
                source = "host_ncbi"
                picks = [(acc, "", genes) for acc, genes in got]

        if not picks:
            audit["no_assembly"].append((taxid, name, n_runs))
            continue

        for rank, (acc, asm, genes) in enumerate(picks, 1):
            if not acc:
                continue
            rows.append({
                "taxid": taxid,
                "organism_name": name,
                "kingdom": "host",
                "source": source,
                "accession": acc,
                "assembly_level": asm or "",
                "release_date": "",
                "has_annotation": "True",
                "country": "",
                "protein_coding_genes": str(genes) if genes else "",
                "scaffold_n50_kb": "",
                "total_length_mb": "",
                "fasta_type": "cds",
                "busco_lineage": "",   # hosts are not BUSCO-screened
                "selection_rank": rank,
                "selection_reason": f"host_floor:{floor}",
            })
    return rows, audit


def print_audit(runs_total, rows, audit):
    taxa = {r["taxid"] for r in rows}
    print(f"\n── Host selection ───────────────────────────────────────────────")
    print(f"  Host taxa fetchable : {len(taxa)}")
    print(f"  Assemblies selected : {len(rows)}")
    by_source = defaultdict(int)
    for r in rows:
        by_source[r["source"]] += 1
    for s, n in sorted(by_source.items()):
        print(f"    {s:<16} {n:>4}")
    for label, key in (("not a plant host", "not_host"),
                       ("no fetchable assembly", "no_assembly")):
        items = audit[key]
        if items:
            runs = sum(n for *_, n in items)
            print(f"  {label}: {len(items)} taxa, {runs} runs")
            for t, nm, n in sorted(items, key=lambda x: -x[2])[:8]:
                print(f"      {t:<8}{nm[:34]:<36}{n} runs")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-list", default=str(RUN_LIST))
    ap.add_argument("--big-floor", type=int, default=5,
                    help="Ensembl cultivar transcriptomes for the dominant hosts "
                         "(wheat/rice/barley); 1 reference for every other host")
    ap.add_argument("--audit", action="store_true",
                    help="report the host list and exit, fetch nothing")
    ap.add_argument("--download", action="store_true",
                    help="download CDS for the selected hosts")
    ap.add_argument("--reselect", action="store_true",
                    help="re-run host selection even if host_candidates.tsv exists. "
                         "Off by default so a Setonix run downloads the committed, "
                         "reviewed selection rather than re-querying Ensembl/NCBI, "
                         "whose catalogues move — the same reason kraken_db_search.py "
                         "uses --from-table.")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    # Download the committed table by default; only re-select when asked or when no
    # table exists yet. Selection is a reviewed local step, like the pathogen side.
    if HOST_TSV.exists() and not args.reselect and not args.audit:
        with open(HOST_TSV) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        print(f"Using committed selection: {HOST_TSV} ({len(rows)} rows). "
              f"--reselect to regenerate.")
    else:
        rows, audit = select_hosts(Path(args.run_list), args.big_floor)
        runs_total = sum(load_host_runs(Path(args.run_list)).values())
        print_audit(runs_total, rows, audit)
        if args.audit:
            return
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        with open(HOST_TSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=S.CANDIDATE_COLS, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {HOST_TSV} ({len(rows)} rows)")

    if args.download:
        print(f"\nDownloading host CDS → {HOST_CDS}")
        HOST_CDS.mkdir(parents=True, exist_ok=True)
        n_ok = n_fail = 0

        def work(row):
            acc = row["accession"]
            dest = HOST_CDS / acc
            if row["source"] == "host_ensembl":
                # assembly_level holds "ensembl:<release>" for these rows.
                release = int(row["assembly_level"].split(":", 1)[1])
                return acc, download_ensembl_cds(acc, release, dest)
            return acc, S.download_cds(acc, dest)

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(work, r): r["accession"] for r in rows}
            for done, fut in enumerate(as_completed(futs), 1):
                acc, fnas = fut.result()
                if fnas:
                    n_ok += 1
                else:
                    n_fail += 1
                if done % 5 == 0 or done == len(rows):
                    print(f"  [{done}/{len(rows)}] ok={n_ok} fail={n_fail}",
                          end="\r", flush=True)
        print(f"\n  Host CDS: {n_ok} ok, {n_fail} failed")
    else:
        print("Next: python kraken/db/kraken_db_hosts.py --download")


if __name__ == "__main__":
    main()
