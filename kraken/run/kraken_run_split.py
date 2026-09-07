#!/usr/bin/env python3
"""
kraken_run_split.py — host genome resolution + host-read removal via
per-taxid BBMap indices.
Submodule 2, step 2 of 3 (select -> split -> assign).

Owns host reference genome resolution/download as of the 2026-09-07 refactor
(moved out of kraken_run_select.py — see that script's docstring for why:
`select` should stay a cheap, fast, safe-to-rerun filtering step, not do
multi-GB network I/O as a side effect. `split` is the natural owner since it
needs the genome on disk to build its BBMap index anyway).

Design (agreed 2026-09-01, see kraken/README.md's "kraken_run_split.py design"
section for the full rationale): ONE bbmap.sh index PER distinct host taxid, not
one combined multi-reference BBSplit index. The 2,719-sample field/aerial cohort
names 116 candidate host taxids; 94 have an NCBI genomic assembly totalling
~275Gb, dominated by outliers (Pinus radiata 22.4Gb, Allium sativum 16.5Gb,
Triticum aestivum 14.5Gb, Avena insularis 14.2Gb, ...) that can't just be
excluded — they're the dominant crops in the cohort, not obscure edge cases.
A single combined index at that scale isn't buildable on any Setonix node;
indexing each genome separately is trivial even for the largest one alone.

Three stages, each independently resumable:

  0. RESOLVE + DOWNLOAD HOST GENOMES — for every distinct host taxid named
     across run_list.tsv's candidate_host_taxids column (kraken_run_select.py's
     output), pick and fetch its single best genomic assembly (not CDS —
     BBMap aligns reads rather than doing k-mer LCA, so it has no use for
     CDS/annotation, and NCBI's plant gene-annotation coverage is patchy
     enough that requiring it would exclude most hosts) via
     kraken_db_search.download_cds(..., include="genome") into
     kraken_db_search/data/cds/host/ — the same shared pool
     kraken_db_build.py's pathogen fetch uses (cds/pathogen/ there). Persists
     the resolved taxid->accession map so re-runs skip already-fetched taxids.

  1. BUILD INDEX — one bbmap.sh index per distinct host taxid resolved above.
     The genomic FASTA is already extracted (kraken_db_search.download_cds()
     unzips on download) — no separate "extract" stage needed. Skips any
     taxid whose index already exists. Embarrassingly parallel across taxids,
     but kept to a modest default worker count (--workers) since a handful of
     these genomes are tens of Gb and each bbmap.sh index build gets its own
     -Xmx allocation — too much concurrency risks overcommitting node memory.

  2. SPLIT — for each run in run_list.tsv, align its reads against EVERY named
     candidate's index separately (not one combined pass). The candidate with
     the highest mapped fraction is the confirmed host — an independent,
     read-level confirmation, often more reliable than metadata alone. Its
     unmapped reads are the pathogen-enriched output for kraken_run_assign.py.
     Other candidates' mapped fractions are kept as QC/confirmation metadata
     only — no cross-candidate intersection of unmapped reads. Mapped fraction
     is computed directly from read counts (input reads vs. unmapped-output
     reads), not parsed from BBMap's statsfile text, so it doesn't depend on
     guessing an exact key format.

**Two aligners available (--aligner {bbmap,hisat2}, default bbmap for now)**:
BBMap has a hard, undocumented-until-hit 500Mbp-per-chromosome limit — confirmed
2026-09-07 against real data: Triticum aestivum (wheat, ~55% of the whole cohort)
fails outright with "AssertionError ... reference file appears empty" (a misleading
message; the real cause, per BBMap's author on the project tracker, is simply that
wheat's chromosomes exceed the 500Mbp ceiling BBMap doesn't support, full stop, no
flag or workaround). HISAT2 uses an FM-index (BWT-based), not BBMap's in-memory
per-chromosome array approach, so it doesn't share this ceiling, and is the
standard splice-aware choice for RNA-seq anyway (BBMap is a generic aligner, not
RNA-seq-specific). Added 2026-09-07 to A/B test against BBMap and a raw-reads
(no split at all) Kraken2 baseline, per Leon's question: does host-read removal
even have a material effect on Kraken2's output, given the current Kraken2 DB
(db_v2) is pathogen-only (host genomes deliberately excluded, see
kraken_db_build.py) — so a host read can only produce a false pathogen hit if it
accidentally shares a k-mer with a real pathogen sequence, not from "competing"
for classification the way it might with a combined DB. Each aligner's indices
live under their own subdirectory so both can be tested side by side without
collision.

Run on Setonix (requires the bbmap and/or hisat2 module + NCBI datasets CLI):
    module load bbmap/38.96--h5c4e2a8_0        # for --aligner bbmap (default)
    module load hisat2/2.2.1-w7a5u7v            # for --aligner hisat2
    python kraken/run/kraken_run_split.py --build-index                    # stages 0+1 only
    python kraken/run/kraken_run_split.py --aligner hisat2                 # all three stages, HISAT2
    python kraken/run/kraken_run_split.py --limit 5                        # smoke test

Output:
    kraken/output/run/split/data/host_taxid_to_accession.json  (tracked — every
        candidate taxid's downloaded accession; owned here now, not by select)
    kraken/output/run/split/data/index/{aligner}/{taxid}/     (gitignored —
        one index dir per (aligner, host taxid) pair)
    kraken/output/run/split/data/reads/{run}_{1,2}.fastq.gz    (gitignored —
        confirmed-host-removed, pathogen-enriched reads; NOTE: shared output
        path regardless of --aligner used to produce it — don't run both
        aligners' split stage back to back without moving/renaming results
        in between, the second run will silently overwrite the first's output)
    kraken/output/run/split/data/split_results.tsv             (tracked — one
        row per run: confirmed host taxid + mapped %, every candidate's mapped
        % for QC, agreement with meta_classify.py's llm_host_resolved)
"""

import argparse
import csv
import gzip
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _util import _Tee, make_log_dir, link_latest, load_json, save_json

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "db"))
from kraken_db_search import download_cds, datasets_query, quality_key, scaffold_plus

# ── paths ──────────────────────────────────────────────────────────────────────

SELECT_DATA     = Path("kraken/output/run/select/data")
RUN_LIST        = SELECT_DATA / "run_list.tsv"
SELECT_READS_DIR = SELECT_DATA / "reads"
HOST_CDS_DIR    = Path("kraken/output/db/search/data/cds/host")

OUT_DIR        = Path("kraken/output/run/split")
DATA_DIR       = OUT_DIR / "data"
INDEX_DIR      = DATA_DIR / "index"
SPLIT_READS_DIR = DATA_DIR / "reads"
SPLIT_RESULTS  = DATA_DIR / "split_results.tsv"
# Every candidate taxid -> its downloaded accession. Owned by this script now
# (moved from kraken_run_select.py's data dir 2026-09-07) since this is the
# script that actually needs and fetches the genome.
HOST_TAXID_MAP = DATA_DIR / "host_taxid_to_accession.json"

SPLIT_RESULTS_COLS = [
    "Run", "BioSample", "BioProject", "llm_host_resolved",
    "confirmed_host_taxid", "confirmed_host_mapped_pct", "agrees_with_llm",
    "candidate_results", "status",
]

# BBMap needs its whole index (roughly proportional to reference size) resident
# in memory. 48g comfortably covers even the largest single host genome in this
# cohort (~22Gb Pinus radiata) with headroom — sized for ONE genome at a time,
# since indices are per-taxid now, not combined. BBMap CANNOT build for genomes
# with a chromosome >500Mbp regardless of memory (see module docstring) — this
# affects wheat and likely other cohort outliers; that's a hard limit, not a
# memory-tuning problem, and raising BBMAP_XMX will not fix it.
BBMAP_XMX = "48g"

_FASTA_EXCLUDE_SUFFIXES = ("_clean.fna", "_combined.fna", ".tagged.fna")

ALIGNERS = ("bbmap", "hisat2")


# ── stage 0: resolve + download host genomes (moved from kraken_run_select.py
# 2026-09-07 — no BUSCO screening needed here, unlike kraken_db_search.py's
# pathogen selection: host reference genomes are well-annotated model/crop
# species, not the strain-diversity problem BUSCO filtering solves) ──────────

def best_host_assembly(taxid: int) -> dict:
    """Pick the single best genomic assembly for a host taxid. Unlike
    kraken_db_search.py's pathogen selection, does NOT require annotation —
    BBMap aligns reads rather than doing k-mer LCA, so it has no use for
    CDS/annotation, and requiring it would exclude most plant hosts (NCBI's
    plant gene-annotation pipeline coverage is much patchier than fungi/
    vertebrates — even well-studied species like Nicotiana benthamiana often
    have zero NCBI-annotated assemblies despite having good genomic ones).
    Falls back from scaffold-plus down to any assembly — always take
    something over nothing."""
    assemblies = datasets_query(taxid)
    for pool in [
        [a for a in assemblies if scaffold_plus(a)],
        assemblies,
    ]:
        if pool:
            break
    if not pool:
        return None
    pool.sort(key=quality_key, reverse=True)
    return pool[0]


def ensure_host_cds(taxid: int) -> str:
    """Fetch the single best genomic assembly for a host taxid if not already
    present. Returns the accession, or '' on failure."""
    best = best_host_assembly(taxid)
    if not best:
        return ""
    accession = best.get("accession", "")
    if not accession:
        return ""
    dest = HOST_CDS_DIR / accession
    fnas = download_cds(accession, dest, include="genome")
    return accession if fnas else ""


# ── reference lookup ──────────────────────────────────────────────────────────

def find_host_fasta(accession: str) -> Path | None:
    """Locate the already-extracted genomic .fna for an accession.
    kraken_db_search.download_cds() (called by kraken_run_select.py) unzips on
    download, so this is just a lookup, never a fetch."""
    d = HOST_CDS_DIR / accession
    if not d.is_dir():
        return None
    fnas = [f for f in d.glob("**/*.fna") if not f.name.endswith(_FASTA_EXCLUDE_SUFFIXES)]
    return fnas[0] if fnas else None


# ── stage 1: build one index per taxid ────────────────────────────────────────

def _index_dir(aligner: str, taxid: str) -> Path:
    return INDEX_DIR / aligner / taxid


def _index_ready_bbmap(taxid: str) -> bool:
    """BBMap writes ref/genome/1/summary.txt on a successful index build."""
    return (_index_dir("bbmap", taxid) / "ref" / "genome" / "1" / "summary.txt").exists()


def _build_index_bbmap(taxid: str, accession: str) -> str:
    """Build a BBMap index for one host taxid. Resumable: skips if already
    built. Returns 'ok', 'cached', or 'failed'."""
    if _index_ready_bbmap(taxid):
        return "cached"
    fasta = find_host_fasta(accession)
    if not fasta:
        return "failed"
    idx_dir = _index_dir("bbmap", taxid)
    idx_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["bbmap.sh", f"-Xmx{BBMAP_XMX}",
           f"ref={fasta}", f"path={idx_dir}", "build=1", "overwrite=t"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0 or not _index_ready_bbmap(taxid):
        print(f"  taxid {taxid}: BBMap index build FAILED\n{r.stderr[-1500:]}", flush=True)
        return "failed"
    return "ok"


# HISAT2 uses an FM-index (Burrows-Wheeler transform), not BBMap's in-memory
# per-chromosome array approach — no known chromosome-length ceiling, and it's
# the standard splice-aware aligner for RNA-seq specifically (BBMap is generic).
# --large-index forced explicitly rather than relying on hisat2-build's own
# >4Gb auto-detection, since every genome in this cohort that matters is well
# past that threshold anyway — no downside to being explicit.
def _index_ready_hisat2(taxid: str) -> bool:
    """hisat2-build writes <base>.1.ht2 (or .1.ht2l for --large-index) on success."""
    base = _index_dir("hisat2", taxid) / taxid
    return base.with_suffix(".1.ht2l").exists() or base.with_suffix(".1.ht2").exists()


def _build_index_hisat2(taxid: str, accession: str, threads: int = 8) -> str:
    """Build a HISAT2 index for one host taxid. Resumable: skips if already
    built. Returns 'ok', 'cached', or 'failed'."""
    if _index_ready_hisat2(taxid):
        return "cached"
    fasta = find_host_fasta(accession)
    if not fasta:
        return "failed"
    idx_dir = _index_dir("hisat2", taxid)
    idx_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["hisat2-build", "--large-index", "-p", str(threads),
           str(fasta), str(idx_dir / taxid)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if r.returncode != 0 or not _index_ready_hisat2(taxid):
        print(f"  taxid {taxid}: HISAT2 index build FAILED\n{r.stderr[-1500:]}", flush=True)
        return "failed"
    return "ok"


def index_ready(taxid: str, aligner: str) -> bool:
    return _index_ready_hisat2(taxid) if aligner == "hisat2" else _index_ready_bbmap(taxid)


def build_index(taxid: str, accession: str, aligner: str, threads: int = 8) -> str:
    if aligner == "hisat2":
        return _build_index_hisat2(taxid, accession, threads)
    return _build_index_bbmap(taxid, accession)


# ── stage 2: split ─────────────────────────────────────────────────────────────

def _count_fastq_reads(path: Path) -> int:
    """Count reads in a (possibly gzipped) FASTQ file — 4 lines per read."""
    if not path.exists():
        return 0
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(path, "rt") as f:
        for _ in f:
            n += 1
    return n // 4


def _align_against_taxid_bbmap(taxid: str, run: str, r1: Path, r2: Path | None,
                               work_dir: Path) -> dict:
    idx_dir = _index_dir("bbmap", taxid)
    out_u1 = work_dir / f"{run}__bbmap__{taxid}__unmapped_1.fastq.gz"
    cmd = ["bbmap.sh", f"-Xmx{BBMAP_XMX}", f"path={idx_dir}", "build=1",
           "overwrite=t", "statsfile=stderr"]
    if r2 is not None:
        out_u2 = work_dir / f"{run}__bbmap__{taxid}__unmapped_2.fastq.gz"
        cmd += [f"in={r1}", f"in2={r2}", f"outu={out_u1}", f"outu2={out_u2}"]
    else:
        out_u2 = None
        cmd += [f"in={r1}", f"outu={out_u1}"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        return {}
    return {"unmapped_r1": out_u1, "unmapped_r2": out_u2}


def _align_against_taxid_hisat2(taxid: str, run: str, r1: Path, r2: Path | None,
                                work_dir: Path, threads: int = 8) -> dict:
    idx_dir = _index_dir("hisat2", taxid)
    prefix = work_dir / f"{run}__hisat2__{taxid}__unmapped"
    cmd = ["hisat2", "-p", str(threads), "-x", str(idx_dir / taxid),
           "-S", "/dev/null"]
    if r2 is not None:
        # hisat2's --un-conc-gz <prefix> writes concordantly-unmapped pairs to
        # <prefix>.1 / <prefix>.2 — content IS gzip-compressed (confirmed via a
        # live `gzip -c >` pipe on the actual running process, 2026-09-07), but
        # the FILENAME itself does NOT get a .gz suffix despite the flag name.
        # Rename to proper .gz-suffixed names immediately so downstream code
        # (_count_fastq_reads()'s suffix-based gzip detection, and split_run()'s
        # final promoted filenames) doesn't need special-casing per aligner.
        cmd += ["-1", str(r1), "-2", str(r2), "--un-conc-gz", str(prefix)]
        raw_u1, raw_u2 = Path(f"{prefix}.1"), Path(f"{prefix}.2")
        out_u1, out_u2 = Path(f"{prefix}.1.gz"), Path(f"{prefix}.2.gz")
    else:
        # Single-end: --un-gz <path> writes the exact filename given, still with
        # the same no-.gz-suffix behavior as above.
        raw_u1 = Path(f"{prefix}")
        out_u1 = Path(f"{prefix}.gz")
        raw_u2 = out_u2 = None
        cmd += ["-U", str(r1), "--un-gz", str(raw_u1)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        print(f"  {run}/{taxid}: hisat2 align FAILED\n{r.stderr[-1500:]}", flush=True)
        return {}
    if not raw_u1.exists():
        print(f"  {run}/{taxid}: hisat2 completed but expected output "
              f"{raw_u1} missing", flush=True)
        return {}
    raw_u1.rename(out_u1)
    if raw_u2 is not None:
        if not raw_u2.exists():
            print(f"  {run}/{taxid}: hisat2 completed but expected output "
                  f"{raw_u2} missing", flush=True)
            return {}
        raw_u2.rename(out_u2)
    return {"unmapped_r1": out_u1, "unmapped_r2": out_u2}


def align_against_taxid(taxid: str, run: str, r1: Path, r2: Path | None,
                        work_dir: Path, aligner: str, threads: int = 8) -> dict:
    """Align one run's reads against one taxid's pre-built index. Returns
    {'mapped_pct': float, 'unmapped_r1': Path, 'unmapped_r2': Path|None,
    'n_input': int, 'n_unmapped': int} or {} on failure. Mapped % always
    computed by direct read counting (input vs. unmapped-output), not by
    parsing either aligner's own stats text — keeps the two aligners directly
    comparable and doesn't depend on guessing an exact key format for either."""
    if not index_ready(taxid, aligner):
        return {}
    raw = (_align_against_taxid_hisat2(taxid, run, r1, r2, work_dir, threads)
           if aligner == "hisat2" else
           _align_against_taxid_bbmap(taxid, run, r1, r2, work_dir))
    if not raw:
        return {}
    out_u1, out_u2 = raw["unmapped_r1"], raw["unmapped_r2"]
    n_input = _count_fastq_reads(r1) + (_count_fastq_reads(r2) if r2 else 0)
    n_unmapped = _count_fastq_reads(out_u1) + (_count_fastq_reads(out_u2) if out_u2 else 0)
    mapped_pct = (n_input - n_unmapped) / n_input * 100 if n_input else 0.0
    return {"mapped_pct": mapped_pct, "unmapped_r1": out_u1, "unmapped_r2": out_u2,
            "n_input": n_input, "n_unmapped": n_unmapped}


def split_run(row: dict, work_dir: Path, aligner: str, threads: int = 8) -> dict:
    """Run stage 2 for one Run: align against every candidate host's index,
    pick the highest-mapping candidate as confirmed host, keep its unmapped
    reads as the final pathogen-enriched output. Returns a SPLIT_RESULTS row."""
    run = row["Run"]
    r1 = SELECT_READS_DIR / f"{run}_1.fastq.gz"
    r2 = SELECT_READS_DIR / f"{run}_2.fastq.gz"
    se = SELECT_READS_DIR / f"{run}.fastq.gz"
    paired = r1.exists() and r2.exists()
    if not paired and not se.exists():
        return {**{k: "" for k in SPLIT_RESULTS_COLS}, "Run": run,
                "BioSample": row.get("BioSample", ""), "status": "no_reads"}

    candidates = sorted({t.strip() for t in (row.get("candidate_host_taxids", "") or "").split(";") if t.strip()})
    results = {}
    for taxid in candidates:
        if paired:
            res = align_against_taxid(taxid, run, r1, r2, work_dir, aligner, threads)
        else:
            res = align_against_taxid(taxid, run, se, None, work_dir, aligner, threads)
        if res:
            results[taxid] = res

    if not results:
        return {**{k: "" for k in SPLIT_RESULTS_COLS}, "Run": run,
                "BioSample": row.get("BioSample", ""),
                "BioProject": row.get("BioProject", ""),
                "llm_host_resolved": row.get("llm_host_resolved", ""),
                "status": "no_index_available"}

    winner = max(results, key=lambda t: results[t]["mapped_pct"])
    win = results[winner]

    # Promote the winner's unmapped reads to the final output location.
    SPLIT_READS_DIR.mkdir(parents=True, exist_ok=True)
    final_r1 = SPLIT_READS_DIR / f"{run}_1.fastq.gz"
    win["unmapped_r1"].rename(final_r1)
    final_r2 = None
    if win.get("unmapped_r2"):
        final_r2 = SPLIT_READS_DIR / f"{run}_2.fastq.gz"
        win["unmapped_r2"].rename(final_r2)

    # Clean up the non-winning candidates' unmapped FASTQ (kept only stats).
    for t, res in results.items():
        if t == winner:
            continue
        for key in ("unmapped_r1", "unmapped_r2"):
            p = res.get(key)
            if p and p.exists():
                p.unlink()

    candidate_summary = "; ".join(
        f"{t}:{results[t]['mapped_pct']:.1f}%" for t in sorted(results))
    resolved_taxid = (row.get("host_taxid", "") or "").strip()
    agrees = "" if not resolved_taxid else str(resolved_taxid == winner)

    return {
        "Run": run, "BioSample": row.get("BioSample", ""),
        "BioProject": row.get("BioProject", ""),
        "llm_host_resolved": row.get("llm_host_resolved", ""),
        "confirmed_host_taxid": winner,
        "confirmed_host_mapped_pct": f"{win['mapped_pct']:.2f}",
        "agrees_with_llm": agrees,
        "candidate_results": candidate_summary,
        "status": "ok",
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-index", action="store_true",
                    help="Stop after resolving genomes + building indices; "
                         "don't run the split stage")
    ap.add_argument("--limit", type=int, default=None,
                    help="Process at most N runs in the split stage (for testing)")
    ap.add_argument("--aligner", choices=ALIGNERS, default="bbmap",
                    help="Alignment tool for host-read removal (default: bbmap). "
                         "bbmap CANNOT build for genomes with a chromosome "
                         ">500Mbp (wheat and likely other cohort outliers — see "
                         "module docstring); hisat2 has no such limit and is the "
                         "standard splice-aware RNA-seq choice. Indices for each "
                         "aligner live in separate subdirectories so both can be "
                         "tested without collision — but the final split reads "
                         "output path is shared regardless of --aligner, so "
                         "don't run both aligners' split stage back to back "
                         "without moving results in between.")
    ap.add_argument("--workers", type=int, default=4,
                    help="Parallel jobs for genome download, index build, and "
                         "split stages alike (default 4 — kept modest since "
                         "some host genomes are tens of Gb and each aligner job "
                         "gets a large memory/thread allocation of its own; too "
                         "much concurrency risks overcommitting node memory)")
    ap.add_argument("--threads", type=int, default=8,
                    help="Threads per index-build/alignment job (hisat2 only; "
                         "bbmap doesn't take a comparable flag here). Default 8.")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_run_split.log")
    link_latest(logs_base, log_dir / "kraken_run_split.log")
    sys.stdout = log

    try:
        if not RUN_LIST.exists():
            sys.exit(f"Error: {RUN_LIST} not found — run kraken_run_select.py --download first.")

        # ── stage 0: resolve + download host genomes ────────────────────────
        with open(RUN_LIST, newline="") as fh:
            select_rows = list(csv.DictReader(fh, delimiter="\t"))
        host_taxids_needed = set()
        for row in select_rows:
            for t in (row.get("candidate_host_taxids", "") or "").split(";"):
                t = t.strip()
                if t:
                    host_taxids_needed.add(int(t))
        print(f"Distinct host taxids named across {len(select_rows)} runs: "
              f"{len(host_taxids_needed)}", flush=True)

        taxid_to_accession = {int(k): v for k, v in load_json(HOST_TAXID_MAP).items()}
        already_done = {t for t in host_taxids_needed if t in taxid_to_accession}
        todo = host_taxids_needed - already_done
        print(f"Host genomes: {len(already_done)} already resolved (from a prior "
              f"run), {len(todo)} to fetch …", flush=True)
        if todo:
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                futs = {pool.submit(ensure_host_cds, tid): tid for tid in todo}
                for fut in as_completed(futs):
                    tid = futs[fut]
                    acc = fut.result()
                    taxid_to_accession[tid] = acc
                    print(f"  host taxid {tid}: {'OK ' + acc if acc else 'FAILED'}", flush=True)
            save_json({str(k): v for k, v in taxid_to_accession.items()}, HOST_TAXID_MAP)
        taxid_to_accession = {str(k): v for k, v in taxid_to_accession.items() if v}
        if not taxid_to_accession:
            sys.exit("Error: no host genomes resolved — nothing to index or split against.")

        # ── stage 1: build indices ─────────────────────────────────────────
        print(f"\nBuilding {args.aligner} indices for {len(taxid_to_accession)} "
              f"distinct host taxids ({args.workers} workers) …", flush=True)
        n_ok = n_cached = n_fail = 0
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(build_index, taxid, acc, args.aligner, args.threads): taxid
                    for taxid, acc in taxid_to_accession.items() if acc}
            for done, fut in enumerate(as_completed(futs), 1):
                taxid = futs[fut]
                status = fut.result()
                if status == "ok":
                    n_ok += 1
                elif status == "cached":
                    n_cached += 1
                else:
                    n_fail += 1
                if done % 10 == 0 or done == len(futs):
                    elapsed = time.time() - t0
                    print(f"  [{done}/{len(futs)}] built={n_ok} cached={n_cached} "
                          f"failed={n_fail}  ({elapsed:.0f}s)", flush=True)
        print(f"Index build complete: {n_ok} built, {n_cached} already cached, "
              f"{n_fail} failed", flush=True)

        if args.build_index:
            print("\n--build-index: stopping before the split stage.")
            return

        # ── stage 2: split ──────────────────────────────────────────────────
        rows = select_rows
        if args.limit:
            rows = rows[:args.limit]

        # Resume: skip runs whose split already completed successfully AND
        # whose promoted output still exists on disk. Only "ok" rows are
        # cacheable — failures always get retried, since a failure might get
        # fixed (real example: the 2026-09-07 HISAT2 filename bug) and
        # shouldn't be silently skipped forever as "already done". The disk
        # check (not just trusting the TSV row) protects against exactly the
        # mistake made the same day — deleting _tmp/ output and unknowingly
        # forcing a full, expensive realignment redo when a completed result
        # was already sitting there, only a filename fix away from usable.
        cached_results = {}
        if SPLIT_RESULTS.exists():
            with open(SPLIT_RESULTS, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("status") != "ok":
                        continue
                    final_r1 = SPLIT_READS_DIR / f"{r['Run']}_1.fastq.gz"
                    final_se = SPLIT_READS_DIR / f"{r['Run']}.fastq.gz"
                    if final_r1.exists() or final_se.exists():
                        cached_results[r["Run"]] = r

        rows_todo = [row for row in rows if row["Run"] not in cached_results]
        print(f"\nSplit: {len(cached_results)} already done (resumed, output "
              f"verified on disk), {len(rows_todo)} to process against their "
              f"candidate host {args.aligner} indices ({args.workers} "
              f"workers) …", flush=True)

        work_dir = DATA_DIR / "_tmp"
        work_dir.mkdir(parents=True, exist_ok=True)
        results = list(cached_results.values())
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(split_run, row, work_dir, args.aligner, args.threads): row["Run"]
                    for row in rows_todo}
            for done, fut in enumerate(as_completed(futs), 1):
                results.append(fut.result())
                if done % 10 == 0 or done == len(rows_todo):
                    elapsed = time.time() - t0
                    n_status = {}
                    for r in results:
                        n_status[r["status"]] = n_status.get(r["status"], 0) + 1
                    print(f"  [{done}/{len(rows_todo)}] {n_status}  ({elapsed:.0f}s)", flush=True)

        with open(SPLIT_RESULTS, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=SPLIT_RESULTS_COLS, delimiter="\t")
            w.writeheader()
            w.writerows(results)
        print(f"\nOutput: {SPLIT_RESULTS}")
        print(f"Next: python kraken/run/kraken_run_assign.py")
    finally:
        log.close()


if __name__ == "__main__":
    main()
