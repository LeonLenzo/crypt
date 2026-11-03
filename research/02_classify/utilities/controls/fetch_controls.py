#!/usr/bin/env python3
"""fetch_controls.py — pull a negative-control study's reads from ENA by accession.

    python 02_classify/utilities/controls/fetch_controls.py --accession PRJEB29046 --dry-run
    python 02_classify/utilities/controls/fetch_controls.py --accession PRJEB29046 PRJNA613154
    python 02_classify/utilities/controls/fetch_controls.py --meta-only --accession SRP117269

## Why this exists and is not a flag on fetch_reads.py

`fetch_reads.py` resolves what to download from `run_list.tsv` and `run_depth.tsv`, which are
THE COHORT: 7,688 curated runs whose URLs were cached one request per BioProject during
`03_cohort/resolve.py`. A control study is not in either table and never will be, because it is
not a cohort sample and must not be counted as one. So the plan has to come from the archive at
fetch time, which is the one real difference. Everything after that -- the retry schedule, the
HTTP_RETRY/HTTP_FATAL split, stall detection, verify-then-rename -- is IMPORTED from
`fetch_reads.fetch_run` rather than reimplemented. That code cost a 623-run failure to get
right and a second copy of it would drift.

## Reads land in a SEPARATE tree

`_layout.CONTROL_READS`, not `READS`. See the comment there: `scan_tree(READS)` is what
04_classify and the reads manifest treat as "the cohort on disk", and a control accession
sitting in that tree would enter the denominator of every question about how much of the cohort
is fetched. The layout inside is identical, so `classify.py --reads/--reports` runs over it
with no special case.

## The archive's answer is CAPTURED, not trusted

Every filereport response is written to `CONTROL_META/<accession>.tsv` as served. That is not
politeness about provenance: PRJEB29046's ENA sample title for ERR2822761 reads "4 days post
control inoculation", while the paper's own data-availability section names that run Sample P2,
the 4 DPI *P. tritici-repentis* infection. ENA has the mock's title on the infected run. The
roles in `01_literature/controls/control_template_runs.tsv` follow the PAPER, as the audit
chain requires, and classification settles it independently -- a mock library has no Ptr.
Keeping the capture keeps the disagreement visible. Same lesson as
`01_literature`'s: the archive is a source, not an authority.

## What is checked before anything is downloaded

A control is only useful if it is the chemistry the claim is made in, so the metadata pass
prints `library_strategy`, `library_source` and `library_selection` per run and REFUSES any run
whose `library_source` is not TRANSCRIPTOMIC. That filter is not hypothetical: 245 of
*P. triticina*'s 418 archive runs labelled "RNA-Seq" are `library_source=GENOMIC`, and a
genomic library would read as a false positive factory against a CDS database. Non-Illumina
and Sanger rows are dropped the same way -- PRJNA560644 carries six Sanger sequences beside its
twelve RNA-seq runs.
"""

import argparse
import csv
import subprocess
import sys
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _layout import (CONTROL_FETCH_STATUS, CONTROL_META, CONTROL_READS,  # noqa: E402
                     CONTROL_READS_TMP)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "03_cohort"))
from fetch_reads import fetch_run  # noqa: E402

# One request per accession, exactly as 01_literature's ENA lesson requires: the portal serves
# a whole study in one call, and per-sample requests were ~150x slower and returned less.
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ("run_accession,study_accession,secondary_study_accession,sample_accession,"
          "library_strategy,library_source,library_selection,instrument_platform,"
          "instrument_model,read_count,base_count,fastq_bytes,fastq_ftp,"
          "scientific_name,sample_title,experiment_title")

# ── THE MOLECULE CHECK, AND WHY IT IS NOT library_source ───────────────────────────────────
# This filter answers exactly one question: is this library RNA? A CDS database cannot be
# scored on DNA. It does NOT decide whether a study is an admissible control -- that is a
# methods read, and it lives in `01_literature/controls/control_candidates.tsv`. Keeping the
# two apart is the point: a label can be wrong, a quoted methods sentence cannot.
#
# It keyed on `library_source == TRANSCRIPTOMIC` until 2026-10-27 and that rejected a valid
# control outright. **PRJEB58154 (Child 2023, Z. tritici, the highest undeclared load at 343)
# has all 18 runs labelled `library_source=METAGENOMIC`, `scientific_name=plant metagenome`,
# with `library_strategy=RNA-Seq`.** The submitter described the library's CONTENT -- a wheat
# leaf carrying a fungus is two organisms -- not its molecule. The paper is unambiguous: sample
# labels are strain x dpi x replicate ("B6_2" = Dztbck1, 6 dpi, rep 2) for three single strains
# inoculated separately, and libraries were TruSeq HT stranded mRNA, so poly(A) and matching
# the cohort's chemistry. Rejecting it lost the single most valuable control to a metadata
# convention.
#
# STRATEGY ALONE IS ALSO NOT ENOUGH, which is why both are tested. The original justification
# stands and is the opposite error: 245 of P. triticina's 418 runs declare
# `library_strategy=RNA-Seq` with `library_source=GENOMIC`. So the rule is RNA-Seq strategy AND
# not a DNA source.
WANT_STRATEGY = "RNA-SEQ"
DNA_SOURCES = {"GENOMIC", "METAGENOMIC_DNA", "SYNTHETIC", "VIRAL RNA"}
# Admissible, but the second two describe a library that MAY be a community rather than one
# pathogen plus its host. That cannot be settled from a label, so they are admitted and FLAGGED
# for the methods read instead of being silently taken or silently dropped.
RNA_SOURCES = {"TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "METAGENOMIC"}
FLAG_SOURCES = {"METATRANSCRIPTOMIC", "METAGENOMIC"}
WANT_PLATFORM = "ILLUMINA"

STATUS_COLS = ["Run", "BioProject", "status", "seconds", "bytes", "n_files", "note"]

# Modest, and deliberately below fetch_reads' 24. The controls are small and are fetched
# ALONGSIDE or BEHIND a cohort fetch that already holds 24 streams; job 51341245 ended with EBI
# refusing connections after ~13 h of sustained concurrency, and the cause was never pinned
# down. A control batch is tens of GB, so there is nothing to gain by crowding that.
DEFAULT_WORKERS = 6


def filereport(accession: str, refresh: bool = False) -> list[dict]:
    """The archive's own rows for one accession, captured to CONTROL_META and reused."""
    CONTROL_META.mkdir(parents=True, exist_ok=True)
    cache = CONTROL_META / f"{accession}.tsv"
    if cache.exists() and cache.stat().st_size > 0 and not refresh:
        text = cache.read_text()
    else:
        q = urllib.parse.urlencode({"accession": accession, "result": "read_run",
                                    "fields": FIELDS, "format": "tsv"})
        r = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "120",
                            f"{ENA_API}?{q}"], capture_output=True, text=True)
        if r.returncode != 0 or not r.stdout.strip():
            sys.exit(f"{accession}: ENA filereport failed -- "
                     f"exit {r.returncode}: {(r.stderr or '').strip()[:200]}")
        text = r.stdout
        # Write only AFTER a non-empty response, so a failed call never leaves a cache that
        # looks like "this study has no runs". An accession that genuinely has none is a real
        # finding (C. zeae-maydis / PRJNA1463572) and must come back from the API every time,
        # not from a file.
        cache.write_text(text)
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        sys.exit(f"{accession}: ENA returned no runs. If the paper cites this accession, that "
                 f"is an ABSENCE_EVIDENCE finding, not a fetch problem.")
    return rows


def admissible(row: dict) -> str:
    """"" if the run may be fetched, else the reason it may not. See THE MOLECULE CHECK."""
    strat = (row.get("library_strategy") or "").upper()
    src = (row.get("library_source") or "").upper()
    if strat != WANT_STRATEGY:
        return f"library_strategy={row.get('library_strategy') or 'empty'}"
    if src in DNA_SOURCES:
        return f"library_source={src} is DNA"
    if src not in RNA_SOURCES:
        return f"library_source={src or 'empty'} is not a known RNA source"
    if (row.get("instrument_platform") or "").upper() != WANT_PLATFORM:
        return f"platform={row.get('instrument_platform') or 'empty'}"
    if not row.get("fastq_ftp"):
        return "no ENA fastq"
    return ""


def flagged(row: dict) -> str:
    """"" unless the source label leaves the single-pathogen question open."""
    src = (row.get("library_source") or "").upper()
    if src in FLAG_SOURCES:
        return (f"{src}: RNA, so fetchable, but the label does not say whether this is ONE "
                f"pathogen plus its host or a community -- confirm against the methods")
    return ""


def plan(accession: str, refresh: bool) -> tuple[list, list, set]:
    """(todo, rejected, notes) for one accession. todo rows are (run, urls, sizes)."""
    todo, rejected, notes = [], [], set()
    for r in filereport(accession, refresh):
        why = admissible(r)
        if why:
            rejected.append((r["run_accession"], why))
            continue
        if flagged(r):
            notes.add(flagged(r))
        urls = [u for u in r["fastq_ftp"].split(";") if u]
        sizes = [int(s) for s in (r.get("fastq_bytes") or "").split(";") if s.strip().isdigit()]
        sizes += [0] * (len(urls) - len(sizes))
        todo.append((r["run_accession"], urls, sizes))
    return sorted(todo), rejected, notes


def describe(accession: str, refresh: bool) -> None:
    """Print what the archive says, run by run. The pass to read BEFORE spending the stream."""
    rows = filereport(accession, refresh)
    print(f"\n{accession}: {len(rows)} runs")
    print(f"  {'run':14s} {'source':18s} {'strategy':10s} {'GB':>6s} {'Mreads':>7s}  title")
    for r in sorted(rows, key=lambda x: x["run_accession"]):
        gb = sum(int(s) for s in (r.get("fastq_bytes") or "").split(";")
                 if s.strip().isdigit()) / 1e9
        mr = int(r["read_count"]) / 1e6 if (r.get("read_count") or "").isdigit() else 0
        why, flag = admissible(r), flagged(r)
        print(f"  {r['run_accession']:14s} {(r.get('library_source') or '-'):18s} "
              f"{(r.get('library_strategy') or '-'):10s} {gb:6.2f} {mr:7.1f}  "
              f"{(r.get('sample_title') or r.get('experiment_title') or '')[:48]}"
              f"{'   REJECTED: ' + why if why else ('   FLAG' if flag else '')}")
    notes = {flagged(r) for r in rows if flagged(r)}
    for n in sorted(notes):
        print(f"  NOTE  {n}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", nargs="+", required=True,
                    help="ENA/SRA study accessions. PRJ, SRP and ERP all resolve")
    ap.add_argument("--meta-only", action="store_true",
                    help="query and print, download nothing. Read this before every batch")
    ap.add_argument("--refresh", action="store_true", help="re-query ENA, ignore the capture")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if args.meta_only:
        for acc in args.accession:
            describe(acc, args.refresh)
        return

    CONTROL_READS_TMP.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.time()
    for acc in args.accession:
        todo, rejected, notes = plan(acc, args.refresh)
        tb = sum(sum(s) for _, _, s in todo) / 2**40
        print(f"\n--- {acc}: {len(todo)} runs, {tb:.4f} TiB"
              f"{'  (DRY RUN)' if args.dry_run else ''} ---", flush=True)
        for run, why in rejected:
            print(f"  {run}: REJECTED, {why}")
        for n in sorted(notes):
            print(f"  NOTE  {n}")
        if not todo:
            print(f"  {acc}: nothing admissible")
            continue
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(fetch_run, run, acc, urls, sizes, args.dry_run,
                              CONTROL_READS, CONTROL_READS_TMP)
                    for run, urls, sizes in todo]
            for f in as_completed(futs):
                rows.append(f.result())

    if not args.dry_run and rows:
        CONTROL_FETCH_STATUS.parent.mkdir(parents=True, exist_ok=True)
        new = not CONTROL_FETCH_STATUS.exists()
        with open(CONTROL_FETCH_STATUS, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=STATUS_COLS, delimiter="\t", lineterminator="\n")
            if new:
                w.writeheader()
            w.writerows(rows)

    ok = sum(1 for r in rows if r["status"] in ("ok", "cached"))
    bad = [r for r in rows if r["status"] not in ("ok", "cached", "would-fetch")]
    gb = sum(r["bytes"] for r in rows) / 2**30
    print(f"\n{ok}/{len(rows)} ok, {len(bad)} failed, {gb:.1f} GiB, "
          f"{(time.time() - t0) / 60:.1f} min")
    for r in bad[:10]:
        print(f"  FAILED {r['Run']}: {r['status']} {r['note']}")
    parked = sum(f.stat().st_size for f in CONTROL_READS_TMP.rglob("*") if f.is_file())
    if parked:
        print(f"  {parked / 2**30:.1f} GiB parked in {CONTROL_READS_TMP.name}/ -- a resubmit "
              f"RESUMES these")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
