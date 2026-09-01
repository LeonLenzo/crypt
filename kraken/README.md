# kraken — Orthogonal Kraken2 Species-level Validation

## Rationale

The STAT-based screening pipeline relies on pre-computed k-mer profiles built from the entire SRA at NCBI scale. While this enables whole-corpus screening, it introduces a known reliability problem: k-mers shared across closely related species are promoted to the lowest common ancestor (LCA) node in the STAT database, preventing species-level resolution for taxonomically dense clades. A kallisto pilot on 45 MAL runs dominated by *Puccinia striiformis* f. sp. *tritici* (PST, wheat yellow rust) confirmed this directly — STAT assigned zero eukaryotic percentage to all 15 PST-dominated runs, while both kallisto and Kraken2 detected PST at 65–68% of classified reads. STAT resolves PST reads to the Dikarya node (kingdom-level), not species.

This raises a general concern: STAT detections are reliable where pathogen k-mer sets are sufficiently unique at SRA scale, but blind spots exist for organisms that share abundant k-mers with close relatives in the database. Kraken2 with a purpose-built, curated reference addresses this by controlling the database composition explicitly.

## Module layout

Two submodules, both run on Setonix (large data, needs HPC compute + storage):

```
kraken/
├── db/                    Submodule 1 — reference DB build (search → busco → build)
│   ├── kraken_db_search.py   Step 1/3 — select candidate assemblies (query NCBI,
│   │                         PHI-base seed pan-genome + genus fill-in) AND download
│   │                         their CDS FASTA. The ONE place in the DB pipeline that
│   │                         downloads.
│   ├── kraken_db_busco.py    Step 2/3 — BUSCO-score every candidate (reads CDS
│   │                         already on disk, never downloads), apply completeness
│   │                         thresholds.
│   ├── kraken_db_build.py    Step 3/3 — build the Kraken2 DB from BUSCO-selected
│   │                         assemblies (reads CDS already on disk, never downloads).
│   └── figures/               compare_host_pathogen.py, host_breakdown.py,
│                               busco_completeness.R — all DB-comparison scoped
├── run/                   Submodule 2 — read classification (select → split → assign)
│   ├── kraken_run_select.py  Step 1/3 — select target BioSamples from samples.tsv,
│   │                         download reads (prefetch/fasterq-dump) + host CDS.
│   ├── kraken_run_split.py   Step 2/3 — NOT YET BUILT. BBSplit host-read removal via
│   │                         a combined multi-reference index.
│   ├── kraken_run_assign.py  Step 3/3 — NOT YET BUILT. Kraken2 classification.
│   └── classify.py           Kraken2 classification (--run-list/--runs-tsv +
│                             --reads-dir or ENA streaming) — planned basis for
│                             kraken_run_assign.py (adapt, not rewrite).
├── benchmark_download.py Diagnostic — ENA FTP vs HTTPS vs S3 vs prefetch throughput
│                         (not submodule-specific, stays at kraken/ root).
├── legacy/                excluded from git — download.py + old SLURM wrapper (both
│                          superseded by kraken_run_select.py's built-in downloading)
└── slurm/                 SLURM wrappers (all currently db/-scoped)
```

`manifest.py` (repo root, shared — see `_util.py`'s `build_manifest()`/`upload_to_acacia()`)
records what's in a gitignored `data/` dir (path/size/mtime, +md5 for `.k2d` files) so huge
Setonix data stays visible from the repo without being tracked; kraken/ is its first user
but it's module-agnostic, not kraken-specific.

Output convention matches `stat/` and `metadata/`, with the db/run split mirrored:
`kraken/output/{db,run}/{search,busco,build,select,classify}/{data,logs}/`. Figure
scripts are the one exception — output stays in a shared `kraken/output/figures/`
collection dir regardless of which submodule the script lives in, matching how
`stat/output/figures/` and `metadata/output/figures/` work. Large data (CDS downloads,
BUSCO lineage caches, Kraken2 DBs, downloaded FASTQ) lives gitignored under `data/`;
small tracked TSVs (`ref_candidates.tsv`, `busco_scores.tsv`, `manifest.tsv`) live
alongside it in the same directory.

**2026-09-01 restructure**: split the flat `kraken/` script layout into `db/` + `run/`
subdirectories matching the two submodules (Leon's call — "let's have separate
kraken/db and kraken/run dirs"). Filenames kept their full `kraken_db_`/`kraken_run_`
prefixes (searchable/greppable as-is) even though now slightly redundant with the
parent directory name. Output directories were moved to mirror the split
(`kraken/output/kraken_db_search/` → `kraken/output/db/search/`, etc.) — this was
done locally and needs the equivalent move applied on Setonix once it's back from
maintenance (see the prepared move commands in memory/kraken_restructure_plan.md).

## Database design

The Kraken2 database contains only PHI-base eukaryotic pathogen CDS sequences (fungi and oomycetes). Host sequences are deliberately excluded — because the DB contains only pathogen sequences, `pct_classified` in Kraken2 output directly represents pathogen burden without requiring normalisation against host signal. Host reads simply go unclassified. Including host CDS would require the same Viridiplantae gate logic used in STAT MAL, reintroducing the ambiguity the orthogonal approach is meant to resolve.

### Genus fill-in (no masking)

Rather than masking to create species-diagnostic k-mers, the DB uses **genus-level completeness** to achieve specificity. For each PHI-base genus detected in STAT screening, the best-annotated assembly for every non-PHI-base species within that genus is also included. This genus fill-in means that k-mers shared between species within a genus resolve — via Kraken2's LCA algorithm — to the genus node rather than being spuriously assigned to the wrong species. Species-level detections therefore represent reads for which no within-genus congener claims the k-mer, providing a natural specificity filter without manual masking.

### Coverage and gaps

CDS-based sequences are used throughout rather than whole-genome sequences: RNA-seq reads derive from spliced transcripts and align poorly to genomic sequence across introns, and CDS sequences are shorter and more specific. Seeds lacking NCBI annotation (no gene model) are excluded — without CDS coordinates there is no transcript-compatible sequence to add.

## Building the DB (submodule 1)

Run from `crypt/` on Setonix, in order:

```bash
python kraken/db/kraken_db_search.py --scope          # optional: pangenome.tsv / genus_fill.tsv report
python kraken/db/kraken_db_search.py --download        # select candidates + download CDS
python kraken/db/kraken_db_busco.py                    # BUSCO score + threshold + fallback selection
python kraken/db/kraken_db_build.py                    # build the Kraken2 DB
```

Or via SLURM: `sbatch kraken/slurm/kraken_db_search.slurm` → `kraken_db_busco.slurm` →
`kraken_db_build.slurm`. Each step is resumable (accession-level caches); re-running
after a threshold change only needs `kraken_db_busco.py --finalize-only` (no re-scan).

**Current DB**: `db_v2` (`kraken/output/db/build/data/db_v2/`, ~20GB) — built
2026-08-16 from 1,017 BUSCO-selected assemblies (thresholds: fungal ≥ 50%, oomycete ≥
65%). The older `db_pathogens` pilot DB (pre-BUSCO-rebuild) was dropped 2026-08-28 once
`db_v2` was confirmed current and working.

## Running submodule 2 (select → split → assign)

```bash
python kraken/run/kraken_run_select.py --limit N --download   # select + download reads + host CDS
python kraken/run/kraken_run_split.py                         # NOT YET BUILT — BBSplit host removal
python kraken/run/kraken_run_assign.py                         # NOT YET BUILT — Kraken2 classification
```

Only step 1 exists so far. `kraken/run/classify.py` (`--run-list`/`--runs-tsv` +
`--reads-dir`, confidence=0.15, min-hit-groups=3, results append to
`kraken/output/run/classify/data/kraken_cache.jsonl`) is the planned basis for
`kraken_run_assign.py` — kept in the active tree for that reason, not yet wired into
the new flow.

**Validation target**: ALL 2,719 field/aerial BioSamples from
`metadata/figures/sample_funnel_v3.py` — corrected 2026-08-31 from an earlier draft
scope of "the 473 already flagged cryptic-co-infected" (restricting to already-flagged
samples would defeat the point of a validation pass meant to catch what STAT *missed*,
e.g. the documented PST/rust blind spot, which shows up as "clean" per STAT). 467 of
the 2,719 happen to already be flagged. Each comes with `pathogen_match_status`, named
vs STAT-detected pathogens, and a full manuscript.

`kraken/legacy/download.py` (the old ENA-curl-based downloader) is superseded by
`kraken_run_select.py`'s built-in prefetch/fasterq-dump downloading.

## Pilot results

A 45-run pilot (2026-08-05) followed by a 32-run high-confidence validation set (2026-08-07) confirmed:

- STAT correctly detects the majority of HAL-mode eukaryotic co-infections where the secondary organism has a sufficiently distinct k-mer profile
- The PST blind spot is confined to rust fungi (Pucciniales); non-rust fungal detections in the pilot set showed broad STAT–Kraken2 agreement
- The pathogen-only DB (no host sequences) produces cleaner `pct_classified` values than earlier host-inclusive DB pilots, where host background introduced noise at low pathogen burden thresholds

## Limitations

**CDS-only database misses splicing junctions.** Kraken2 k-mers span 35 bp; reads that bridge exon–exon junctions (absent from CDS sequences) will not match database k-mers, slightly reducing sensitivity. In practice this effect is small for highly expressed genes, where most reads fall within exons.

**Database completeness.** Seeds lacking annotation are not in the database. Co-infections involving unannotated or poorly assembled pathogens are invisible to Kraken2 but detectable by STAT (which uses the full SRA k-mer set including environmental and low-coverage organisms).

**Rust blind spot not fully characterised.** The STAT PST blind spot was confirmed by the pilot, but the extent to which other rust clades (e.g., *Melampsora*, *Hemileia*) share this problem has not been systematically tested.

## History

**2026-08-28 restructure**: consolidated `scope_db.py` + `screen_refs.py` + the CDS-download
step (previously duplicated across `busco_screen.py`/`build.py`) into `kraken_db_search.py`;
`busco_screen.py` + `filter_refs.py` into `kraken_db_busco.py` (now producing one merged
output table instead of two); `build.py` into `kraken_db_build.py` (fetch logic removed —
reads CDS that `kraken_db_search.py` already downloaded). Manifest generation and Acacia
upload (previously `kraken_manifest.py` and a function inside `build.py`) moved to the
shared `_util.py` + repo-root `manifest.py` — both are generic, not kraken-specific.
Old scripts, their old outputs, and dead debug/experiment leftovers (an abandoned
"unmasked" DB variant) were fully deleted rather than kept in a `legacy/` dir, once the
migration to the new scripts was verified to reproduce `db_v2`'s exact composition with
no data loss. Also dropped: the 2,473-run stratified SRA control set (`control/sample.py`)
and in silico controls (`control/insilico.py`, stratum J) — abandoned per Leon 2026-08-28,
those runs turned out to have unknown ground truth and weren't suitable as controls; the
2,719-BioSample field/aerial target (above) replaces them as the validation set.

## Key output files

| File | Contents |
|------|----------|
| `kraken/output/db/search/data/ref_candidates.tsv` | Candidate assemblies (seed pan-genome + genus fill-in), pre-BUSCO |
| `kraken/output/db/busco/data/busco_scores.tsv` | Merged: candidate metadata + BUSCO score + pass/fail + final `selected` decision |
| `kraken/output/db/build/data/db_v2/` | Current Kraken2 DB (gitignored; see `manifest.tsv` alongside it) |
| `kraken/output/run/classify/data/kraken_cache.jsonl` | Append-only classification cache; one JSON object per run |
| `kraken/output/run/select/data/run_list.tsv` | Submodule 2 target BioSamples/Runs + host resolution + download status |
