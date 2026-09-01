# Pipeline Schematic — crypt co-infection mining pipeline

> Regenerated 2026-09-01 (previous version 2026-08-17 was stale across nearly every
> module — metadata rename, stat rename, kraken restructure submodules 1+2, dropped
> control-set subpipeline, and the 2026-08-31 code-not-data gitignore policy all
> postdate it). All paths relative to `crypt/` root unless stated. Run all scripts
> from the `crypt/` root directory.

---

## Shared root

### `_util.py`
- **Inputs:** none
- **Outputs:** none (library only)
- **Provides:** `_Tee` (stdout tee to log file), `make_log_dir` / `link_latest` (timestamped log management), `http_get` (GET with backoff), `load_json` / `save_json` (atomic JSON I/O), `build_manifest` (data-dir → manifest.tsv), `resolve_taxon_name` (+ `_clean_taxon_name`/`_taxon_candidates` helpers, `HOST_NAME_ALIASES` ~50-entry common→scientific name table), `upload_to_acacia` (S3 sync)
- **Used by:** every Python script in stat/, kraken/, and metadata/

### `manifest.py`
- **Inputs:** `--data-dir PATH` (required), `--out PATH` (default `<data-dir>/manifest.tsv`), `--checksum-suffix` (repeatable, default `.k2d`)
- **Outputs:** `<data-dir>/manifest.tsv` — path/size/mtime, +md5 for checksum-suffix matches
- **Calls:** `_util.build_manifest()`
- **Purpose:** run against any gitignored bulk-data dir so the repo shows what exists (on local disk or Setonix scratch) without holding the data itself. Every `stat/`, `metadata/`, and `kraken/` output `data/` dir with real gitignored content now has one (added 2026-09-01 for stat/metadata — kraken/ already had this from the 2026-08-28 restructure).

---

## Module: stat/

Three scripts in strict sequential order. All paths use the `stat/output/` tree.

---

### `stat/stat_build.py`
- **Inputs:**
  - PHI-base CSV (`--phibase`, default `stat/output/stat_build/data/phi-base_current.csv`, auto-downloaded if absent/`--fetch`)
  - Local `ete3.NCBITaxa` SQLite DB
- **Outputs:**
  - `stat/output/stat_build/data/phibase_db.json` — kingdom-separated taxid maps (gitignored — code-not-data policy)
  - `stat/output/stat_build/logs/history/<timestamp>/{build.log,build_summary.txt}` (+ `.latest` symlinks — now tracked)
- **Calls:** `ete3.NCBITaxa` only, no subprocess
- **Feeds into:** `stat/stat_fetch.py`, `stat/stat_filter.py`, `kraken/db/kraken_db_search.py`, `metadata/meta_classify.py`, `kraken/db/figures/compare_host_pathogen.py`
- **Note:** `ictv_vmr.xlsx` still sits in `data/` but is now orphaned — only the old `stat/legacy/build.py` read it (viruses dropped from scope in the 2026-08-03 eukaryotic-only pivot)

---

### `stat/stat_fetch.py`
- **Inputs:** `--mode {mal|hal}` (required); `stat/output/stat_build/data/phibase_db.json`; NCBI Entrez (SRA RunInfo); NCBI STAT endpoint; resumes from its own checkpoints
- **Outputs:**
  - `stat/output/stat_fetch/data/stat_cache.jsonl` — unified RunInfo+STAT cache, append-only (2.4GB, gitignored)
  - `stat_cache_index.txt`, `{mode}_uid_checkpoint.txt`, `{mode}_accessions.txt`
  - `stat/output/stat_fetch/logs/history/<timestamp>/{mode}.log` (+ summary, now tracked)
- **Calls:** `urllib`-based HTTP only (via `_util.http_get`), no subprocess
- **Caution:** do NOT run MAL and HAL simultaneously — both write `stat_cache.jsonl`
- **Feeds into:** `stat/stat_filter.py`, `stat/figures/prep_scatter.py`, `kraken/db/figures/compare_host_pathogen.py`

---

### `stat/stat_filter.py`
- **Inputs:** `--mode {mal|hal|both}`, `--skip-validate`; `stat/output/stat_build/data/phibase_db.json`; `stat/output/stat_fetch/data/stat_cache.jsonl`; `{mode}_accessions.txt`
- **Outputs:**
  - `stat/output/stat_filter/data/runs.tsv` — 10,995 rows (MAL 5,223 + HAL 5,772 gate-pass); key cols per CLAUDE.md's Output schemas section (gitignored — code-not-data policy)
  - `{mode}_kingdom_dist.tsv`, `{mode}_species_dist.tsv` (unless `--skip-validate`)
  - `stat/output/stat_filter/logs/history/<timestamp>/...` (now tracked)
- **Calls:** none (pure Python over already-local data)
- **Key algorithm:** `specific_hits()` — leaf-level species detection via count nesting; `KINGDOM_THRESHOLDS`: Fungi ≥0.5%, Oomycota ≥0.5%, Nematoda ≥1.0%; MAL gate: Viridiplantae ≥1%; HAL gate: any euk pathogen ≥1%
- **Feeds into:** the single most-consumed file in the pipeline — `stat/figures/prep_scatter.py`, `metadata/meta_search.py`, `metadata/meta_classify.py`, `kraken/db/kraken_db_search.py`, `kraken/run/kraken_run_select.py`, `kraken/run/classify.py` (via `--runs-tsv`)

---

### `stat/figures/prep_scatter.py`
- **Inputs:** `stat/output/stat_fetch/data/stat_cache.jsonl`, `stat/output/stat_filter/data/runs.tsv`, `stat/output/stat_fetch/data/{mal,hal}_accessions.txt` — genuinely 3 distinct sources, not reducible to `runs.tsv` alone: `runs.tsv` only holds gate-*pass* runs, so the plot's background/grey (gate-*fail*) points require the raw `stat_cache` + accession files to know each non-gate-pass run's mode. Even for gate-pass rows, host_pct/euk_pct are deliberately recomputed fresh from `stat_cache` rather than trusted from `runs.tsv`, so the gate lines are guaranteed consistent between foreground and background.
- **Outputs:** `stat/output/figures/scatter/scatter_data.tsv` (gitignored; layer/mode/coinf/host_pct/euk_pct)
- **Feeds into:** `stat/figures/scatter.R`

### `stat/figures/scatter.R`
- **Inputs:** `stat/output/figures/scatter/scatter_data.tsv`
- **Outputs:** `stat/output/figures/scatter/{scatter.pdf,scatter.png,scatter_caption.txt}` (tracked)
- **Calls:** ggplot2, dplyr, patchwork

### `stat/diagnostics/` — deleted 2026-09-01
One-off investigation scripts (`diag_mal_gate.py`, `benchmark_strategies.py`,
`check_primary_alignment.py`, `check_refseq_coverage.py`) that had already served
their purpose and referenced stale pre-rename paths. Deleted outright, not moved to
legacy — no archival value, findings already baked into the pipeline's design.

### `stat/migrate_cache.py` — deleted 2026-09-01
One-off migration script (old `stat_cache.jsonl` array format → new object format),
already run and complete. Its own docstring flagged the fallback code in
`stat_fetch.py` as "can be deleted once stat_cache has been fully rebuilt in the new
format" — that fallback function was removed too.

### `stat/legacy/` — superseded scripts, excluded from git entirely (see `.gitignore`)

---

## Module: metadata/

**The whole four-script pipeline this doc previously described
(`ncbi_metadata.py`/`web_metadata.py`/`undermind.py`/`classify_metadata.py`) is gone
from `metadata/` root** — all moved to `metadata/legacy/`, replaced by a new
three-script pipeline below (rebuilt 2026-08-18–26, committed 2026-08-31). Undermind
was dropped entirely, not just moved ("not going to be used moving forwards").

---

### `metadata/meta_search.py`
- **Inputs:** `stat/output/stat_filter/data/runs.tsv`; NCBI Entrez (BioProject/PubMed/PMC), EBI EuropePMC, CrossRef, Serper (needs `SERPER_API_KEY`), EBI BioSamples API (ENA/DDBJ `SAME*`/`SAMD*`, empty from NCBI efetch)
- **Outputs:** `metadata/output/meta_search/data/{bioprojects.json, biosamples.json}` (gitignored, code-not-data) + `{bp_cache.json, bs_cache.json, serper_cache.json}` (gitignored, API caches)
- **Replaces:** `ncbi_metadata.py` + `web_metadata.py`
- **Feeds into:** `metadata/meta_text.py`, `metadata/meta_classify.py`, `metadata/figures/prep_lit_resolution.py`

### `metadata/meta_text.py`
- **Inputs:** `metadata/output/meta_search/data/bioprojects.json`; Unpaywall, PMC OA full-text, manual PDF drop (`metadata/output/meta_text/data/manual_pdfs/` — 715MB copyrighted journal PDFs, never committed)
- **Outputs:** `text_cache.jsonl` (gitignored), `failed_dois.tsv` (tracked); `--apply` writes `full_text` back into `bioprojects.json`
- **Calls:** `pdfminer.high_level.extract_text` (hard dependency)
- **Feeds into:** `metadata/figures/prep_lit_resolution.py`, `metadata/meta_classify.py`

### `metadata/meta_classify.py`
- **Inputs:** `metadata/output/meta_search/data/{bioprojects.json,biosamples.json}`; `stat/output/stat_filter/data/runs.tsv`; `stat/output/stat_build/data/phibase_db.json`; OpenAI API (`gpt-4o-mini`, `OPENAI_API_KEY`)
- **Outputs:** `metadata/output/meta_classify/data/samples.tsv` — **PRIMARY ANALYSIS INPUT**, one row per `biosample_representative` BioSample, full-text BioProjects only (gitignored); `classify_cache.jsonl`, `host_disambig_cache.jsonl` (gitignored)
- **Replaces:** `classify_metadata.py`
- **CLI:** `--focus {stress|setting|tissue|coinfection|hostpath|extract}`, `--disambiguate-hosts`, `--rerun-all`, `--workers N`
- **Feeds into:** `metadata/figures/sample_funnel_v3.py`, `kraken/run/kraken_run_select.py`
- **Note (2026-09-01):** `host_disambig_cache.jsonl` records now store the exact `named_hosts` candidate list each `host_resolved` value was picked from — a cached resolution only counts as done if the current candidates still match (fixed a real bug where an improved hostpath prompt's better scientific-name candidates silently went unused for already-resolved BioSamples: 46% of 422 previously-attempted entries were stale on inspection).

---

### `metadata/legacy/classify_metadata.py` — moved to legacy 2026-09-01
Older 3-call classifier (stress/setting/tissue only, no full-text gate, no host/pathogen
extraction). Superseded by `meta_classify.py`. Its only consumer, `sample_funnel_v2.py`,
moved to legacy alongside it (see below).

### `metadata/figures/sample_funnel_v3.py` — current primary figure
- **Inputs:** `metadata/output/meta_classify/data/samples.tsv`
- **Outputs:** `metadata/output/figures/sankey/sample_funnel_v3.{html,svg,png}`
- **Calls:** plotly (+kaleido)
- **Supersedes:** `sample_funnel.py` (already gone before this survey), `sample_funnel_v2.py` (moved to `metadata/legacy/figures/` 2026-09-01, read from `classify_metadata.py`'s dead-end output)

### `metadata/figures/prep_lit_resolution.py` + `metadata/figures/lit_resolution_alluvial.R` — current literature-resolution figure
- **Inputs:** `bioprojects.json`, `text_cache.jsonl` → `metadata/output/figures/sankey/lit_resolution_data.tsv` → `lit_resolution_alluvial.{png,svg}` (ggplot2 + ggalluvial)
- **Supersedes:** `metadata/figures/lit_resolution_sankey.py` (Plotly) — explicitly called "retired" in `prep_lit_resolution.py`'s own docstring; moved to `metadata/legacy/figures/` 2026-09-01. CLAUDE.md's active task list previously listed building the Plotly version as an open TODO — fixed 2026-09-01, it was stale (the R version already superseded it).

### `metadata/tools/export_review_lists.py`, `metadata/tools/review_designs.py` — deleted 2026-09-01
Both hardcoded a pre-restructure `output/`/`scripts/` path tree that no longer existed
anywhere in the repo — genuinely broken (hard `FileNotFoundError`), not just stale.
Deleted rather than fixed; would have needed a full rewrite either way.

### `metadata/legacy/` — excluded from git entirely, confirmed present locally (ncbi_metadata.py, web_metadata.py, undermind.py, fetch_xml.py, fetch_lit.py, filter_kw.py, llm_classify.py, serper_dump/resolve/scrape.py, backfill_doi.py, tissue_vocab.py, + legacy figures)

---

## Module: kraken/

**Two submodules.** Submodule 1 (DB build) is fully restructured and complete.
Submodule 2 (run/select/split/assign) is at step 2 of 3 — `kraken_run_select.py`
is built and (per 2026-09-01 testing) runs; `kraken_run_split.py` is written but not
yet tested against real data (Setonix down for maintenance); `kraken_run_assign.py`
does not exist yet. `kraken/control/` (the old stratified control-set subpipeline) is
confirmed **fully removed** — dropped 2026-08-28, not just deprecated.

---

### `kraken/db/kraken_db_search.py` — submodule 1, step 1/3
- **Inputs:** `stat/output/stat_build/data/phibase_db.json`, `stat/output/stat_filter/data/runs.tsv`; NCBI `datasets` CLI
- **Outputs:** `kraken/output/db/search/data/{ref_candidates.tsv, pangenome.tsv, genus_fill.tsv (--scope only)}` (tracked) + `cds/pathogen/{accession}/` (gitignored, `--download` only)
- **Calls:** `datasets summary/download` (subprocess), `unzip`
- **Feeds into:** `kraken/db/kraken_db_busco.py`, `kraken/db/kraken_db_build.py`, `kraken/run/kraken_run_select.py` (imports functions directly)

### `kraken/db/kraken_db_busco.py` — submodule 1, step 2/3
- **Inputs:** `ref_candidates.tsv`; CDS already on disk at `cds/pathogen/`; BUSCO Singularity + lineage DBs (Setonix)
- **Outputs:** `kraken/output/db/busco/data/{busco_scores.tsv, busco_scan_cache.tsv}` — thresholds fungi ≥50%, oomycete ≥65% (per-taxid fallback)
- **Calls:** BUSCO via `subprocess.Popen`
- **Feeds into:** `kraken/db/kraken_db_build.py`, `kraken/db/figures/busco_completeness.R`

### `kraken/db/kraken_db_build.py` — submodule 1, step 3/3
- **Inputs:** `busco_scores.tsv` (rows `selected=True`); CDS from `cds/pathogen/`; NCBI taxdump (auto-downloaded)
- **Outputs:** Kraken2 DB directory — current production `db_v2` (20GB, gitignored, `kraken/output/db/build/data/db_v2/`); optional `--upload-to-acacia`
- **Calls:** `kraken2-build --add-to-library` / `--build`; `_util.upload_to_acacia` → `aws s3 sync`
- **Feeds into:** `kraken/run/classify.py` (via `--db` at run time, no code coupling)

---

### `kraken/run/kraken_run_select.py` — submodule 2, step 1/3 (only step built)
- **Inputs:** `metadata/output/meta_classify/data/samples.tsv`, `stat/output/stat_filter/data/runs.tsv`; imports from `kraken_db_search.py`; requires `prefetch`/`fasterq-dump` (sra-tools) + `datasets` CLI
- **Outputs:** `kraken/output/run/select/data/run_list.tsv` (tracked), `.../data/host_taxid_to_accession.json` (tracked — every candidate taxid's downloaded accession, added 2026-09-01 so `kraken_run_split.py` can find/build an index for every named candidate, not just the resolved host), `.../data/reads/{run}_{1,2}.fastq.gz` (gitignored), `kraken/output/db/search/data/cds/host/{accession}/` (gitignored, shared pool with pathogen CDS)
- **CLI:** `--setting`, `--aerial-only`/`--no-aerial-only`, `--cryptic-only`, `--limit N`, `--download`, `--workers N`
- **Feeds into:** `kraken/run/kraken_run_split.py`
- **Status (2026-09-01):** smoke-tested with `--limit 2 --download` on Setonix; host resolution worked, but the host genome fetch step hit a 300s timeout — leading theory is it ran on the Setonix login node rather than through SLURM (untested at time of writing, Setonix down for maintenance)

### `kraken/run/kraken_run_split.py` — written 2026-09-01, NOT YET TESTED on Setonix
Design agreed with Leon (full rationale in kraken/README.md's "kraken_run_split.py
design" section): one `bbmap.sh` index per host taxid (94 individual builds, not one
combined index — ~275Gb total across all host genomes, dominated by outliers like
wheat/pine/oat that can't just be excluded since they're the dominant crops in the
cohort — even excluding only genomes >10Gb would drop coverage for 45% of the
cohort). Per run, each candidate host is aligned separately; highest mapping rate
(computed by direct read counting, not by parsing BBMap's statsfile text) wins as
confirmed host, its unmapped reads become the final output; other candidates kept
as QC metadata only.
- **Inputs:** `kraken/output/run/select/data/{run_list.tsv, host_taxid_to_accession.json, reads/}`, `kraken/output/db/search/data/cds/host/{accession}/*.fna`
- **Outputs:** `kraken/output/run/split/data/index/{taxid}/` (gitignored), `.../data/reads/{run}_{1,2}.fastq.gz` (gitignored, non-host reads), `.../data/split_results.tsv` (tracked)
- **CLI:** `--build-index` (stop after indexing), `--limit N`, `--workers N` (default 4 — deliberately modest, each `bbmap.sh` job gets its own `-Xmx48g`)
- **Calls:** `bbmap.sh` (module `bbmap/38.96--h5c4e2a8_0` on Setonix, no bootstrap needed)
- **Feeds into:** `kraken/run/kraken_run_assign.py` (not yet built)
- **Not yet verified**: never run against real data — Setonix down for maintenance when written. The read-counting mapped-% logic and the winner-selection behavior need a real smoke test (`--limit N`) before trusting them at scale.

### `kraken/run/kraken_run_assign.py` — not built
Kraken2 classification on `kraken_run_split.py`'s non-host reads; planned as a light
adaptation of `kraken/run/classify.py` (already supports `--reads-dir` for
pre-downloaded local reads), not a rewrite.

---

### `kraken/run/classify.py`
- **Inputs:** `--mode {mal|hal|both}` / `--runs-tsv` / `--run-list`; `--reads-dir` (pre-downloaded) or streams from ENA FTP; `--db PATH` (required)
- **Outputs:** `kraken/output/run/classify/data/{kraken_cache.jsonl, kraken_cache_index.txt}` (gitignored), optional per-run `.kreport` via `--reports-dir`
- **Calls:** `kraken2` (subprocess, confidence=0.15, min-hit-groups=3)
- **Feeds into:** nothing currently — no figure script reads `kraken_cache.jsonl` yet; `kraken/db/figures/compare_host_pathogen.py` instead expects a manually-scp'd `/tmp/setonix_kraken_metrics.json`
- **Status (2026-09-01):** requires `--run-list`/`--runs-tsv` explicitly now — the dead `IN_DIR` fallback landmine was removed. Kept in the active tree (not legacy) since it's the planned basis for `kraken_run_assign.py` ("light adaptation, not a rewrite" — already supports `--reads-dir` for pre-downloaded local reads, which is exactly what `kraken_run_select.py` produces).

### `kraken/legacy/download.py` — moved to legacy 2026-09-01
Superseded by `kraken_run_select.py`'s built-in downloading (prefetch/fasterq-dump)
— nothing in the current select→split→assign flow calls this anymore. Its SLURM
wrapper (`download_runs.slurm`) moved to legacy alongside it. Excluded from git
per the `kraken/legacy/` convention (matches `stat/legacy/`, `metadata/legacy/`).

### `kraken/benchmark_download.py`
- **Inputs:** `--methods {ena-ftp,ena-https,sra-s3,prefetch}`, `--run-list` (else falls back to dead `kraken/control/output/data/control_runs.tsv` constant)
- **Outputs:** results TSV per `--out`
- **Purpose:** throughput benchmarking, diagnostic only

---

### `kraken/db/figures/` — paths fixed, scripts relocated here 2026-09-01

Originally at `kraken/figures/` with hardcoded pre-rename paths (fixed same session);
moved into `kraken/db/figures/` when the db/run submodule split happened, since all
three are DB-comparison scoped, not run/-scoped:

| Script | Was-broken input path (now fixed) |
|---|---|
| `compare_host_pathogen.py` | `stat/output/build/data/phibase_db.json` → `stat/output/stat_build/data/phibase_db.json` |
| `host_breakdown.py` | `stat/output/filter_runs/data/runs.tsv` → `stat/output/stat_filter/data/runs.tsv` |
| `busco_completeness.R` | `kraken/busco_scores.tsv` (relative, didn't exist) → `kraken/output/db/busco/data/busco_scores.tsv` |

`compare_host_pathogen.py` also depends on `/tmp/setonix_kraken_metrics.json`, a
manually-scp'd file not produced by any tracked script. Figure *output* stays at the
shared `kraken/output/figures/{db_comparison,busco,host_breakdown}/` collection dir
(same convention as `stat/output/figures/`, `metadata/output/figures/`) — only the
script locations moved, not where they write.

---

### `kraken/slurm/` — cleaned up 2026-09-01

| Script | Calls | Status |
|---|---|---|
| `kraken_db_search.slurm` | `kraken_db_search.py --download --workers 16` | current (16 CPU/32G/8h) |
| `kraken_db_busco.slurm` | `kraken_db_busco.py` | current (128 CPU/230G/24h) |
| `kraken_db_build.slurm` | `kraken_db_build.py` → `db_v2` | current, only slurm file referencing `db_v2` |
| `benchmark_download.slurm`, `benchmark_workers.slurm` | `benchmark_download.py` | diagnostic, current |

Deleted (all pointed at dropped resources — `db_pathogens`, the removed `control/`
subpipeline, or a dead pre-restructure `runs.tsv` path): `kraken_classify.slurm`,
`kraken_classify_test.slurm`, `kraken_classify_pathogens_test.slurm`,
`smoke_classify.slurm`. A fresh wrapper will get written once `kraken_run_assign.py`
exists. `download_runs.slurm` moved to `kraken/legacy/` alongside `download.py`
(superseded by `kraken_run_select.py`'s built-in downloading).

### `kraken/pilot/` — legacy, excluded from git (kallisto + kraken pilot scripts, historical validation)
### `kraken/legacy/` — excluded from git (`download.py`, `download_runs.slurm` — superseded by kraken_run_select.py)

---

## Cross-module data flows (current, 2026-09-01)

```
stat/output/stat_build/data/phibase_db.json
  stat/stat_build.py  →→  stat/stat_fetch.py
                 →→  stat/stat_filter.py
                 →→  kraken/db/kraken_db_search.py
                 →→  metadata/meta_classify.py
                 →→  kraken/db/figures/compare_host_pathogen.py

stat/output/stat_fetch/data/stat_cache.jsonl  (unified: RunInfo + STAT per run)
  stat/stat_fetch.py  →→  stat/stat_filter.py
                      →→  stat/figures/prep_scatter.py

stat/output/stat_fetch/data/{mal,hal}_accessions.txt  (mode membership)
  stat/stat_fetch.py  →→  stat/stat_filter.py
                      →→  stat/figures/prep_scatter.py

stat/output/stat_filter/data/runs.tsv  — the single most-consumed file in the pipeline
  stat/stat_filter.py →→  metadata/meta_search.py
                      →→  metadata/meta_classify.py
                      →→  stat/figures/prep_scatter.py
                      →→  kraken/db/kraken_db_search.py
                      →→  kraken/run/kraken_run_select.py
                      →→  kraken/run/classify.py  (via --runs-tsv)
                      →→  kraken/db/figures/host_breakdown.py

metadata/output/meta_search/data/bioprojects.json
  meta_search.py  →→  meta_text.py
                  →→  meta_classify.py
                  →→  figures/prep_lit_resolution.py

metadata/output/meta_search/data/biosamples.json
  meta_search.py  →→  meta_classify.py

metadata/output/meta_text/data/text_cache.jsonl
  meta_text.py  →→  meta_classify.py (via bioprojects.json.full_text after --apply)
                →→  figures/prep_lit_resolution.py

metadata/output/meta_classify/data/samples.tsv  — PRIMARY ANALYSIS INPUT
  meta_classify.py  →→  figures/sample_funnel_v3.py
                    →→  kraken/run/kraken_run_select.py

kraken/output/db/search/data/ref_candidates.tsv
  kraken_db_search.py  →→  kraken_db_busco.py

kraken/output/db/busco/data/busco_scores.tsv
  kraken_db_busco.py  →→  kraken_db_build.py
                      →→  kraken/db/figures/busco_completeness.R

kraken/output/db/build/data/db_v2/
  kraken_db_build.py  →→  kraken/run/classify.py  (via --db, no code coupling)

kraken/output/run/select/data/run_list.tsv
  kraken_run_select.py  →→  (nothing yet — split/assign not built)
```

---

## Issues found (2026-09-01 survey)

### Resolved 2026-09-01 (Leon reviewed each item individually)
- `kraken/db/figures/compare_host_pathogen.py`, `kraken/db/figures/host_breakdown.py`,
  `kraken/db/figures/busco_completeness.R` — hardcoded pre-rename paths fixed.
- `metadata/tools/export_review_lists.py`, `metadata/tools/review_designs.py` —
  deleted (broken beyond docstring drift, hardcoded a dead `output/`/`scripts/` tree,
  would have needed a full rewrite).
- `kraken/slurm/kraken_classify.slurm`, `kraken_classify_test.slurm`,
  `kraken_classify_pathogens_test.slurm`, `smoke_classify.slurm` — deleted (all
  pointed at dropped resources; a fresh wrapper gets written once
  `kraken_run_assign.py` exists).
- `kraken/run/kraken_run_select.py`'s `reads/` output dir — `.gitignore` gap closed.
- `kraken/run/classify.py`'s dead `IN_DIR` fallback landmine — removed; now requires
  `--run-list`/`--runs-tsv` explicitly and fails loudly otherwise. Kept in the active
  tree (still the planned basis for `kraken_run_assign.py`).
- `kraken/legacy/download.py` — confirmed genuinely superseded by `kraken_run_select.py`'s
  built-in downloading; moved to `kraken/legacy/` along with its SLURM wrapper
  (`download_runs.slurm`).

### Also resolved 2026-09-01
- `metadata/classify_metadata.py`, `metadata/figures/sample_funnel_v2.py`,
  `metadata/figures/lit_resolution_sankey.py` (+ their tracked output figures) moved
  to `metadata/legacy/` — all three superseded (by `meta_classify.py`,
  `sample_funnel_v3.py`, `prep_lit_resolution.py`+`lit_resolution_alluvial.R`
  respectively). `sample_funnel.py` (v1) had already been removed before this survey.
- `download_strategy.md` (repo root) — deleted. Built on `sample_funnel_v2.py`'s
  classification, the dropped control-set, `download.py`, `db_pathogens`, and a
  500k-read subsampling strategy superseded by `kraken_run_select.py`'s full-file
  prefetch/fasterq-dump approach — every reference in it was dead.
- `stat/diagnostics/` (4 scripts) and `stat/migrate_cache.py` — deleted (one-off
  investigation/migration tools, already served their purpose; the latter's own
  docstring flagged its `stat_fetch.py` fallback as removable once done, so that
  fallback function was removed too).
- `metadata/tools/` — now gone entirely: the 2 broken scripts were deleted (see
  above), and `review_lists/` (8 .docx outputs, orphaned once their generator was
  deleted) moved to `metadata/legacy/tools/review_lists/`.

### CLAUDE.md drift beyond the above
- CLAUDE.md's "Module structure" section for metadata/ still describes the dropped
  four-script pipeline (`ncbi_metadata.py` etc.) as current, while its own "Output
  schemas" and "Active task list" sections further down already correctly reference
  `meta_classify.py`/`samples.tsv` — internally inconsistent, not just outdated.
- CLAUDE.md's active task #1 (`kraken_select.py`, marked ACTIVE/to-build) is done —
  `kraken_run_select.py` exists, works, and isn't mentioned in `kraken/README.md`'s
  module layout diagram at all.
- `stat/stat_build.py`'s own docstring lists stale consumer filenames
  (`stat/fetch_runs.py`, `stat/filter_runs.py`, `kraken/screen_refs.py`,
  `metadata/filter_kw.py`) — all four renamed or moved to legacy.
- The 2026-08-31 code-not-data gitignore policy excludes `phibase_db.json`,
  `runs.tsv`, `bioprojects.json`, `biosamples.json`, `samples.tsv`, `text_cache.jsonl`,
  `scatter_data.tsv` from git — CLAUDE.md's "Output / cache layout" section presents
  all of these without that caveat.

### Confirmed gone (not stale — actually removed)
- `kraken/control/` (whole control-set subpipeline: `sample.py`, `insilico.py`,
  `benchmark_download.py`'s control-set default, associated slurm scripts) — dropped
  2026-08-28, confirmed absent.
- The `metadata/figures/mal_guilds.*`, `field_hc_guilds.*`, `crypt_host_tree.*`,
  `coinf_rate.R`, `kingdom_comp.R`, `novel_heatmap.R` guild/network figure family from
  the previous draft — confirmed absent; these were explicitly deferred pending Kraken2
  validation (CLAUDE.md "guild figures deferred" note), not yet rebuilt against the
  current `samples.tsv` schema.

---

## Suggestions for Figma diagram

1. Three vertical swim lanes: **stat/** | **metadata/** | **kraken/** (kraken/ now
   splits into two visually distinct sub-lanes: DB-build submodule 1, complete;
   run/ submodule 2, only step 1 of 3 built)
2. Within stat/: three sequential boxes, stat_build → stat_fetch → stat_filter
3. `phibase_db.json` and `runs.tsv` are the two cross-module arteries — bold
   horizontal arrows crossing lanes, both explicitly marked "gitignored — not in git"
4. Within metadata/: linear chain meta_search → meta_text → meta_classify,
   `samples.tsv` labeled as the primary analysis input feeding the sankey figure and
   kraken/'s run/ submodule
5. Within kraken/: DB-build sub-lane (kraken_db_search → kraken_db_busco →
   kraken_db_build, ending at `db_v2`) parallel to run/ sub-lane (kraken_run_select
   done; split/assign shown as dashed/greyed "not built yet" boxes)
6. Figures collapse into one "Figures" box per module; mark the remaining
   superseded-but-live metadata figure scripts (`sample_funnel.py`/`_v2.py`,
   `lit_resolution_sankey.py`) with a distinct "stale" style rather than omitting
   them (useful for a future cleanup pass — see "Still open" above)
7. Slurm wrappers as a dashed cluster boundary around their target script (now just
   `kraken_db_search.slurm` → `kraken_db_busco.slurm` → `kraken_db_build.slurm` in
   sequence, plus the two benchmark wrappers — the stale classify variants are gone)
8. Grey/out-of-scope boxes: `stat/legacy/`, `metadata/legacy/`, `kraken/legacy/`,
   `kraken/pilot/` (`stat/diagnostics/` deleted 2026-09-01, no longer needs a box)
