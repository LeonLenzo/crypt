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
- **Feeds into:** `stat/stat_fetch.py`, `stat/stat_filter.py`, `kraken/kraken_db_search.py`, `metadata/meta_classify.py`, `kraken/figures/compare_host_pathogen.py` (path broken there, see Issues), `stat/diagnostics/*`
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
- **Feeds into:** `stat/stat_filter.py`, `stat/figures/prep_scatter.py`, `kraken/figures/compare_host_pathogen.py`

---

### `stat/stat_filter.py`
- **Inputs:** `--mode {mal|hal|both}`, `--skip-validate`; `stat/output/stat_build/data/phibase_db.json`; `stat/output/stat_fetch/data/stat_cache.jsonl`; `{mode}_accessions.txt`
- **Outputs:**
  - `stat/output/stat_filter/data/runs.tsv` — 10,995 rows (MAL 5,223 + HAL 5,772 gate-pass); key cols per CLAUDE.md's Output schemas section (gitignored — code-not-data policy)
  - `{mode}_kingdom_dist.tsv`, `{mode}_species_dist.tsv` (unless `--skip-validate`)
  - `stat/output/stat_filter/logs/history/<timestamp>/...` (now tracked)
- **Calls:** none (pure Python over already-local data)
- **Key algorithm:** `specific_hits()` — leaf-level species detection via count nesting; `KINGDOM_THRESHOLDS`: Fungi ≥0.5%, Oomycota ≥0.5%, Nematoda ≥1.0%; MAL gate: Viridiplantae ≥1%; HAL gate: any euk pathogen ≥1%
- **Feeds into:** the single most-consumed file in the pipeline — `stat/figures/prep_scatter.py`, `metadata/meta_search.py`, `metadata/meta_classify.py`, `kraken/kraken_db_search.py`, `kraken/kraken_run_select.py`, `kraken/classify.py`/`kraken/download.py` (via `--runs-tsv`)

---

### `stat/figures/prep_scatter.py`
- **Inputs:** `stat/output/stat_fetch/data/stat_cache.jsonl`, `stat/output/stat_filter/data/runs.tsv`, `stat/output/stat_fetch/data/{mal,hal}_accessions.txt` — genuinely 3 distinct sources, not reducible to `runs.tsv` alone: `runs.tsv` only holds gate-*pass* runs, so the plot's background/grey (gate-*fail*) points require the raw `stat_cache` + accession files to know each non-gate-pass run's mode. Even for gate-pass rows, host_pct/euk_pct are deliberately recomputed fresh from `stat_cache` rather than trusted from `runs.tsv`, so the gate lines are guaranteed consistent between foreground and background.
- **Outputs:** `stat/output/figures/scatter/scatter_data.tsv` (gitignored; layer/mode/coinf/host_pct/euk_pct)
- **Feeds into:** `stat/figures/scatter.R`

### `stat/figures/scatter.R`
- **Inputs:** `stat/output/figures/scatter/scatter_data.tsv`
- **Outputs:** `stat/output/figures/scatter/{scatter.pdf,scatter.png,scatter_caption.txt}` (tracked)
- **Calls:** ggplot2, dplyr, patchwork

### `stat/diagnostics/` — excluded from detail, still present (`diag_mal_gate.py`, `benchmark_strategies.py`, `check_primary_alignment.py`, `check_refseq_coverage.py`), not part of the main pipeline

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
- **Feeds into:** `metadata/figures/sample_funnel_v3.py`, `kraken/kraken_run_select.py`
- **Note (2026-09-01):** `host_disambig_cache.jsonl` records now store the exact `named_hosts` candidate list each `host_resolved` value was picked from — a cached resolution only counts as done if the current candidates still match (fixed a real bug where an improved hostpath prompt's better scientific-name candidates silently went unused for already-resolved BioSamples: 46% of 422 previously-attempted entries were stale on inspection).

---

### `metadata/classify_metadata.py` — superseded, not yet moved to `legacy/`
Older 3-call classifier (stress/setting/tissue only, no full-text gate, no host/pathogen
extraction). Output `metadata/output/classify_metadata/data/samples.tsv` is gitignored
specifically to avoid it looking like a live duplicate of `meta_classify.py`'s current
output. Only consumer left: `metadata/figures/sample_funnel_v2.py` (also superseded,
see below). Low-priority cleanup: move both to `legacy/`.

### `metadata/figures/sample_funnel_v3.py` — current primary figure
- **Inputs:** `metadata/output/meta_classify/data/samples.tsv`
- **Outputs:** `metadata/output/figures/sankey/sample_funnel_v3.{html,svg,png}`
- **Calls:** plotly (+kaleido)
- **Supersedes:** `sample_funnel.py`, `sample_funnel_v2.py` (both stale, read from `classify_metadata.py`'s dead-end output — kept for reference, not moved to legacy)

### `metadata/figures/prep_lit_resolution.py` + `metadata/figures/lit_resolution_alluvial.R` — current literature-resolution figure
- **Inputs:** `bioprojects.json`, `text_cache.jsonl` → `metadata/output/figures/sankey/lit_resolution_data.tsv` → `lit_resolution_alluvial.{png,svg}` (ggplot2 + ggalluvial)
- **Supersedes:** `metadata/figures/lit_resolution_sankey.py` (Plotly) — explicitly called "retired" in `prep_lit_resolution.py`'s own docstring, though never moved to legacy and still present at top level. **CLAUDE.md's active task list still lists building the Plotly version as an open TODO — it's stale; the R version already superseded it.**

### `metadata/tools/export_review_lists.py`, `metadata/tools/review_designs.py` — broken, not just stale
Both hardcode a pre-restructure `output/`/`scripts/` path tree (e.g.
`output/05_llm_classify/data/bioproject_llm.tsv`) that **no longer exists anywhere in
the repo**. Neither can run without a rewrite — this isn't docstring drift, it's a
hard `FileNotFoundError`.

### `metadata/legacy/` — excluded from git entirely, confirmed present locally (ncbi_metadata.py, web_metadata.py, undermind.py, fetch_xml.py, fetch_lit.py, filter_kw.py, llm_classify.py, serper_dump/resolve/scrape.py, backfill_doi.py, tissue_vocab.py, + legacy figures)

---

## Module: kraken/

**Two submodules.** Submodule 1 (DB build) is fully restructured and complete.
Submodule 2 (run/select/split/assign) exists only as step 1 of 3 — `kraken_run_select.py`
is built and (per 2026-09-01 testing) runs, but `kraken_run_split.py`/`kraken_run_assign.py`
do not exist yet. `kraken/control/` (the old stratified control-set subpipeline) is
confirmed **fully removed** — dropped 2026-08-28, not just deprecated.

---

### `kraken/kraken_db_search.py` — submodule 1, step 1/3
- **Inputs:** `stat/output/stat_build/data/phibase_db.json`, `stat/output/stat_filter/data/runs.tsv`; NCBI `datasets` CLI
- **Outputs:** `kraken/output/kraken_db_search/data/{ref_candidates.tsv, pangenome.tsv, genus_fill.tsv (--scope only)}` (tracked) + `cds/pathogen/{accession}/` (gitignored, `--download` only)
- **Calls:** `datasets summary/download` (subprocess), `unzip`
- **Feeds into:** `kraken/kraken_db_busco.py`, `kraken/kraken_db_build.py`, `kraken/kraken_run_select.py` (imports functions directly)

### `kraken/kraken_db_busco.py` — submodule 1, step 2/3
- **Inputs:** `ref_candidates.tsv`; CDS already on disk at `cds/pathogen/`; BUSCO Singularity + lineage DBs (Setonix)
- **Outputs:** `kraken/output/kraken_db_busco/data/{busco_scores.tsv, busco_scan_cache.tsv}` — thresholds fungi ≥50%, oomycete ≥65% (per-taxid fallback)
- **Calls:** BUSCO via `subprocess.Popen`
- **Feeds into:** `kraken/kraken_db_build.py`, `kraken/figures/busco_completeness.R` (path broken there, see Issues)

### `kraken/kraken_db_build.py` — submodule 1, step 3/3
- **Inputs:** `busco_scores.tsv` (rows `selected=True`); CDS from `cds/pathogen/`; NCBI taxdump (auto-downloaded)
- **Outputs:** Kraken2 DB directory — current production `db_v2` (20GB, gitignored, `kraken/output/kraken_db_build/data/db_v2/`); optional `--upload-to-acacia`
- **Calls:** `kraken2-build --add-to-library` / `--build`; `_util.upload_to_acacia` → `aws s3 sync`
- **Feeds into:** `kraken/classify.py` (via `--db` at run time, no code coupling)

---

### `kraken/kraken_run_select.py` — submodule 2, step 1/3 (only step built)
- **Inputs:** `metadata/output/meta_classify/data/samples.tsv`, `stat/output/stat_filter/data/runs.tsv`; imports from `kraken_db_search.py`; requires `prefetch`/`fasterq-dump` (sra-tools) + `datasets` CLI
- **Outputs:** `kraken/output/kraken_run_select/data/run_list.tsv` (tracked), `.../data/reads/{run}_{1,2}.fastq.gz` (**gitignore gap — not yet excluded, unlike every other kraken large-data subdir; fix before running at scale**), `kraken/output/kraken_db_search/data/cds/host/{accession}/` (gitignored, shared pool with pathogen CDS)
- **CLI:** `--setting`, `--aerial-only`/`--no-aerial-only`, `--cryptic-only`, `--limit N`, `--download`, `--workers N`
- **Feeds into:** nothing yet wired downstream — `kraken_run_split.py`/`kraken_run_assign.py` not built
- **Status (2026-09-01):** smoke-tested with `--limit 2 --download` on Setonix; host resolution worked, but the host genome fetch step hit a 300s timeout — leading theory is it ran on the Setonix login node rather than through SLURM (untested at time of writing, Setonix down for maintenance)

### `kraken/kraken_run_split.py`, `kraken/kraken_run_assign.py` — **not built**
Design agreed (BBSplit host-read removal via a combined multi-reference index;
Kraken2 classification adapted from `classify.py`) but zero code written.

---

### `kraken/classify.py`
- **Inputs:** `--mode {mal|hal|both}` / `--runs-tsv` / `--run-list`; `--reads-dir` (pre-downloaded) or streams from ENA FTP; `--db PATH` (required)
- **Outputs:** `kraken/output/classify/data/{kraken_cache.jsonl, kraken_cache_index.txt}` (gitignored), optional per-run `.kreport` via `--reports-dir`
- **Calls:** `kraken2` (subprocess, confidence=0.15, min-hit-groups=3)
- **Feeds into:** nothing currently — no figure script reads `kraken_cache.jsonl` yet; `kraken/figures/compare_host_pathogen.py` instead expects a manually-scp'd `/tmp/setonix_kraken_metrics.json`
- **Note:** dead fallback constant `IN_DIR = stat/output/fetch_runs/data` (pre-rename path) only triggers if neither `--run-list` nor `--runs-tsv` is passed — a landmine, not a hard break, for anyone running the bare documented invocation

### `kraken/download.py`
- **Inputs:** `--run-list` or `--runs-tsv` + `--biosample-rep`/`--hc`; `--reads-dir` (required)
- **Outputs:** `{reads-dir}/{run}_{1,2}.fastq.gz`; `download_index.txt` + `errors.log`
- **Calls:** curl (subprocess, ENA HTTPS)
- **Feeds into:** `kraken/classify.py --reads-dir`
- **Note:** docstring usage examples still reference `kraken/control/output/data/run_ids.txt` — dead, control subpipeline dropped

### `kraken/benchmark_download.py`
- **Inputs:** `--methods {ena-ftp,ena-https,sra-s3,prefetch}`, `--run-list` (else falls back to dead `kraken/control/output/data/control_runs.tsv` constant)
- **Outputs:** results TSV per `--out`
- **Purpose:** throughput benchmarking, diagnostic only

---

### `kraken/figures/` — all three currently broken (hardcoded stale paths, not just stale docstrings)

| Script | Broken input path | Correct path |
|---|---|---|
| `compare_host_pathogen.py` | `stat/output/build/data/phibase_db.json` | `stat/output/stat_build/data/phibase_db.json` |
| `host_breakdown.py` | `stat/output/filter_runs/data/runs.tsv` | `stat/output/stat_filter/data/runs.tsv` |
| `busco_completeness.R` | `kraken/busco_scores.tsv` (relative, doesn't exist) | `kraken/output/kraken_db_busco/data/busco_scores.tsv` |

`compare_host_pathogen.py` also depends on `/tmp/setonix_kraken_metrics.json`, a
manually-scp'd file not produced by any tracked script.

---

### `kraken/slurm/` — several stale, some pointing at dropped resources

| Script | Calls | Status |
|---|---|---|
| `kraken_db_search.slurm` | `kraken_db_search.py --download --workers 16` | current (16 CPU/32G/8h) |
| `kraken_db_busco.slurm` | `kraken_db_busco.py` | current (128 CPU/230G/24h) |
| `kraken_db_build.slurm` | `kraken_db_build.py` → `db_v2` | current, only slurm file referencing `db_v2` |
| `kraken_classify.slurm` | `classify.py --runs-tsv output/02_filter_runs/...` | **stale** — dead pre-restructure path + generic (non-`db_v2`) DB dir |
| `kraken_classify_test.slurm` | `classify.py` → `$SCRATCH/kraken/db` | **stale** — not `db_v2` |
| `kraken_classify_pathogens_test.slurm` | `classify.py` → `db_pathogens` | **stale** — `db_pathogens` dropped 2026-08-28 |
| `smoke_classify.slurm` | `classify.py --run-list kraken/control/...` | **stale on two counts** — dropped control subpipeline + `db_pathogens` |
| `download_runs.slurm` | `download.py` | current; comment references a not-yet-existing `host_split.py (TBD)` |
| `benchmark_download.slurm`, `benchmark_workers.slurm` | `benchmark_download.py` | diagnostic, current |

### `kraken/pilot/` — legacy, excluded from git (kallisto + kraken pilot scripts, historical validation)

---

## Cross-module data flows (current, 2026-09-01)

```
stat/output/stat_build/data/phibase_db.json
  stat/stat_build.py  →→  stat/stat_fetch.py
                 →→  stat/stat_filter.py
                 →→  kraken/kraken_db_search.py
                 →→  metadata/meta_classify.py
                 →→  kraken/figures/compare_host_pathogen.py  (broken path — see Issues)

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
                      →→  kraken/kraken_db_search.py
                      →→  kraken/kraken_run_select.py
                      →→  kraken/classify.py / kraken/download.py  (via --runs-tsv)
                      →→  kraken/figures/host_breakdown.py  (broken path)

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
                    →→  kraken/kraken_run_select.py

kraken/output/kraken_db_search/data/ref_candidates.tsv
  kraken_db_search.py  →→  kraken_db_busco.py

kraken/output/kraken_db_busco/data/busco_scores.tsv
  kraken_db_busco.py  →→  kraken_db_build.py
                      →→  kraken/figures/busco_completeness.R  (broken path)

kraken/output/kraken_db_build/data/db_v2/
  kraken_db_build.py  →→  kraken/classify.py  (via --db, no code coupling)

kraken/output/kraken_run_select/data/run_list.tsv
  kraken_run_select.py  →→  (nothing yet — split/assign not built)
```

---

## Issues found (2026-09-01 survey)

### Genuinely broken, not just stale docstrings
- `kraken/figures/compare_host_pathogen.py`, `kraken/figures/host_breakdown.py`,
  `kraken/figures/busco_completeness.R` — hardcoded pre-rename paths in actual code,
  will error/produce nothing if run today (table above).
- `metadata/tools/export_review_lists.py`, `metadata/tools/review_designs.py` —
  hardcode a pre-restructure `output/`/`scripts/` tree that no longer exists anywhere
  in the repo. Not runnable without a rewrite.
- `kraken/slurm/kraken_classify.slurm`, `kraken_classify_test.slurm`,
  `kraken_classify_pathogens_test.slurm`, `smoke_classify.slurm` — point at dropped
  resources (`db_pathogens`, the removed `kraken/control/`, or a dead pre-restructure
  `runs.tsv` path). Only `kraken_db_build.slurm`/`kraken_db_busco.slurm` reference
  current `db_v2`.
- `kraken/kraken_run_select.py`'s `reads/` output dir has no `.gitignore` coverage —
  every other kraken large-data subdir does. Worth fixing before running at scale,
  given last session's near-misses on committing large data.

### Landmines (work today, break on the "obvious" invocation)
- `kraken/classify.py` and `kraken/download.py` both fall back to a dead
  `IN_DIR = stat/output/fetch_runs/data/{mode}_runs.json` path (format changed to
  `stat_cache.jsonl` at rename) if neither `--run-list` nor `--runs-tsv` is passed.

### Superseded but not moved to legacy/ (dead-end branches still runnable)
- `metadata/classify_metadata.py` → superseded by `meta_classify.py`
- `metadata/figures/sample_funnel.py`, `sample_funnel_v2.py` → superseded by `sample_funnel_v3.py`
- `metadata/figures/lit_resolution_sankey.py` (Plotly) → superseded by
  `prep_lit_resolution.py` + `lit_resolution_alluvial.R` (R/ggalluvial). **CLAUDE.md's
  active task list still lists building the Plotly version as an open TODO — stale,
  the R version already exists and superseded it before that TODO was ever removed.**

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
6. Figures collapse into one "Figures" box per module; mark the 3 broken
   `kraken/figures/*` scripts and the 2 superseded-but-live metadata figure scripts
   with a distinct "stale" style rather than omitting them (useful for a cleanup pass)
7. Slurm wrappers as a dashed cluster boundary around their target script; mark the
   4 stale kraken_classify* slurm variants distinctly
8. Grey/out-of-scope boxes: `stat/legacy/`, `metadata/legacy/`, `stat/diagnostics/`,
   `kraken/pilot/`, `metadata/tools/review_designs.py`,
   `metadata/tools/export_review_lists.py` (both broken)
