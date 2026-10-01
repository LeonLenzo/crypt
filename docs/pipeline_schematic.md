# Pipeline schematic

What each step consumes, what it produces, and what depends on it. Regenerated
2026-10-01 against the numbered layout.

**This document deliberately does not repeat each script's options or output columns.**
Those live in the script's own docstring, where they cannot drift out of step with the
code. The previous version of this file duplicated them and went stale across two
restructures. What belongs here is the chain: the order, the hand-offs, and where the
data physically sits.

Run every script from the repository root (`phd/01-review`). A numbered directory is a
pipeline step; an unnumbered one (`slurm/`, `utilities/`, `data/`, `logs/`, `figures/`)
is support. Each step keeps its own `data/` and `logs/`, so a step's outputs sit beside
the code that wrote them.

---

## The chain

```
01_stat/01_build      PHI-base CSV + NCBI taxonomy
                         └─> phibase_db.json ─────────────────────┐
                                                                  │
01_stat/02_fetch      SRA RunInfo + NCBI STAT profiles            │
                         └─> stat_cache.jsonl ───────┐            │
                                                     │            │
01_stat/03_filter     gate → detect → dedup  <───────┴────────────┤
                         └─> runs.tsv ──────────────┬─────────────┤
                                                    │             │
02_literature/01_search   BioProject/BioSample +    │             │
                          literature resolution <───┤             │
                         └─> bioprojects.json ──┐   │             │
                         └─> biosamples.json ───┤   │             │
                                                │   │             │
02_literature/02_text     full-text retrieval <─┤   │             │
                         └─> full_text written back into          │
                             bioprojects.json (--apply)           │
                                                │   │             │
02_literature/03_classify  LLM study design <───┴───┴─────────────┤
                         └─> samples.tsv ───────────┬─────────────┘
                                                    │
03_kraken/01_search      reference selection +      │
                         CDS download  <────────────┴── phibase_db.json, runs.tsv
                         └─> ref_candidates.tsv
                         └─> data/cds/{pathogen,host}/ ──┐
                                                         │
     utilities/busco      completeness QC (not a gate) <─┤
                         └─> busco_scores.tsv ───────┐   │
                                                     │   │
03_kraken/02_build       kraken2-build  <────────────┴───┘
                         └─> data/db_v3/ + db_v3_inspect.tsv ──┐
                                                               │
03_kraken/03_select      target BioSamples + read download     │
                         <── samples.tsv, runs.tsv             │
                         └─> run_list.tsv, data/reads/ ──┐     │
                                                         │     │
03_kraken/04_assign      kraken2 classify  <─────────────┴─────┤
                         └─> data/reports/*.txt ──┐            │
                         host/ ─> host_calls.tsv ─┤            │
                                                  │            │
03_kraken/05_filter      artefact score + coverage │            │
                         + co-occurrence  <───────┴────────────┘
                         └─> data/_classified.tsv ──┐
                                                    │
04_align                 competitive read assignment│
                         <──────────────────────────┘ + cds/pathogen/
                         └─> (specified, not yet run)
```

`01_stat/01_build/data/phibase_db.json` is the one artefact consumed by three different
modules (`02_fetch`, `03_filter`, `02_literature/03_classify`, `03_kraken/01_search`). It
is the pipeline's single reference anchor; nothing downstream resolves taxa independently.

---

## Step by step

Paths are relative to the repository root. "Large" marks data that is gitignored and
lives on Setonix scratch, represented in the repo by a `manifest.tsv`.

### 01_stat — screen the archive

| Step | Consumes | Produces |
|---|---|---|
| `01_build/stat_build.py` | PHI-base CSV (`data/phi-base_current.csv`), local NCBI taxonomy via `ete3` | `01_build/data/phibase_db.json` |
| `02_fetch/stat_fetch.py` | `phibase_db.json` | `02_fetch/data/stat_cache.jsonl` (large), `stat_cache_index.txt`, `{mode}_uid_checkpoint.txt`, `{mode}_accessions.txt` |
| `03_filter/stat_filter.py` | `phibase_db.json`, `stat_cache.jsonl`, `{mode}_accessions.txt` | `03_filter/data/runs.tsv` |

Two things here are load-bearing and easy to break:

- `_expand_taxids()` in `stat_build.py` must pass `intermediate_nodes=True`. The default
  returns only leaf taxa and silently drops species nodes that have named strains as
  children, which hits the rust formae speciales directly.
- `stat_fetch.py` writes a single shared `stat_cache.jsonl` for both modes, so MAL and
  HAL cannot run concurrently. Chain them.

`specific_hits()` in `stat_filter.py` is the non-trivial part: it finds leaf-level species
detections by locating k-mer counts not nested under any more specific count. Genus-level
STAT counts reflect LCA promotion of shared k-mers and are not usable for co-infection.

### 02_literature — establish study design

| Step | Consumes | Produces |
|---|---|---|
| `01_search/meta_search.py` | `01_stat/03_filter/data/runs.tsv` | `01_search/data/{bioprojects.json, biosamples.json}` + `bp_cache.json`, `bs_cache.json`, `serper_cache.json` |
| `02_text/meta_text.py` | `01_search/data/bioprojects.json` | `02_text/data/{text_cache.jsonl, failed_dois.tsv}`; `--apply` writes `full_text` back into `bioprojects.json` |
| `03_classify/meta_classify.py` | `phibase_db.json`, `runs.tsv`, `bioprojects.json`, `biosamples.json` | `03_classify/data/samples.tsv` + `classify_cache.jsonl`, `host_disambig_cache.jsonl` |

The entry gate is the module's defining constraint: **only BioProjects with retrievable
full text are classified at all.** `samples.tsv` is therefore already the full-text-gated
population, and no separate evidence stage is needed downstream. The predecessor
classified everything including bare titles, guessing from near-empty context.

Two deliberate divisions of labour:

- Author-named pathogens and hosts are resolved to NCBI taxids in plain Python
  (`_util.resolve_taxon_name()`), never asked of the LLM. An LLM recalling taxids from
  memory produces confident wrong numbers.
- The downstream fields (`pathogen_match_status`, the agreement flags,
  `metadata_disagreement_flag`) are computed in Python, not LLM judgments, so the logic
  can be revised without re-running paid calls.

`--disambiguate-hosts` is a separate per-BioSample pass, cached by BioSample rather than
BioProject, for the studies whose `named_hosts` holds more than one species (a rust's
alternate-host life cycle, for instance). Without it a BioProject-constant host would be
duplicated across every sample, which is wrong for per-sample analysis.

### 03_kraken — detect from the reads

| Step | Consumes | Produces |
|---|---|---|
| `01_search/kraken_search.py` | `phibase_db.json`, `runs.tsv`, holobase (`_refdata/holobase.db`) | `01_search/data/ref_candidates.tsv` (symlink to the versioned table), `data/cds/pathogen/` (large) |
| `01_search/kraken_hosts.py` | `03_select/data/run_list.tsv` | `01_search/data/host_candidates.tsv`, `data/cds/host/` (large) |
| `utilities/busco/kraken_busco.py` | `data/cds/` already on disk (never downloads) | `utilities/busco/data/busco_scores.tsv`, `busco_scan_cache.tsv` (both tracked) |
| `02_build/kraken_build.py` | `busco_scores.tsv`, `data/cds/{pathogen,host}/` | `02_build/data/db_v3/` (large), `db_v3_inspect.tsv` (tracked) |
| `03_select/kraken_select.py` | `02_literature/03_classify/data/samples.tsv`, `runs.tsv` | `03_select/data/run_list.tsv`, `data/reads/` (large) |
| `04_assign/kraken_assign.py` | `data/reads/`, `db_v3/` | `04_assign/data/reports/*.txt`, `kraken_cache.jsonl` |
| `04_assign/host/kraken_host_call.py` | `04_assign/data/reports/` | `04_assign/host/data/host_calls.tsv` |
| `05_filter/detect.py` | `04_assign/data/reports/`, `host_calls.tsv`, `db_v3_inspect.tsv` | `05_filter/data/_classified.tsv` |

Ordering notes that are not obvious from the tree:

- `kraken_hosts.py` reads `03_select/data/run_list.tsv`, so host CDS selection runs
  **after** read selection even though both sit under `01_search`. The host set is
  defined by which runs were actually chosen.
- `01_search` is the only step in the database chain that downloads. BUSCO and the build
  both read CDS already on disk.
- `busco_scan_cache.tsv` is tracked in git (1,082 accessions with full scores) and
  `run_scan()` builds its skip list from it, so a rebuild does not rescan.
- A threshold change needs no rescan: `run_finalize()` recomputes `busco_status` from
  `complete_pct`, so `--finalize-only` is enough.

BUSCO sits in `utilities/`, not in the numbered chain, because it is a QC record rather
than a filter. Of 2,063 candidates, 1,986 pass, 7 fall below the completeness bar and 70
have no CDS; a per-taxid fallback then adds 34 that failed but are their taxon's only
representative. The screen cannot remove a taxon, only prune surplus assemblies within
one. The attrition worth watching is the 70 missing annotations, not the 7 failures.

`04_assign` requires `--reads-dir`; the ENA streaming path was removed once every read
was local. Report read **counts**, not percentages: a percentage moves whenever the
dominant organism's abundance moves, which conflates two different things.

### 04_align — resolve what k-mers cannot

| Step | Consumes | Produces |
|---|---|---|
| `align.py` | `03_kraken/05_filter/data/_classified.tsv`, `03_kraken/01_search/data/cds/pathogen/` | specified, not yet run |

Promoted out of `03_kraken` because it answers a different question by a different method
(kallisto competitive EM, not k-mer counts), and because its output calibrates the
`05_filter` thresholds. Nesting it under kraken inverted that dependency. See
[../04_align/README.md](../04_align/README.md).

---

## Figures, and which step owns them

A figure script lives with the step that produced the newest data it reads.

| Figure script | Reads | Lives in |
|---|---|---|
| `prep_scatter.py` + `scatter.R` | `runs.tsv` | `01_stat/03_filter/figures/` |
| `prep_lit_resolution.py` + `lit_resolution_alluvial.R` | `text_cache.jsonl` | `02_literature/02_text/figures/` |
| `sample_funnel_v3.py` | `samples.tsv` | `02_literature/03_classify/figures/` |
| `crypt_host_tree.py` + `.R` | `samples.tsv` | `02_literature/03_classify/figures/` |
| `host_breakdown.py` | `runs.tsv` | `03_kraken/03_select/figures/` |
| `prep_taxonomy.py`, `composite.R`, `prep_network.py`, `network.R`, `prep_slopes.py`, `slopes.R`, `prep_signal_split.py`, `signal_split.R`, `model2d.R`, `example_detections.R`, `prep_method_comparison.py`, `method_comparison.R` | `_classified.tsv` | `03_kraken/05_filter/figures/` |
| `busco_completeness.R` | `busco_scores.tsv` | `03_kraken/utilities/busco/figures/` |
| `saturation_curves.py` | saturation sweep output | `03_kraken/utilities/saturation/figures/` |

`host_breakdown.py` is the one awkward case: it is a kraken-side figure built on STAT
output. It sits with `03_select` because that is the kraken step whose cohort it describes.

Every figure follows the same two-stage convention: a `prep_*.py` that writes a plotting
table, then an `.R` that renders it. The split exists so the expensive data assembly is
not repeated on every styling change.

---

## Shared root

| File | Provides |
|---|---|
| `_util.py` | `_Tee`, `make_log_dir`, `link_latest`, `http_get`, `load_json`, `save_json`, `build_manifest`, `resolve_taxon_name` (+ `HOST_NAME_ALIASES`), `upload_to_acacia` |
| `manifest.py` | CLI over `build_manifest()`: run against any gitignored bulk-data directory so the repo records what exists without holding it |

Scripts reach `_util` by inserting the repository root on `sys.path` with a
`Path(__file__).resolve().parents[N]` walk, never by package import. This is why numbered
directory names are safe: a leading digit would break `import`, but nothing imports by
module name.

---

## Where the data actually lives

Three locations, by design:

- **In the repo.** Code, small tracked TSVs, figures, `manifest.tsv` files, and log
  history (`logs/history/`, `*.latest`). Logs are tracked deliberately, for transparency.
- **Setonix scratch.** Everything marked large above: CDS, reads, the Kraken2 database,
  the per-run reports. Subject to the 21-day purge, which is why `04_assign` output (66 MB
  against 14 TB of input reads) is copied off promptly.
- **Acacia.** `s3://pawsey1168-kraken-db` on the project allocation, not the personal
  one. The `--endpoint-url` flag is mandatory; omitting it returns a misleading
  `InvalidAccessKeyId`.

No research data is committed. Verify with `git check-ignore` before staging, not after.
