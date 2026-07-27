# crypt

Screens NCBI SRA plant RNA-seq runs for cryptic co-infections using NCBI STAT pre-computed k-mer taxonomy — no read re-alignment required.

## Background

Plant disease studies deposited in SRA are designed around a single target pathogen. Field-collected samples, however, routinely harbour additional co-infecting organisms that go undetected and unreported under single-target study designs. We hypothesise that a substantial fraction of publicly available plant RNA-seq data contains secondary pathogen signal sufficient for detection via k-mer taxonomy, representing a largely untapped resource for co-infection epidemiology.

NCBI STAT provides pre-computed 32-mer taxonomy profiles for all SRA runs. By cross-referencing STAT outputs against PHI-base (plant–pathogen interactions) and the ICTV plant virus master species list, secondary pathogens can be identified in runs where they were not the study target — without downloading or re-aligning raw data.

## Two screening modes

| Mode | Abbrev | Library organism | Retention gate | Signal |
|------|--------|-----------------|----------------|--------|
| Microbe-as-library | MAL | PHI-base plant pathogen | Viridiplantae ≥ 1% of STAT reads | Secondary pathogens co-infecting the host plant |
| Host-as-library | HAL | PHI-base plant host species | Any PHI-base pathogen ≥ 1% of STAT reads | Primary and secondary pathogens in host transcriptomes |

MAL and HAL are complementary and largely non-overlapping. MAL targets pathogen-focused sequencing (population genomics, disease surveys, pathotype characterisation); HAL targets host-focused sequencing (resistance transcriptomics, field transcriptomics).

## Dependencies

```
python >= 3.11      stdlib only for steps 01–03; miniconda compatible
system python3      required for 00_build.py (ete3/sqlite3 ABI conflict with miniconda)
ete3                NCBI taxonomy resolution (00_build.py only)
openpyxl            ICTV VMR Excel parsing (00_build.py only)

R (figure scripts):
  ggraph, igraph, ggplot2, dplyr, tibble, ggtree, ape, cowplot
```

Set `NCBI_API_KEY` in `~/.bashrc` for 10 req/s Entrez access (2.5 req/s without).
Set `S2_API_KEY` for Semantic Scholar access in 03_find.py (skipped without key).

## Quick start

```bash
# Build PHI-base + ICTV reference DB (system python3 required)
python3 00_build.py

# Fetch SRA run IDs + STAT taxonomy profiles
# Run MAL before HAL — shared stat cache; use tmux for long runs
python 01_fetch.py --mode mal
python 01_fetch.py --mode hal

# Apply retention gate, calibrate thresholds, classify co-infections
python 02_filter.py               # both modes, unified output
python 02_filter.py --skip-validate   # skip threshold tables after review

# Fetch BioProject metadata and link publications
python 03_find.py
python 03_find.py --hc            # restrict to same_genus_secondary=False
```

Steps 01–03 are resumable via their respective caches. Long-running STAT fetches should be run in tmux.

## Pipeline

### `00_build.py` — reference database

Builds `output/00_build/data/phibase_db.json`, the reference database used by all downstream steps. Sources PHI-base (fungi, bacteria, oomycetes, nematodes) and the ICTV Virus Metadata Resource (plant viruses), downloading both automatically if absent. Uses `ete3 NCBITaxa` to resolve names to NCBI taxids and expand each seed species to all descendant strains, subspecies, and formae speciales. The expansion uses `intermediate_nodes=True`, which is critical: many plant pathogen species (e.g. *Potato virus Y*, rust f. sp. taxa) sit as internal nodes in the NCBI taxonomy with named strains as children — without this flag they would be silently excluded.

The output DB stores pathogen taxids grouped by kingdom (fungi/bacteria/oomycetes/nematodes/viruses), host taxids, and a bidirectional host–pathogen interaction map derived from PHI-base. Must be run with system `python3` due to an ete3/miniconda sqlite3 ABI incompatibility.

### `01_fetch.py` — SRA run fetch + STAT taxonomy

Fetches SRA run accessions and NCBI STAT k-mer taxonomy profiles for either MAL or HAL mode. No retention gate is applied here — all runs are fetched and screened regardless.

**SRA fetch**: queries NCBI Entrez using `txid{taxid}[Organism:exp]` (library source organism, not the `[Host]` field, which is free text and returns nothing useful). For MAL, this means 205 PHI-base pathogen seed taxids; for HAL, 180 plant host seed taxids. `[Organism:exp]` expands to all descendant strains automatically. Run metadata (platform, BioProject, BioSample, scientific name) is fetched in parallel batches via POST to avoid HTTP 414 errors.

**STAT fetch**: queries `trace.ncbi.nlm.nih.gov` for each run's pre-computed k-mer taxonomy profile. Results are written to a shared append-only cache (`stat_cache.jsonl`) so MAL and HAL share data without duplication. A `stat_cache_index.txt` file (accessions only) is read at startup for fast resume. Do not run MAL and HAL simultaneously — the cache is shared and only one writer is safe.

### `02_filter.py` — gate, validate, and classify

Three sequential phases operating on the STAT cache:

**Phase 1 — retention gate**: keeps runs where the primary host signal confirms the sample is biologically relevant. MAL requires ≥ 1% Viridiplantae reads (confirming in planta context); HAL requires ≥ 1% of any known PHI-base pathogen or plant virus (confirming pathogen presence in a host transcriptome).

**Phase 2 — validate**: generates per-kingdom read percentage distributions as TSVs and prints breakpoint tables and ASCII histograms. Used to review and adjust detection thresholds (`KINGDOM_THRESHOLDS`) before committing to Phase 3. Skip with `--skip-validate` once thresholds are set.

**Phase 3 — crypt**: the core co-infection detection step. For each retained run, `specific_hits()` identifies leaf-level species detections in STAT by finding counts that are not nested under any other count (i.e., genuinely species-diagnostic signal, not genus-level aggregates). Detected species are cross-referenced against the PHI-base/ICTV DB. Each secondary pathogen is classified as `known` (interaction in PHI-base), `novel_host_range` (pathogen known but not on this host), `novel_combination` (both organisms known but interaction not recorded), or `unresolved` (taxid lookup failed). BioSamples with multiple runs are deduplicated and a `biosample_representative` flag marks the single highest-coverage run per sample.

Output is a unified `crypt.tsv` with a `mode` column covering both MAL and HAL.

### `03_find.py` — BioProject metadata and publication linking

Fetches title, submission date, study design, and linked publications for all 1,797 BioProjects in `crypt.tsv` — including single-infection projects, which serve as the denominator for field co-infection prevalence estimates.

For each BioProject, PMID search uses a short-circuit strategy: six sources are tried in order and the search stops at the first hit. The order was determined empirically by benchmarking discovery yield on unlabelled entries (BioProjects with no known PMID):

| Strategy | Yield | Notes |
|----------|-------|-------|
| BioProject XML | 15% | Free — extracted from the title fetch XML |
| PMC full-text | 40% | Highest yield; searches full text including methods |
| Europe PMC | 35% | Catches preprints and supplementary-method mentions |
| ENA XML | — | Instant; only useful for PRJEB accessions |
| Semantic Scholar | — | Requires `S2_API_KEY`; skipped without key |
| PubMed | — | Last resort; lowest yield |

`study_design` is inferred from title, description, and abstract keywords: `coinf_experiment` (intentional mixed-infection design — exclude from novel interaction counts), `field_survey` (field-collected samples), or `unclear`.

Results are cached in `meta_cache.json` (versioned; v4). Entries with no PMID are retried on each run in case new publications have appeared. The output `bioproject_meta.tsv` includes `n_coinf`, `n_single`, and `coinf_rate` per BioProject for downstream prevalence analysis.

## Output

`output/02_filter/data/crypt.tsv` — one row per confirmed run:

| Column | Description |
|--------|-------------|
| `mode` | `mal` or `hal` |
| `host` | Top Viridiplantae leaf in STAT (MAL) or library organism (HAL) |
| `primary_pathogen` | Library organism (MAL) or top STAT-detected pathogen (HAL) |
| `secondary_pathogens` | Additional PHI-base/ICTV pathogens detected; semicolon-separated |
| `co_infection_flag` | `single` / `multi_species` / `multi_kingdom` |
| `same_genus_secondary` | `True` if any secondary shares genus with primary |
| `interaction_status` | `known` / `novel_host_range` / `novel_combination` / `unresolved` |
| `biosample_representative` | `True` for one run per BioSample (use for sample-level statistics) |

`novel_host_range` — pathogen confirmed in PHI-base but not recorded on this host species — is the primary signal of interest.

Filter to `co_infection_flag != "single"` for co-infected runs; additionally `same_genus_secondary == "False"` for high-confidence detections; `biosample_representative == "True"` for sample-level statistics.

`output/03_find/data/bioproject_meta.tsv` — one row per BioProject, including `study_design` (`coinf_experiment` / `field_survey` / `unclear`) inferred from BioProject title, description, and linked publication abstracts. `coinf_experiment` projects require manual methods review before inclusion in novel interaction counts.

## Results (MAL + HAL complete — 2026-07)

**Combined** (607,746 runs screened):
- 14,335 confirmed runs (11,853 BioSample-representative); 1,797 BioProjects
- co_infection_flag (biosample_representative): single 10,078 | known co-infection 1,196 | novel_host_range 579
- 436 BioProjects with co-infected runs; 1,361 single-infection BioProjects

**MAL** (48,418 runs): 6,852 confirmed (14.2% gate pass); top secondaries: *Puccinia graminis*, *Zymoseptoria tritici*, *Pepino mosaic virus*, *Botrytis cinerea*

**HAL** (559,328 runs): ~7,483 confirmed (~1.3% gate pass); complementary biology — legume viruses, apple and grapevine pathogens, corn and beet pathogens

## Interpreting `same_genus_secondary`

STAT uses an LCA-based k-mer merging strategy: k-mers shared between sibling species are promoted to the genus node at database build time and are not assigned to either species. Species-level detections therefore represent genuinely species-diagnostic signal — inter-species k-mer bleed within a genus is largely prevented by design. Same-genus secondaries are lower confidence not because of cross-mapping, but because closely related species retain fewer unique diagnostic k-mers after LCA merging, pushing their signal closer to the detection threshold. Cross-kingdom co-detections are the highest-confidence signal and biologically unambiguous.

## Reference databases

- **PHI-base**: [phi-base.org](https://phi-base.org)
- **ICTV VMR**: [ictv.global/vmr](https://ictv.global/vmr/current)
- **NCBI STAT**: Katz et al. (2021) *J Bioinform Comput Biol* — [PMC8450716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8450716/)

## Author

Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)
