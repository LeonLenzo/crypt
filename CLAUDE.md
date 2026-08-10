# crypt — Cryptic Co-infection Mining from SRA RNA-seq

## What this is
A Python pipeline that mines NCBI SRA for RNA-seq runs containing mixed
host/pathogen libraries, then uses NCBI STAT pre-computed k-mer taxonomy
to detect cryptic secondary co-infecting organisms not reported by the
original study.

**Scope: plant hosts only** (Viridiplantae / PHI-base plant entries).

**Core hypothesis**: field-collected plant samples sequenced for a single
pathogen (or as host transcriptomics) contain unreported co-infecting
organisms detectable via STAT k-mer taxonomy without re-aligning reads.

## Two search strategies

| Mode | Abbrev | Library organism | STAT retention gate |
|---|---|---|---|
| microbe-as-library | MAL | PHI-base plant pathogen | Viridiplantae reads ≥ 1% (confirms in planta) |
| host-as-library | HAL | PHI-base plant host species | ≥1 PHI-base pathogen or plant virus detected in STAT |

**Key design principles**:
- PHI-base anchors both the SRA query and the STAT retention gate; ICTV adds plant viruses
- No YAML config, no BioProject filtering — queries are hardcoded; STAT is the only filter
- Every run is screened independently (field/experimental variability means runs from the same BioProject can differ)
- Species-level PHI-base confirmation is sufficient (f. sp. / cultivar not required)
- STAT uses organism **names not taxids** — resolved via the name lookup table in phibase_db.json
- `txid{n}[Host]` does NOT work in SRA esearch (text field, not taxonomy-linked); in-planta confirmation is done via STAT host_pct instead

## Module structure

Three modules, each with scripts in the module root, figure scripts in `module/figures/`,
and all generated output under `module/output/{data,figures}/`. Run all scripts from `crypt/` root.

`_util.py` at crypt/ root — shared: `_Tee`, `http_get`, `load_json`, `save_json`.

### stat/  — STAT k-mer screening

```
stat/build.py        Build PHI-base + ICTV reference DB. MUST use system python3 (ete3).
                     Output: stat/output/data/phibase_db.json
stat/fetch_runs.py   Fetch SRA RunInfo + NCBI STAT taxonomy. --mode {mal|hal}
                     Output: stat/output/data/{mode}_runs.json, stat_cache.jsonl
stat/filter_runs.py  Gate → validate → crypt detection. KINGDOM_THRESHOLDS here.
                     --mode {mal|hal|both}  --skip-validate after threshold review
                     Output: stat/output/data/runs.tsv

stat/figures/
  prep_scatter.py    Prepare scatter_data.tsv for scatter.R
  scatter.R          MAL and HAL scatter panels
stat/output/figures/scatter/
```

### kraken/  — Kraken2 species-level validation

```
kraken/build.py      Download PHI-base euk pathogen CDS; BBDuk masking; build Kraken2 DB.
                     Input: stat/output/data/phibase_db.json
                     Output: kraken/output/logs/; DB on Setonix scratch
kraken/classify.py   Stream reads from ENA FTP → kraken2 → parse species-level %.
                     Confidence=0.15, min-hit-groups=3.
                     Input: stat/output/data/{mode}_runs.json  OR  --run-list PATH
                     Output: kraken/output/data/kraken_cache.jsonl
kraken/screen_refs.py   Screen reference assemblies for best accession per taxid.
kraken/control/
  sample.py          Build 2,473-run stratified validation set (strata A–I).
                     Input: stat/output/data/, metadata/output/data/
                     Output: kraken/control/output/data/

kraken/figures/
  compare_host_pathogen.py   DB comparison figures (masked vs unmasked, STAT vs kraken)
kraken/output/figures/db_comparison/
```

### metadata/  — BioProject/BioSample/literature enrichment + LLM classification

```
metadata/fetch_xml.py    Fetch BioProject XML + BioSample XML attributes (POST efetch).
                         Input: stat/output/data/runs.tsv
                         Output: metadata/output/data/{bioprojects,biosamples}.json
metadata/fetch_lit.py    Fetch PMIDs + PMC methods (6 strategies, short-circuit at first hit).
                         --retry to re-attempt no-PMID entries.
                         Output: metadata/output/data/literature.json
metadata/filter_kw.py    Keyword classification (treatment × study_setting). --hc for diff-genus only.
                         Output: metadata/output/data/biosample_kw.tsv
metadata/llm_classify.py LLM BioProject classification via OpenAI gpt-4o-mini. OPENAI_API_KEY required.
                         --rerun-all to re-classify; --workers N (default 8).
                         Output: metadata/output/data/bioproject_llm.tsv

metadata/figures/
  mal_guilds.py + .R              Co-infection guild network (MAL+HAL)
  field_hc_guilds.py + .R         Field HC guild network + hull variants
  crypt_host_tree.py + .R         Plant host fan tree
  sample_funnel.py                5-column Sankey funnel
  coinf_rate.R / kingdom_comp.R / novel_heatmap.R / study_design_figures.R
metadata/output/figures/{guilds,host_tree,sankey,coinf_rate,kingdom_comp,novel_heatmap,study_design}/
```

### Running the pipeline

```bash
# stat module — run from crypt/
python3 stat/build.py                        # system python3 required (ete3)
python stat/fetch_runs.py --mode mal
python stat/fetch_runs.py --mode hal
python stat/filter_runs.py                   # both modes; --skip-validate after review

# metadata module
python metadata/fetch_xml.py
python metadata/fetch_lit.py
python metadata/fetch_lit.py --retry         # retry no-PMID entries
python metadata/filter_kw.py
python metadata/filter_kw.py --hc           # restrict to same_genus_secondary=False
python metadata/llm_classify.py
python metadata/llm_classify.py --rerun-all

# kraken module (Setonix)
python kraken/build.py --db-dir /scratch/... --genomes-dir /scratch/...
python kraken/classify.py --mode both --db /scratch/...
python kraken/control/sample.py             # build validation set first
python kraken/classify.py --run-list kraken/control/output/data/run_ids.txt --db /scratch/...
```

## Reference database (`stat/output/data/phibase_db.json`)

Built by `stat/build.py` (system python3 required). PHI-base: 205 pathogen seeds + 180 host
seeds. ICTV VMR: ~2,630 plant viruses. 15,372 pathogen taxids total; 21,352 name entries.
See `stat/build.py` for JSON schema; `PhibaseDB` class in `stat/filter_runs.py` for usage.
`stat/output/data/phibase_db.json` is an input to `kraken/build.py` and `metadata/filter_kw.py`.

CRITICAL: `_expand_taxids()` MUST use `intermediate_nodes=True` — default silently drops
species nodes that have named strains as children (PVY, ToBRFV, f.sp. taxa etc.).

## STAT approach

See `stat/fetch_runs.py` for fetch implementation and `stat/filter_runs.py` for detection logic.
Key algorithm: `specific_hits()` — leaf-level species detection via count nesting.
Cache: `stat/output/data/stat_cache.jsonl` (append-only, one line per run) + `stat_cache_index.txt`.

## STAT detection thresholds (KINGDOM_THRESHOLDS in stat/filter_runs.py)

Eukaryotic pathogens only (bacteria/viruses excluded — see Notes):

| Kingdom  | Detection threshold |
|----------|---------------------|
| Fungi    | ≥ 0.5%              |
| Oomycota | ≥ 0.5%              |
| Nematoda | ≥ 1.0%              |

MAL in-planta gate: `MAL_MIN_HOST_PCT = 1.0` (Viridiplantae reads ≥ 1%).
HAL pathogen gate:  `HAL_MIN_PATHOGEN_PCT = 1.0` (eukaryotic PHI-base pathogen ≥ 1%).

LibrarySource pre-filter: `_RNA_SOURCES = {TRANSCRIPTOMIC, TRANSCRIPTOMIC SINGLE CELL,
  METATRANSCRIPTOMIC, VIRAL RNA}` — applied before the gate; excludes GENOMIC / METAGENOMIC /
  OTHER runs that appear in RNA-Seq[Strategy] results due to submitter labelling errors.
  Removed 1,419 MAL + 7,343 HAL runs (e.g. PRJNA1111826 — Puccinia triticina genome sequencing
  mislabelled as RNA-Seq; `LibraryStrategy=RNA-Seq` but `LibrarySource=GENOMIC`).

## Output schemas

See file headers or script docstrings for full column definitions.

- `stat/output/data/runs.tsv` (stat/filter_runs.py): key cols — Run, mode, BioSample, BioProject,
  host, host_pct, stat_pathogens, interaction_status, co_infection_flag, same_genus_secondary,
  biosample_representative, fungi_pct/oomycete_pct/nematode_pct, analyzed.
  `interaction_status == "novel_host_range"` is the primary signal.
  Filter `biosample_representative == True` for one row per biological sample.

- `metadata/output/data/biosample_kw.tsv` (metadata/filter_kw.py): one row per BioSample;
  BioProject metadata repeated. Join with `bioproject_llm.tsv` on BioProject for LLM columns.
  Primary analysis input. BioSample XML coverage: geo_loc_name 61%, tissue 56%, collection_date 42%.
  Use `llm_treatment` preferentially over keyword `treatment` (kw has ~51% contamination).

- `stat/output/data/phibase_db.json` (stat/build.py): kingdom-separated taxid maps.

## Entrez / external API credentials
- `NCBI_API_KEY` env var → 10 req/s (9 used for Entrez); without: 2.5 req/s
- `S2_API_KEY` env var → Semantic Scholar 10 req/s; without: 1 req/s
  Apply free at semanticscholar.org/product/api (Academic Graph API)
- STAT endpoint rate: empirically ~2–10 req/s; RATE=15 + 32 workers configured
- EBI (EuropePMC, ENA): shared 2 req/s rate limiter (_ext_wait, gap=0.5s)
- http_get accepts no_retry_429=True — fail fast on 429 for external APIs
- User-Agent: `crypt/fetch_runs (leon.lenzo@curtin.edu.au)` etc.

## Output / cache layout

```
stat/output/data/        phibase_db.json, stat_cache.jsonl, runs.tsv, {mal,hal}_runs.json
stat/output/figures/     scatter/

metadata/output/data/    bioprojects.json, biosamples.json, literature.json,
                         biosample_kw.tsv, bioproject_llm.tsv
metadata/output/figures/ guilds/, host_tree/, sankey/, coinf_rate/, kingdom_comp/,
                         novel_heatmap/, study_design/

kraken/output/data/      kraken_cache.jsonl, kraken_cache_index.txt
kraken/output/figures/   db_comparison/
kraken/control/output/   data/run_ids.txt, control_runs.tsv, manifest.tsv
```

## Actual run results (2026-08-03, eukaryotic-only scope)

**stat/fetch_runs.py / stat/filter_runs.py (euk-only pivot applied 2026-08-03):**

| Mode | Runs fetched | Non-RNA excluded | Screened | Gate pass | Gate pass rate |
|------|-------------|-----------------|---------|-----------|----------------|
| MAL  | 48,418      | 1,419           | 46,315  | 5,223     | 11.3%          |
| HAL  | 559,950     | 7,343           | 546,816 | 5,772     | 1.1%           |

- stat/output/data/runs.tsv — 10,995 rows; BioProjects: 1,285; 608,368 STAT cache entries
- metadata/output/data/biosamples.json — 9,002 entries; geo 61%, tissue 56%, collection_date 42%
- metadata/output/data/biosample_kw.tsv — 9,002 rows; 90% named_host resolved
- metadata/output/data/bioproject_llm.tsv — 1,285 BPs
  - llm_treatment: single=895, host_study=219, abiotic_stress=94, coinf_experiment=21, unclear=19, surveillance=6
  - llm_setting: lab=701, unclear=431, field=122

**Headline co-infection rates (field vs controlled, LLM-classified):**
- Field: 22.7% (277/1,219 classified)  |  Controlled: 11.4% (531/4,650 classified)

## Kraken2 validation pipeline

STAT has a confirmed blind spot for PST (euk_pct=0% for all 15 PST runs in pilot;
Kraken2+kallisto show 65-68%). Kraken2 uses CDS k-mers; species-diagnostic masking via
BBDuk removes cross-species noise before DB build.

**BBDuk masking**: pathogen-vs-pathogen only — each pathogen masked against k-mers shared
with any other pathogen (kmercountexact mincount=2). Every k-mer in DB is species-diagnostic
among pathogens. Eliminates same-order (Melampsora in PST runs), same-family cross-hits.

**Pathogen-only DB (Setonix, 2026-08-07):**
- `/scratch/pawsey1168/llenzo/kraken/db_pathogens/` — 1.9G hash, fungi + oomycetes only
- No host sequences — host identity from SRA metadata; pct_classified = pathogen burden directly

**Control validation set (2026-08-08):** 2,473 stratified runs across strata B–I (stratum A
stub pending non-plant source). Run list at `kraken/control/output/data/run_ids.txt`.
Submit on Setonix: `git pull && sbatch kraken/slurm/kraken_classify_control.slurm`

## Performance notes

- **stat/fetch_runs SRA RunInfo fetch**: parallel POST requests (ThreadPoolExecutor, MAX_WORKERS=8).
  GET requests with 500 UIDs cause HTTP 414; POST avoids this. ~20 req/s achieved.
- **stat/fetch_runs STAT fetch**: STAT_RATE=15 req/s configured; apparent variability (~2–10 req/s
  observed) was caused by thread hangs (workers blocking on network I/O), not server-side
  time-of-day variation. MAL (~48k): ~1.5 hrs. HAL (~559k): ~3–4 days with hangs factored in.
- **Shared stat_cache**: do NOT run MAL and HAL stat/fetch_runs simultaneously — last writer
  wins on cache saves. Chain HAL after MAL:
  `until ls stat/output/logs/mal_summary.latest 2>/dev/null; do sleep 30; done`

## Notes

- `specific_hits()` leaf detection in `stat/filter_runs.py` is the core non-trivial algorithm
- `_expand_taxids()` MUST use `intermediate_nodes=True` — default returns only leaf taxa,
  silently dropping species nodes that have named strains as children (e.g. PVY, f.sp. taxa).
- MAL `host` column: `detect_host_species()` uses a 305,355-name Viridiplantae allowlist
  (`viridiplantae_names` in phibase_db.json). Broad-clade hits (e.g. "Viridiplantae",
  "Triticinae") filtered via `BROAD_CLADE_NAMES` in guild figure scripts; `named_host`
  used as fallback. `Aegilops tauschii` remapped to `Triticum aestivum` in R hull layer.
- **Fusarium blind spot**: 94% of Fusarium co-infections are same-genus → removed by HC filter.
  Only 6 Fusarium appearances in field HC. Consider separate same-genus-allowed figure.
- `interaction_status == "novel_host_range"` is the key signal; `llm_treatment == "coinf_experiment"`
  flags intentional co-infection designs for manual exclusion.
- `metadata/filter_kw.py --hc` restricts to same_genus_secondary=False; default includes all.
- Use `llm_treatment` (`metadata/llm_classify.py`) over keyword `treatment` — kw has ~51%
  contamination in host_study category (pathogen/disease language co-occurs).
- `metadata/fetch_lit.py` PMID search: 6 strategies short-circuit at first hit.
  S2_API_KEY required for Semantic Scholar (1 req/s); skips without key.
- `stat/fetch_runs.py` STAT fetch hangs overnight: SIGINT doesn't exit cleanly — requires SIGKILL.
  Workaround: tmux + manual restart. Run in background with `tee -a`.
- **Eukaryotic pathogens only (2026-08-03 pivot)**: bacteria/viruses removed from detection.
  PolyA+ selection depletes bacterial mRNA and most plant virus families.
  KINGDOM_THRESHOLDS: Fungi ≥0.5%, Oomycota ≥0.5%, Nematoda ≥1.0%.
- **STAT blind spot**: STAT euk_pct = 0% for ALL PST runs in pilot (15/15 runs).
  Kraken2 + kallisto agree at 65-68% for high-tier runs. See memory/kraken_pilot_results.md.

## Active task list (2026-08-10)

1. [ ] **Kraken2 control validation** (ACTIVE) — 2,473-run stratified set built; job needs
       `git pull && sbatch kraken/slurm/kraken_classify_control.slurm` on Setonix.

2. [ ] **Kraken2 production run** — update `kraken/slurm/kraken_classify.slurm` to use
       `db_pathogens` and `--mem=16G`; run after control validation confirms DB.

3. [ ] **same_family_secondary flag** — add to runs.tsv; precompute in stat/build.py via ete3.
       Three confidence tiers: same_genus > same_family > diff-family.

4. [ ] **Submission timeline + PMID figure** — `bp_submission_date` on x-axis.
       `metadata/output/figures/timeline/`.

5. [ ] **Systematic no-PMID BioProject review** — ~886 no-PMID BioProjects. Try CrossRef title search.

6. [ ] **Tissue normalization** — normalise Leaf/leaf/leaves case variation in BioSample XML.

## Deferred analysis ideas

1. **Field prevalence analysis** — DONE (2026-07-29, metadata/output/figures/study_design/).
   Field single-pathogen BPs: 21% overall / 12% hc coinf rate.
3. **Host susceptibility landscape** — co-infection rate per host on crypt_host_tree.
4. **Gate failure characterisation** — what are the ~85% MAL runs failing Viridiplantae gate?
5. **Viral threshold sensitivity** — re-run stat/filter_runs at 5% vs 10% threshold.

## Target journal

New Phytologist (IF ~10) or ISME Journal (IF ~11). Frame as biology discovery
(novel host-pathogen interactions + phylogenetic host susceptibility structure),
not as a methodology paper.
