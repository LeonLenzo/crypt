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

## Module architecture

Fully standalone scripts (stdlib + json only, except ete3/openpyxl in 00_build.py).
`_util.py` — shared: `_Tee`, `http_get`, `load_json`, `save_json`.

```
00_build.py          Build PHI-base + ICTV reference DB. MUST use system python3 (ete3).
01_fetch_runs.py     Fetch SRA RunInfo + NCBI STAT taxonomy. --mode {mal|hal}
02_filter_runs.py    Gate → validate → crypt detection. KINGDOM_THRESHOLDS here.
                     --mode {mal|hal|both}  --skip-validate after threshold review
03a_fetch_xml.py     Fetch BioProject XML + BioSample XML attributes (POST efetch).
03b_fetch_literature.py  Fetch PMIDs + PMC methods (6 strategies, short-circuit at first hit).
                     --retry to re-attempt no-PMID entries.
04_filter_kw.py      Keyword classification (treatment × study_setting). --hc for diff-genus only.
05_llm_classify.py   LLM BioProject classification via OpenAI gpt-4o-mini. OPENAI_API_KEY required.
                     --rerun-all to re-classify; --workers N (default 8).
```

### Running the pipeline

```bash
# Step 0: build reference DB (system python3 required)
python3 00_build.py

# Step 1: fetch SRA runs + STAT (one mode at a time — shared stat cache)
python 01_fetch_runs.py --mode mal
python 01_fetch_runs.py --mode hal

# Step 2: filter + detect (both modes by default, unified output)
python 02_filter_runs.py                     # validate + gate + crypt, both modes
python 02_filter_runs.py --skip-validate     # after manual threshold review

# Step 3a: BioProject XML + BioSample XML attributes
python 03a_fetch_xml.py

# Step 3b: BioProject PMIDs + PMC methods sections
python 03b_fetch_literature.py
python 03b_fetch_literature.py --retry       # retry no-PMID entries after new strategy

# Step 4: keyword classification + BioSample join → biosample_kw.tsv
python 04_filter_kw.py
python 04_filter_kw.py --hc                 # restrict to same_genus_secondary=False

# Step 5: LLM classification → bioproject_llm.tsv (requires OPENAI_API_KEY)
python 05_llm_classify.py
python 05_llm_classify.py --rerun-all       # re-classify all, ignore cache
```

## Reference database (`output/00_build/phibase_db.json`)

Built by `00_build.py` (system python3 required). PHI-base: 205 pathogen seeds + 180 host
seeds. ICTV VMR: ~2,630 plant viruses. 15,372 pathogen taxids total; 21,352 name entries.
See `00_build.py` for JSON schema; `PhibaseDB` class in `02_filter_runs.py` for usage.

CRITICAL: `_expand_taxids()` MUST use `intermediate_nodes=True` — default silently drops
species nodes that have named strains as children (PVY, ToBRFV, f.sp. taxa etc.).

## STAT approach

See `01_fetch_runs.py` for fetch implementation and `02_filter_runs.py` for detection logic.
Key algorithm: `specific_hits()` — leaf-level species detection via count nesting.
Cache: `stat_cache.jsonl` (append-only, one line per run) + `stat_cache_index.txt` (fast startup).

## STAT detection thresholds (KINGDOM_THRESHOLDS in 02_filter_runs.py)

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

- `runs.tsv` (02_filter_runs.py): key cols — Run, mode, BioSample, BioProject, host, host_pct,
  stat_pathogens, interaction_status, co_infection_flag, same_genus_secondary,
  biosample_representative, fungi_pct/oomycete_pct/nematode_pct, analyzed.
  `interaction_status == "novel_host_range"` is the primary signal.
  Filter `biosample_representative == True` for one row per biological sample.

- `biosample_kw.tsv` (04_filter_kw.py): one row per BioSample; BioProject metadata repeated.
  Join with `bioproject_llm.tsv` on BioProject for LLM columns. Primary analysis input.
  BioSample XML coverage: geo_loc_name 61%, tissue 56%, collection_date 42%.
  Use `llm_treatment` preferentially over keyword `treatment` (kw has ~51% contamination).

- `phibase_db.json` (00_build.py): kingdom-separated taxid maps. See 00_build.py for schema.

## Entrez / external API credentials
- `NCBI_API_KEY` env var → 10 req/s (9 used for Entrez); without: 2.5 req/s
- `S2_API_KEY` env var → Semantic Scholar 10 req/s; without: 1 req/s
  Apply free at semanticscholar.org/product/api (Academic Graph API)
- STAT endpoint rate: empirically ~2–10 req/s; RATE=15 + 32 workers configured
- EBI (EuropePMC, ENA): shared 2 req/s rate limiter (_ext_wait, gap=0.5s)
- http_get accepts no_retry_429=True — fail fast on 429 for external APIs
- User-Agent: `crypt/01_fetch_runs (leon.lenzo@curtin.edu.au)` etc.

## Output / cache layout

Each step: `output/{step}/data/` (outputs + caches) and `output/{step}/logs/` (timestamped subdirs + `latest` symlink).
Key files: `output/01_fetch_runs/data/stat_cache.jsonl`, `runs.tsv`, `biosample_kw.tsv`,
`bioproject_llm.tsv`, `bioprojects.json`, `biosamples.json`, `literature.json`.
Figures: `output/figures/{scatter,guilds,host_tree,sankey,study_design,...}/`.

## Actual run results (2026-08-03, eukaryotic-only scope)

**01_fetch_runs.py / 02_filter_runs.py (euk-only pivot applied 2026-08-03):**

| Mode | Runs fetched | Non-RNA excluded | Screened | Gate pass | Gate pass rate |
|------|-------------|-----------------|---------|-----------|----------------|
| MAL  | 48,418      | 1,419           | 46,315  | 5,223     | 11.3%          |
| HAL  | 559,950     | 7,343           | 546,816 | 5,772     | 1.1%           |

- output/02_filter_runs/data/runs.tsv — 10,995 rows; BioProjects: 1,285; 608,368 STAT cache entries
- output/03a_fetch_xml/data/biosamples.json — 9,002 entries; geo 61%, tissue 56%, collection_date 42%
- output/04_filter_kw/data/biosample_kw.tsv — 9,002 rows; 90% named_host resolved
- output/05_llm_classify/data/bioproject_llm.tsv — 1,285 BPs
  - llm_treatment: single=895, host_study=219, abiotic_stress=94, coinf_experiment=21, unclear=19, surveillance=6
  - llm_setting: lab=701, unclear=431, field=122

**Headline co-infection rates (field vs controlled, LLM-classified):**
- Field: 22.7% (277/1,219 classified)  |  Controlled: 11.4% (531/4,650 classified)

**Figures (output/figures/{subdir}/):**
- scatter/scatter.R: host% vs pathogen%, coloured single/coinf/novel, shape=MAL/HAL
- host_tree/crypt_host_tree.py + .R: fan tree 173 hosts, bars = n_single/n_multi
- guilds/mal_guilds.py + .R: co-occurrence network (MAL+HAL)
- sankey/sample_funnel.py: 5-column funnel; field 22.7%, controlled 11.4%
- study_design/treat_setting_bar.R + single_field_lab.R: field vs lab breakdown

## Kraken2 validation pipeline (`kraken/`)

Replacement for STAT-based detection. STAT has a confirmed blind spot for PST (euk_pct=0%
for all 15 PST runs in pilot; Kraken2+kallisto show 65-68%). Kraken2 uses CDS k-mers;
species-diagnostic masking via BBDuk removes cross-species noise before DB build.

```
kraken/build.py      Download PHI-base euk pathogen CDS + plant host CDS; bidirectional
                     BBDuk masking to species-diagnostic k-mers; build Kraken2 DB.
                     --skip-bbduk to rebuild without masking.
kraken/classify.py   Stream reads from ENA FTP → kraken2 → parse species-level %.
                     Confidence=0.15, min-hit-groups=3.
                     Cache: output/kraken_classify/data/kraken_cache.jsonl
```

**BBDuk masking**: two passes — (1) pathogens masked against shared pathogen k-mers + all
host CDS; (2) hosts masked against shared host k-mers + all pathogen CDS. Result: every
k-mer in the DB is species-diagnostic. Eliminates same-order (Melampsora in PST runs),
same-family (Sorghum in maize), and cross-kingdom (plant CDS in fungal libraries) noise.

**Status (2026-08-06)**: scripts written; 40-run test done (confidence=0.15, mhg=3);
BBDuk masking implementation complete; DB rebuild pending. Full production (~593k runs)
requires Setonix. STAT vs Kraken comparison = paper figure.

## Performance notes

- **01_fetch_runs SRA RunInfo fetch**: parallel POST requests (ThreadPoolExecutor, MAX_WORKERS=8).
  GET requests with 500 UIDs cause HTTP 414; POST avoids this. ~20 req/s achieved.
- **01_fetch_runs STAT fetch**: STAT_RATE=15 req/s configured; apparent variability (~2–10 req/s
  observed) was caused by thread hangs (workers blocking on network I/O), not server-side
  time-of-day variation. MAL (~48k): ~1.5 hrs. HAL (~559k): ~3–4 days with hangs factored in.
- **Shared stat_cache**: do NOT run MAL and HAL 01_fetch_runs simultaneously — last writer
  wins on cache saves. Chain HAL after MAL:
  `until ls output/01_fetch_runs/logs/mal_summary.latest 2>/dev/null; do sleep 30; done`

## Notes

- `specific_hits()` leaf detection in 02_filter_runs.py is the core non-trivial algorithm
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
- 04_filter_kw.py `--hc` restricts to same_genus_secondary=False; default includes all.
- Use `llm_treatment` (05_llm_classify.py) over keyword `treatment` — kw has ~51% contamination
  in host_study category (pathogen/disease language co-occurs).
- 03b_fetch_literature.py PMID search: 6 strategies short-circuit at first hit.
  S2_API_KEY required for Semantic Scholar (1 req/s); skips without key.
- 01_fetch_runs.py STAT fetch hangs overnight: SIGINT doesn't exit cleanly — requires SIGKILL.
  Workaround: tmux + manual restart. Run in background with `tee -a`.
- **Eukaryotic pathogens only (2026-08-03 pivot)**: bacteria/viruses removed from detection.
  PolyA+ selection depletes bacterial mRNA and most plant virus families.
  KINGDOM_THRESHOLDS: Fungi ≥0.5%, Oomycota ≥0.5%, Nematoda ≥1.0%.
- **STAT blind spot**: STAT euk_pct = 0% for ALL PST runs in pilot (15/15 runs).
  Kraken2 + kallisto agree at 65-68% for high-tier runs. See memory/kraken_pilot_results.md.

## Active task list (2026-08-06)

1. [ ] **Kraken2 replacement for STAT** (ACTIVE) — `kraken/build.py` + `kraken/classify.py` written.
       40-run local test done (confidence=0.15, min-hit-groups=3). BBDuk bidirectional masking
       (species-diagnostic k-mers) implemented in build.py; DB rebuild pending. Full production
       run (~593k runs) requires Setonix. See memory/kraken_pipeline.md.

2. [ ] **same_family_secondary flag** — add to runs.tsv; precompute in 00_build.py via ete3.
       Three confidence tiers: same_genus > same_family > diff-family.

3. [ ] **Submission timeline + PMID figure** — `bp_submission_date` on x-axis. `output/figures/timeline/`.

4. [ ] **Systematic no-PMID BioProject review** — ~886 no-PMID BioProjects. Try CrossRef title search.

5. [ ] **Tissue normalization** — normalise Leaf/leaf/leaves case variation in BioSample XML.

## Deferred analysis ideas

1. **Field prevalence analysis** — DONE (2026-07-29, output/figures/study_design/). Field single-pathogen
   BPs: 21% overall / 12% hc coinf rate. Lab single-pathogen: 10% / 2% hc. See study design
   figures for full treatment × setting breakdown.
3. **Host susceptibility landscape** — co-infection rate per host on crypt_host_tree.
   Are some hosts phylogenetically more susceptible?
4. **Gate failure characterisation** — what are the ~85% MAL runs failing Viridiplantae gate?
5. **Viral threshold sensitivity** — re-run 02_filter_runs at 5% vs 10% threshold.
6. **Kraken2 orthogonal validation** — ACTIVE (2026-08-04). Moved from deferred to priority 2
   above. Assembly scan running; Setonix pipeline design in progress.

## Target journal

New Phytologist (IF ~10) or ISME Journal (IF ~11). Frame as biology discovery
(novel host-pathogen interactions + phylogenetic host susceptibility structure),
not as a methodology paper.
