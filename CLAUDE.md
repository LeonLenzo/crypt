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

All scripts are **fully standalone** — no shared package imports, stdlib + json only
(except ete3 and openpyxl in 00_build.py which requires system python3).
Shared utilities (http_get, load_json, save_json, _Tee) live in `_util.py`.

Legacy scripts (01_sra.py → 05_meta.py) are preserved in `scripts/legacy/` for reference.
The canonical pipeline is the seven scripts below.

```
_util.py      Shared utilities: _Tee (stdout+file tee), http_get, load_json, save_json.

00_build.py   Build PHI-base + ICTV reference DB.
              MUST run with system python3 (ete3/sqlite3 incompatibility
              with miniconda).
              Auto-downloads PHI-base CSV from GitHub if absent.
              Auto-downloads ICTV VMR Excel from ictv.global if absent.
              Classifies PHI-base pathogen seeds into kingdoms via ete3 lineage.
              Plant viruses: Host source in {plants, plants (S),
              invertebrates, plants} → resolved via ete3 to NCBI taxids.

01_fetch_runs.py   Fetch SRA run IDs + NCBI STAT taxonomy for MAL or HAL.
              Merges legacy 01_sra.py + 02_stat.py into a single data-gathering step.
              No retention gate — gate moved to 02_filter_runs.py.
              MAL: batches 205 seed pathogen taxids; HAL: 180 seed host taxids.
              [Organism:exp] lets SRA expand to all descendant strains.
              Parallel fetch: SRA_MAX_WORKERS=8 (RunInfo) + STAT_MAX_WORKERS=32 (STAT).
              Separate rate locks: _ri_lock (SRA) + _stat_lock (STAT).
              CACHE_MIN_PCT=0.1 filter applied at STAT write time.
              Output: output/01_fetch_runs/data/{mode}_runs.json + stat_cache.jsonl

02_filter_runs.py  Three phases in one script — gate → validate → crypt.
                Phase 1: retention gate
                  MAL: keep if Viridiplantae reads ≥ MAL_MIN_HOST_PCT (1%)
                  HAL: keep if any PHI-base pathogen/virus ≥ HAL_MIN_PATHOGEN_PCT (1%)
                Phase 2: validate distributions (always runs; --skip-validate to skip
                         after manual threshold review)
                  Per kingdom: breakpoint table, ASCII histogram, gap detection,
                  absolute count floor table. Outputs R-ready TSVs.
                Phase 3: crypt co-infection detection
                  specific_hits() leaf algorithm + PHI-base/ICTV cross-reference.
                  _interaction_status(): known / novel_host_range / novel_combination / unresolved
                  _annotate_biosamples(): biosample_n_runs + biosample_representative columns.
              --mode {mal|hal|both}  (default: both)
              Reads stat_cache.jsonl and {mode}_runs.json from output/01_fetch_runs/data/.
              KINGDOM_THRESHOLDS in this file — edit then --skip-validate to re-run crypt.
              Unified output: output/02_filter_runs/data/runs.tsv with `mode` column.
              BioSample dedup applied across both modes in the combined output.

03a_fetch_xml.py   Fetch BioProject XML and BioSample XML attributes.
              Reads runs.tsv for BioProject and BioSample accessions.
              BioProjects: bootstrapped from find_cache.json (03b source); only new ones
                fetched from NCBI BioProject XML (title, description, submission_date).
              BioSamples: batch efetch via POST (avoids HTTP 414); 300 per request.
                Attributes: tissue, geo_loc_name, collection_date, lat_lon, host,
                isolation_source, dev_stage. Coverage: geo_loc_name ~61%, tissue ~55%,
                collection_date ~42%.
              Stubs (BioSamples with no public attributes) kept as {accession: acc}.
              Output filtered to accessions in runs.tsv only (NCBI may return extras).
              Output: output/03a_fetch_xml/data/bioprojects.json
                      output/03a_fetch_xml/data/biosamples.json (11,117 entries)

03b_fetch_literature.py   Fetch BioProject PMIDs and PMC full-text methods sections.
              Reuses existing find_cache.json (output/03_fetch_meta/data/) — 1,802 entries
              already cached; only new BioProjects fetched.
              6 PMID strategies (short-circuit, stop at first hit; ordered by yield):
                1. BioProject XML <Publication> (15%)
                2. PMC full-text search (40%)
                3. Europe PMC (preprints, supplementary mentions; 35%)
                4. ENA XML API for PRJEB accessions
                5. Semantic Scholar (S2_API_KEY required; 1 req/s)
                6. PubMed text search (last resort)
              After PMID: fetch PMC ID + methods section (up to 8,000 chars).
              Output: output/03b_fetch_literature/data/literature.json
                      (keyed by BioProject; fields: primary_pmid, primary_pub_date,
                      primary_publication, abstract, methods_text, pmcid,
                      n_papers_found, pmid_source)

04_filter_kw.py   Join BioSample metadata + keyword study design classification.
              BioSample is the unit of organisation — one row per biosample_representative
              run; BioProject metadata joined as repeated columns.
              Named host: BioSample XML host attr (priority) → title/abstract scanning.
              study_setting: geo_loc_name presence in BioSample XML upgrades
                "unclear" → "field" when no lab/field keywords fire.
              Two-axis classification (keyword-based):
                treatment: coinf_experiment > abiotic_stress > host_study > single > unclear
                study_setting: lab > field (or geo_loc_name) > unclear > no_data
              --hc flag restricts to same_genus_secondary=False.
              Output: output/04_filter_kw/data/biosample_kw.tsv (11,117 rows)

05_llm_classify.py   LLM-based BioProject classification via OpenAI API.
              Reads biosample_kw.tsv (grouped by BioProject) + bioprojects.json +
                literature.json.
              Prompt includes keyword classifier output — LLM agrees/disagrees.
              Per BioProject:
                llm_treatment, llm_study_setting (updated vocabulary; adds surveillance)
                llm_named_pathogen, llm_named_host — extracted from text
                llm_tissue — inferred tissue type(s); fallback where XML tissue blank
                llm_kw_treatment_agree, llm_kw_setting_agree — override flags
                llm_confidence, llm_rationale
              Cache: output/05_llm_classify/data/classify_cache.jsonl (append-only).
              Re-runs any entry missing required fields (add field → automatic upgrade).
              Requires: OPENAI_API_KEY env var; model gpt-4o-mini (~$2/full run).
              --rerun-all to ignore cache; --workers N for parallelism (default 8).
              Output: output/05_llm_classify/data/bioproject_llm.tsv (1,754 rows)
              Joins with biosample_kw.tsv on BioProject for downstream analysis.

scripts/diagnostics/benchmark_strategies.py
              Diagnostic: tests each 03b PMID strategy on unlabelled BioProjects.
              Measures discovery yield (% returning any PMID) + avg time per strategy.
              Run from crypt/: python scripts/diagnostics/benchmark_strategies.py [--n N] [--seed N]
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

## PHI-base + ICTV reference database (`output/00_build/phibase_db.json`)

Built by `00_build.py` using ete3 NCBITaxa (local SQLite at `~/.etetoolkit/taxa.sqlite`).

**PHI-base** (fungi, bacteria, oomycetes, nematodes):
- 205 pathogen seed species → ~7,985 expanded taxids (strains, subspecies, f. sp.)
- 180 host seed species → ~1,869 expanded taxids (cultivars, varieties)
- Seeds classified into kingdoms via ete3 lineage (Fungi/Bacteria/Oomycota/Nematoda)
- Interaction map: PHI-base host × pathogen pairs (species-level, bidirectional)

**ICTV VMR** (plant viruses):
- ~2,630 plant virus species (Host source: plants / invertebrates, plants)
- Resolved to NCBI taxids via `ncbi.get_name_translator()` — unresolved names skipped
- Expanded via `_expand_taxids()` using `intermediate_nodes=True` — CRITICAL: without
  this, species-level taxids that have named strains as children (e.g. "Potato virus Y"
  taxid 12216, "Tomato brown rugose fruit virus" taxid 1761477) are treated as internal
  nodes and silently excluded. NCBI/ICTV maintain parallel old-style and new-style taxids
  for the same viruses; the ICTV seed resolves to the new-style taxid (e.g. 3240642 for
  "Potyvirus yituberosi") while SRA metadata uses the old-style name ("Potato virus Y")
  which lives as an intermediate node under the new-style taxid. `intermediate_nodes=True`
  captures both. The same issue applies to f. sp. taxa (e.g. "Puccinia graminis f. sp.
  tritici") — they are intermediate nodes under the species seed.
- Viruses mapped to all 180 PHI-base plant host seeds in `pathogen_to_hosts`
  → `known_interaction(virus, any_plant_host)` returns True
- ICTV species names added to `name_to_taxid` as supplement to NCBI names
- Auto-downloaded from `https://ictv.global/vmr/current` (redirects to latest release)

**Combined DB stats** (post-rebuild 2026-07-25, with intermediate_nodes=True fix):
- 15,372 pathogen taxids (8,250 PHI-base + 7,122 ICTV viruses)
- 21,352 name entries — scientific names + synonyms + ICTV names, lowercase

**JSON structure (kingdom-separated schema):**
```
fungal_to_seed       {str(expanded_taxid): seed_taxid}  — PHI-base fungi
bacterial_to_seed    {str(expanded_taxid): seed_taxid}  — PHI-base bacteria
oomycete_to_seed     {str(expanded_taxid): seed_taxid}  — PHI-base oomycetes
nematode_to_seed     {str(expanded_taxid): seed_taxid}  — PHI-base nematodes
virus_to_seed        {str(expanded_taxid): seed_taxid}  — ICTV plant viruses
host_to_seed         {str(expanded_taxid): seed_taxid}  — plant hosts
contaminant_taxids   [taxid, ...]  — phiX174, E.coli K12, etc.
name_to_taxid        {lowercase_name: taxid}
taxid_to_name        {str(taxid): canonical_name}
pathogen_to_hosts    {str(seed_pathogen_taxid): [seed_host_taxids]}
meta                 build date, counts, scope, vmr_path
```
`PhibaseDB` (in 02_filter_runs.py) derives `pathogen_taxids`, `host_taxids`,
`host_to_pathogens` at load time from the kingdom dicts — they are NOT stored in JSON.

## STAT approach

- **Parallel fetch** with global timestamp rate-lock (not per-worker sleep):
  `_rate_lock` + `_rate_last` ensures ≤ RATE requests/s total across all threads
- **STAT endpoint** (`trace.ncbi.nlm.nih.gov`) is separate from Entrez — STAT_RATE=15,
  MAX_WORKERS=32; effective throughput limited by thread hang bug (see Performance notes)
- **One-pass species resolution**: parse all kingdoms in a single pass per run
- **`specific_hits(table, node, analyzed)`** — leaf-level species detection:
  1. Find all entries with `total_count ≤ node_count` (i.e., nested under this kingdom)
  2. Detect "leaf" counts: a count is a leaf if no other count falls in `[LEAF_FRAC*c, c)`
  3. Prefer 2-word binomials over clade names or strain strings
  4. Returns `[(name, pct), ...]` sorted by pct descending
- **Mycoviruses excluded** from pathogen detection (keywords: mycovirus, mitovirus,
  hypovirus, etc. — they infect the fungus, not the plant)
- STAT shared cache: `stat_cache.jsonl` (one line per run, append-only) +
  `stat_cache_index.txt` (accessions only, for fast startup). Old
  `stat_cache.json` design was abandoned: at 8.6 GB it caused OOM during
  saves (json.dumps built a full in-memory string copy) and couldn't be
  appended to without rewriting the whole file.

## STAT detection thresholds (KINGDOM_THRESHOLDS in 02_filter_runs.py)

| Kingdom  | Detection threshold |
|----------|---------------------|
| Fungi    | ≥ 1.0%              |
| Viruses  | ≥ 10.0%             |
| Bacteria | ≥ 5.0%              |
| Oomycota | ≥ 0.5%              |
| Nematoda | ≥ 1.0%              |

MAL in-planta gate: `MAL_MIN_HOST_PCT = 1.0` (Viridiplantae reads ≥ 1%).
HAL pathogen gate:  `HAL_MIN_PATHOGEN_PCT = 1.0` (any PHI-base pathogen/virus ≥ 1%).

LibrarySource pre-filter: `_RNA_SOURCES = {TRANSCRIPTOMIC, TRANSCRIPTOMIC SINGLE CELL,
  METATRANSCRIPTOMIC, VIRAL RNA}` — applied before the gate; excludes GENOMIC / METAGENOMIC /
  OTHER runs that appear in RNA-Seq[Strategy] results due to submitter labelling errors.
  Removed 1,419 MAL + 7,343 HAL runs (e.g. PRJNA1111826 — Puccinia triticina genome sequencing
  mislabelled as RNA-Seq; `LibraryStrategy=RNA-Seq` but `LibrarySource=GENOMIC`).

## Output schema (02_filter_runs.py → `runs.tsv`)

| Column | Notes |
|---|---|
| `Run` | SRR accession |
| `mode` | `mal` or `hal` |
| `BioSample` | NCBI BioSample accession |
| `BioProject`, `SRAStudy`, `Platform` | RunInfo passthrough |
| `host` | top Viridiplantae leaf species in STAT |
| `host_pct` | Viridiplantae % of analysed reads |
| `library_organism` | library source organism (ScientificName from RunInfo) |
| `library_detected` | True if library organism appears in STAT output |
| `stat_pathogens` | semicolon-separated PHI-base/ICTV pathogens detected in STAT |
| `stat_hosts` | semicolon-separated Viridiplantae host species detected in STAT |
| `n_pathogens` | count of stat_pathogens |
| `interaction_status` | known / novel_host_range / novel_combination / unresolved |
| `co_infection_flag` | single / multi_species / multi_kingdom |
| `same_genus_secondary` | True if any STAT pathogen shares genus with library organism |
| `biosample_n_runs` | runs sharing this BioSample |
| `biosample_representative` | True for highest-analyzed run per BioSample |
| `fungi_pct`, `virus_pct`, `bacteria_pct`, `oomycete_pct`, `nematode_pct` | kingdom % from STAT |
| `analyzed` | total reads analysed by STAT |

`co_infection_flag` priority: `multi_kingdom` > `multi_species` > `single`.
`interaction_status`: `novel_host_range` is the primary signal — pathogen known to PHI-base
  but not recorded on this host. `unresolved` means taxid lookup failed.
  For plant viruses, `known` for any PHI-base plant host (ICTV viruses mapped to all 180 host seeds).
`same_genus_secondary`: True when first word of any `stat_pathogens` entry matches genus
  of `library_organism`. Use `same_genus_secondary == False` for highest-confidence detections.
`biosample_representative`: filter to True for one row per biological sample —
  used for paper-level co-infection rate calculations.

## Output schema (04_filter_kw.py → `biosample_kw.tsv`)

One row per BioSample (biosample_representative rows from runs.tsv), with BioProject metadata
joined as repeated columns. This is the primary analysis input for all downstream figures.

| Column group | Columns |
|---|---|
| BioSample identity | `BioSample`, `BioProject`, `SRAStudy`, `Run`, `mode` |
| BioSample XML | `tissue`, `geo_loc_name`, `collection_date`, `lat_lon`, `bs_host`, `isolation_source`, `dev_stage` |
| Co-infection | `n_runs`, `co_infection_flag`, `same_genus_secondary` |
| Organisms | `library_organism`, `library_detected`, `stat_pathogens`, `stat_hosts`, `n_pathogens` |
| STAT summary | `interaction_status`, `host_pct`, `analyzed` |
| Named host | `named_host`, `named_host_source` (biosample_xml / title / abstract / description / methods) |
| Keyword classification | `treatment`, `treatment_keywords`, `study_setting`, `setting_keywords` |
| BioProject text | `title`, `description`, `abstract` |
| Literature | `primary_pmid`, `primary_pub_date`, `primary_publication`, `n_papers_found`, `pmid_source`, `bp_submission_date` |

BioSample XML coverage (03a_fetch_xml.py results): geo_loc_name 61%, tissue 56%, collection_date 42%.
Named host coverage (04_filter_kw.py): 90% total; 15% from biosample_xml (structured), 75% from text.
Joins with bioproject_llm.tsv on `BioProject` for LLM-classified columns.

## k-mer specificity and same-genus secondaries

STAT uses 32 bp k-mers (PMC8450716). The database is built with a hierarchical LCA-merging
strategy: k-mers shared between sibling species are promoted to their common ancestor node
(e.g. genus) and not assigned to any individual species. Species nodes retain only k-mers
unique to that species within the genus. This means species-level detections in STAT output
represent genuinely species-diagnostic signal — uncontrolled k-mer cross-mapping between
sibling species does not occur at the database level.

Implication for `same_genus_secondary`: the original concern was k-mer bleed (reads from
the primary leaking into a same-genus secondary count). STAT's LCA design means this is
largely prevented — shared k-mers go to the genus node, not the species. However, same-genus
secondaries are still lower confidence for two other reasons:
- The diagnostic k-mer set for a species shrinks as genus-level k-mers are promoted;
  closely related species have fewer unique k-mers and are harder to distinguish
- The `specific_hits()` leaf algorithm uses count nesting to find leaf taxa, which can
  misattribute genus-level signal to one sibling when another is absent from the sample

Analysis of MAL data showed:
- 65% of secondary detections are same-genus — concentrated in Puccinia rusts
- Same-genus secondaries cluster near the detection threshold (<1%), consistent with low
  unique-k-mer counts rather than uncontrolled bleed
- Cross-kingdom detections are biologically immune to ambiguity problems

**High-confidence co-infections**: filter to `same_genus_secondary == False` (diff-genus
secondary). These account for ~31% of co-infected runs. The same-genus set (~69%) may also
contain real co-infections but requires more careful biological interpretation.

## Entrez / external API credentials
- `NCBI_API_KEY` env var → 10 req/s (9 used for Entrez); without: 2.5 req/s
- `S2_API_KEY` env var → Semantic Scholar 10 req/s; without: 1 req/s
  Apply free at semanticscholar.org/product/api (Academic Graph API)
- STAT endpoint rate: empirically ~2–10 req/s; RATE=15 + 32 workers configured
- EBI (EuropePMC, ENA): shared 2 req/s rate limiter (_ext_wait, gap=0.5s)
- http_get accepts no_retry_429=True — fail fast on 429 for external APIs
- User-Agent: `crypt/01_fetch_runs (leon.lenzo@curtin.edu.au)` etc.

## Output / cache layout

Each step directory has `data/` (outputs, caches) and `logs/` subdirectories.
Scripts create them automatically on first run.
**Logs are timestamped**: each run creates `logs/YYYY-MM-DD_HH-MM-SS/` subdir.
A `logs/latest` symlink always points to the most recent run.
Implemented via `_util.make_log_dir(OUT_DIR / "logs")` in all pipeline scripts.

```
output/
├── 00_build/
│   ├── data/
│   │   ├── phi-base_current.csv   PHI-base CSV (fetched from GitHub)
│   │   ├── ictv_vmr.xlsx          ICTV VMR Excel (fetched from ictv.global/vmr/current)
│   │   └── phibase_db.json        reference database
│   └── logs/
│       ├── build.log
│       └── _summary.txt
│
├── 01_fetch_runs/
│   ├── data/
│   │   ├── {mode}_runs.json       RunInfo rows keyed by Run accession
│   │   ├── {mode}_uids.json       fetched UID set (for resumability)
│   │   ├── stat_cache.jsonl       STAT responses — append-only, shared MAL+HAL
│   │   │                          format: accession<TAB>json_data per line
│   │   └── stat_cache_index.txt   accessions only — read at startup to resume
│   └── logs/
│       ├── history/YYYY-MM-DD_HH-MM-SS/  timestamped run dirs
│       │   ├── {mode}.log
│       │   └── {mode}_summary.txt
│       ├── {mode}.log.latest   (symlink → history/.../mode.log)
│       └── {mode}_summary.latest  (symlink → history/.../mode_summary.txt)
│
├── 02_filter_runs/
│   ├── data/
│   │   ├── runs.tsv               unified MAL+HAL output (13,323 rows, mode column)
│   │   ├── {mode}_confirmed.json  gate pass results (for inspection)
│   │   ├── {mode}_kingdom_dist.tsv per-run kingdom pcts (R-ready)
│   │   └── {mode}_species_dist.tsv per-detection species pcts (R-ready)
│   └── logs/
│       ├── history/YYYY-MM-DD_HH-MM-SS/
│       │   ├── filter.log
│       │   ├── filter_summary.txt
│       │   └── {mode}_validate_summary.txt
│       ├── filter.log.latest
│       └── filter_summary.latest
│
├── 03_fetch_meta/                 legacy — find_cache.json bootstraps 03b
│   └── data/
│       └── find_cache.json        v5 cache; read by 03b_fetch_literature.py on first run
│
├── 03a_fetch_xml/
│   ├── data/
│   │   ├── bioprojects.json       title, description, submission_date (1,754 entries)
│   │   └── biosamples.json        BioSample harmonised attrs (11,117 entries)
│   │                              fields: tissue, geo_loc_name, collection_date, lat_lon,
│   │                              host, isolation_source, dev_stage
│   └── logs/
│
├── 03b_fetch_literature/
│   ├── data/
│   │   └── literature.json        PMID + PMC methods per BioProject (resumable)
│   │                              fields: primary_pmid, primary_pub_date,
│   │                              primary_publication, abstract, methods_text,
│   │                              pmcid, n_papers_found, pmid_source
│   └── logs/
│
├── 04_filter_kw/
│   ├── data/
│   │   └── biosample_kw.tsv       one row per biosample_representative BioSample
│   │                              (11,117 rows); BioProject metadata as repeated cols
│   │                              See "Output schema (04_filter_kw.py)" above.
│   └── logs/
│
├── 05_llm_classify/
│   ├── data/
│   │   ├── classify_cache.jsonl   append-only LLM cache (one line per BioProject)
│   │   └── bioproject_llm.tsv     one row per BioProject (1,754 rows)
│   │                              cols: BioProject, llm_treatment, llm_study_setting,
│   │                              llm_named_pathogen, llm_named_host, llm_tissue,
│   │                              llm_kw_treatment_agree, llm_kw_setting_agree,
│   │                              llm_confidence, llm_rationale
│   └── logs/
│
├── analysis/
│   └── primary_alignment.tsv      three-way alignment check (biosample_rep rows)
│
├── legacy/                        old pipeline outputs (03_validate, 04_crypt, 05_meta)
│
└── figures/
    ├── host_tree/      crypt_host_tree.py + .R  →  .nwk, _meta.tsv, .pdf, .png
    ├── guilds/         mal_guilds.py + .R       →  guild nodes/edges.tsv, network.pdf/png
    ├── scatter/        scatter.R                →  scatter.pdf/png
    ├── coinf_rate/     coinf_rate.R             →  coinf_rate.pdf/png
    ├── novel_heatmap/  novel_heatmap.R          →  novel_heatmap.pdf/png
    ├── kingdom_comp/   kingdom_comp.R           →  kingdom_comp.pdf/png
    ├── sankey/         sample_funnel.py         →  sample_funnel.png/.html
    └── study_design/   study_design_figures.R   →  .pdf/png
```

## Actual run results (2026-07-30)

**01_fetch_runs.py / 02_filter_runs.py:**

| Mode | Runs fetched | Non-RNA excluded | Screened | Gate pass | Gate pass rate |
|------|-------------|-----------------|---------|-----------|----------------|
| MAL  | 48,418      | 1,419           | 46,315  | 6,191     | 13.4%          |
| HAL  | 559,950     | 7,343           | 546,816 | 7,118     | 1.3%           |

- output/02_filter_runs/data/runs.tsv — 13,323 rows; BioProjects: 1,754 total; 401 co-infected
- 608,368 STAT cache entries total

**03a_fetch_xml.py:**
- 11,117 BioSamples in biosamples.json; BioSample XML coverage: geo_loc_name 61%, tissue 56%, collection_date 42%

**04_filter_kw.py:**
- output/04_filter_kw/data/biosample_kw.tsv — 11,117 rows
- 63.8% co-infected; 55.2% high-confidence (same_genus_secondary=False); 67.1% geolocated
- 90.3% named_host resolved; 15% from biosample XML, 75% from title/abstract text
- study_setting: field=4,386, lab=5,873, unclear=858

**05_llm_classify.py (OpenAI gpt-4o-mini, 2026-07-30):**
- 1,754/1,754 BPs classified; 0 failures
- llm_treatment: single=747, host_study=403, surveillance=293, abiotic_stress=197, coinf_experiment=100, unclear=14
- KW agreement: treatment 60.1%, setting 92.7%

**Figures (output/figures/{subdir}/):**
- scatter/scatter.R: host% vs pathogen%, coloured single/coinf/novel, shape=MAL/HAL
- host_tree/crypt_host_tree.py + .R: fan tree 173 hosts, bars = n_single/n_multi
- guilds/mal_guilds.py + .R: co-occurrence network (277 nodes / 254 edges, MAL+HAL)
- coinf_rate/coinf_rate.R: top 20 hosts by co-infection count, paired bars log scale
- novel_heatmap/novel_heatmap.R: top 35 pathogens × top 30 hosts (590 novel BioSamples)
- kingdom_comp/kingdom_comp.R: stacked bar co-infection flag + secondary kingdom by mode
- study_design/treat_setting_bar.R + single_field_lab.R: field vs lab field signal figures

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
- `_genus(name)` returns the first word (lowercase) — used to set `same_genus_secondary`
- MAL secondary detection excludes the primary/library organism by comparing
  PHI-base seed taxids (so strains of the same species are also excluded)
- `txid33090[Host]` returns 0 results in SRA — [Host] is a free-text field.
  In-planta confirmation is done via STAT Viridiplantae reads (host_pct).
- **MAL `host` column fix (2026-07-29)**: `detect_host_species()` now uses a 305,355-name
  Viridiplantae allowlist (`viridiplantae_names` in phibase_db.json, built by 00_build.py step 7/7)
  to filter STAT hits. Only names confirmed as Viridiplantae are accepted as plant hosts;
  non-plant organisms (co-infecting fungi, bacteria, etc.) are rejected. Crypt.tsv re-run
  shows 73.2% species-level resolution, 5.5% tribe-level (e.g. Triticinae), 13.9% broad-clade,
  7.4% unresolved ("Viridiplantae"). Broad-clade and unresolved runs are excluded from guild
  network figures via `BROAD_CLADE_NAMES` in `output/figures/guilds/field_hc_guilds.py`.
  `field_hc_nodes.tsv` now includes a `top_host` column (most frequent resolved host per
  pathogen node). `field_hc_guilds.py` also loads `named_host` from bioproject_meta.tsv and
  uses it as a fallback when STAT host is in BROAD_CLADE_NAMES (82 runs recovered, only 3
  excluded after fallback, vs 85 excluded before). `wheat_cluster.R` draws overlapping hulls
  from `field_hc_node_hosts.tsv` (one row per pathogen×host pair). Hulls filtered by
  `HULL_WHITELIST` (explicit biologically credible host list) and `MIN_HOST_NODES=3` (host
  must span ≥3 cluster nodes). `Aegilops tauschii` remapped to `Triticum aestivum` in the R
  hull layer (D-genome donor; STAT k-mers are indistinguishable from wheat D-subgenome reads).
- **Fusarium blind spot in HC network**: 94% of Fusarium co-infections are same-genus (multiple
  Fusarium species co-detected), so nearly all are removed by the `same_genus_secondary=False`
  HC filter. Only 6 Fusarium appearances in field HC (2 as primary, 4 as secondary). This is
  a real limitation of the HC filter for pathogens where intra-genus co-infection is the norm.
  Consider showing Fusarium in a separate same-genus-allowed view for the paper.
- 00_build.py NCBITaxa SQLite schema: `species(taxid, parent, spname, common, rank, track)`
  and `synonym(taxid, spname)` — no `name_class` column (ete3-specific schema)
- `_expand_taxids()` MUST use `intermediate_nodes=True` — default (`False`) returns only
  leaf taxa (no children), silently dropping species-level nodes that have named strains.
  This caused PVY, ToBRFV, CMV, TuMV, TuYV, f.sp. taxa etc. to be missing from
  `name_to_taxid`, breaking the HAL gate for virus-dominated runs (2026-07-25 bug fix).
- `get_descendant_taxa` uses a preorder/postorder traversal cache (`.traverse.pkl`).
  Leaves appear with count==1 in the traversal window; internal nodes with count==2.
  `intermediate_nodes=False` returns only count==1 (leaves), excluding species nodes
  that have children (strains). `intermediate_nodes=True` returns all counts.
- ICTV VMR URL `https://ictv.global/vmr/current` redirects to the versioned file
  (e.g. VMR_MSL41.v1.20260721.xlsx) — the redirect target is stable per release
- Plant viruses in ICTV include insect-vectored species (geminiviruses, tospoviruses)
  via `Host source = "invertebrates, plants"` — included intentionally as they infect
  plants in the field even though the insect is the vector
- `interaction_status == "novel_host_range"` is the key signal; `llm_treatment == "coinf_experiment"`
  in bioproject_llm.tsv flags intentional co-infection designs for manual exclusion
- 04_filter_kw.py two-axis keyword output: `treatment` × `study_setting`.
  treatment priority: coinf_experiment > abiotic_stress > host_study > single > unclear.
  Lab wins over field when both keyword sets fire. geo_loc_name in BioSample XML upgrades
  "unclear" → "field" when no keywords fire. Keyword provenance in
  `treatment_keywords` + `setting_keywords` columns.
- kw treatment categories (04_filter_kw.py keyword counts):
    coinf_experiment  intentional co-infection
    abiotic_stress    drought/heat/salt/cold experiment — no declared pathogen
    host_study        pure host biology, development/assembly/physiology
    single            single-pathogen or host-response study
    unclear           no text available
  NOTE: ~51% of kw host_study BPs also contain pathogen/disease language — keyword contamination.
  Use llm_treatment from bioproject_llm.tsv (05_llm_classify.py) for higher-accuracy classification.
- LLM treatment categories (05_llm_classify.py, gpt-4o-mini):
    single=747, host_study=403, surveillance=293, abiotic_stress=197,
    coinf_experiment=100, unclear=14. KW agreement: treatment 60.1%, setting 92.7%.
- 04_filter_kw.py `--hc` restricts to same_genus_secondary=False; default includes all
- `scripts/review_designs.py` generates an HTML annotation tool; loads biosample_kw.tsv +
  bioproject_llm.tsv; re-run after 04_filter_kw.py or 05_llm_classify.py to refresh
- 03b_fetch_literature.py PMID search uses 6 strategies (short-circuit, stop at first hit):
  BioProject XML → PMC full-text → Europe PMC → ENA XML → Semantic Scholar → PubMed.
  OpenAlex DROPPED: moved to credit-based model (10 credits/req, 100 req/18hr free tier).
  elink DROPPED: 10-30s per call; benchmark showed 0% discovery yield on unlabelled entries.
  Semantic Scholar: S2_API_KEY required (1 req/s); skips entirely without key (was getting
  silent 429s that looked like genuine misses). Key obtained 2026-07-27.
  S2 rate limit is 1 req/s with or without key; _s2_wait() uses 1.1s gap regardless.
  Primary paper = earliest pub date. bp_submission_date from BioProject XML attr.
  Benchmark in scripts/diagnostics/benchmark_strategies.py — measures discovery yield on UNLABELLED
  entries (no ground truth); 20-entry sample showed PMC 40%, EuropePMC 35%, S2 0% (429s).
- 03b_fetch_literature.py fetches PMC full-text methods section (up to 8,000 chars) per paper.
  `_pmid_to_pmcid()` searches PMC with `{pmid}[pmid]`; `_fetch_pmc_methods()` fetches JATS XML
  and extracts the first section with "method" or "material" in sec-type or title.
  Stored in literature.json as `methods_text` and `pmcid`. Used by 04_filter_kw._classify().
- **`_pmcids_to_pmids` bug (fixed in find_cache.json v4)**: original used `elink dbfrom=pmc db=pubmed`
  which by default returns references CITED BY the PMC articles (pmc_refs_pubmed link),
  not the articles themselves. This caused tool papers like BLAST (PMID 2231712, 1990) to
  appear as "earliest publication" because they are widely cited. Fixed by using
  `esummary db=pmc` which returns article metadata including the article's own PMID
  via the `articleids` field.
- find_cache.json v5 adds `pmcid` and `methods_text` per entry (v4 → v5 upgrade re-fetches
  PMC methods for entries that have PMIDs). 03b_fetch_literature.py reads find_cache.json
  as bootstrap and writes to literature.json (separate file; no schema version needed).
- 01_fetch_runs.py STAT fetch hangs overnight: all 32 ThreadPoolExecutor threads block on network I/O;
  SIGINT does not exit cleanly — requires SIGKILL. Happened twice during HAL fetch
  (~13 min and ~70 min idle). Watchdog/auto-restart flagged for future rebuild.

## Active task list (2026-07-30)

Completed in prior sessions (kept for context):

- [x] Two-axis keyword classification (`treatment` × `study_setting`) in `04_filter_kw.py`
- [x] LibrarySource pre-filter (1,419 MAL + 7,343 HAL GENOMIC runs excluded)
- [x] HTML annotation tool (`scripts/review_designs.py`)
- [x] Primary organism alignment check (`scripts/check_primary_alignment.py`)
- [x] MAL host detection fix + broad-clade filter (Viridiplantae allowlist in 02_filter_runs.py)
- [x] Study design figures (`output/figures/study_design/`): field 2× rate vs lab; hc filter collapses lab
- [x] **Pipeline restructure (2026-07-30)** — 7-script architecture:
       03_fetch_meta.py → split into 03a_fetch_xml.py (BioSample XML) + 03b_fetch_literature.py
       04_filter_meta.py → renamed 04_filter_kw.py; output biosample_kw.tsv (BioSample rows)
       05_llm_classify.py — new; OpenAI gpt-4o-mini; 1,754 BPs; output bioproject_llm.tsv
       BioSample XML fixed HTTP 414 bug (POST efetch); 11,117 entries; geo 61%, tissue 56%
       Named host resolved: 90% coverage; BioSample XML attr priority over title/abstract text

Pending (priority order):

1. [x] **Pipeline cleanup** — superseded scripts moved to scripts/legacy/; output/legacy/ deleted;
       figure script column names already updated. Remaining: review_designs.py still loads
       output/04_filter_meta/data/bioproject_meta.tsv (stale path + stale data model — needs
       rewrite to load biosample_kw.tsv + bioproject_llm.tsv; deferred pending user decision).

2. [ ] **Kraken2 replacement for STAT** — STAT shown unreliable for some PST races (2026-08-04
       kallisto pilot: kraken/kallisto_pilot/). Plan: Kraken2 database of all PHI-base euk pathogen
       genomes; screen all ~593k step-01 runs on Setonix; STAT vs Kraken comparison = paper
       figure. Assembly coverage re-check underway (scripts/refseq_coverage_all.log).
       See memory/kraken_pipeline.md for full design.

3. [ ] **same_family_secondary flag** — add to runs.tsv analogous to same_genus_secondary.
       Precompute seed taxid → family mapping in 00_build.py via ete3; store in phibase_db.json;
       check in 02_filter_runs.py at annotation time. Gives three confidence tiers:
         same_genus=True       → lowest confidence (k-mer ambiguity)
         same_family=True      → flag for potential realignment validation
         both False            → highest confidence cross-family diff-genus detection
       Botrytis/Sclerotinia (Sclerotiniaceae) are legitimate co-pathogens — retain regardless.
       Fusarium/Alternaria as saprophytic co-pathogens are valid detections regardless of family.

3. [ ] **Submission timeline + PMID figure** — `bp_submission_date` on x-axis, bars/points
       coloured by PMID available vs missing. New script `output/figures/timeline/`.

4. [ ] **Systematic no-PMID BioProject review** — ~886 no-PMID BioProjects.
       Also explore CrossRef title search for additional PMIDs.

5. [ ] **Tissue normalization** — heavy case variation in BioSample XML tissue values
       (Leaf/leaf/leaves/Leaves as separate entries). Normalise in a post-processing pass.

6. [ ] **Pipeline figures for README** — flowchart for each step. `output/figures/pipeline/`.

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

## Architecture decisions (resolved)

- 01_sra + 02_stat → merged into 01_fetch_runs.py ✓
- MAL + HAL outputs → unified crypt.tsv with `mode` column ✓
- Analysis scripts → `output/figures/{subdir}/` (each chart in its own subdir) ✓
- 03_fetch_meta.py scope → ALL BioProjects (not just co-infected) for prevalence analysis ✓
- Logs → timestamped subdirs under logs/history/ via _util.make_log_dir ✓
- Logs → .latest symlinks at logs/ root via _util.link_latest ✓
- 03_fetch_meta.py short-circuit PMID search (stop at first hit) ✓
- OpenAlex dropped (credit-based); Semantic Scholar added with S2_API_KEY support ✓
- elink + _openalex_search dead code removed from 03_fetch_meta.py ✓
- Legacy output/02_stat/ + output/01_sra/ fallbacks removed from 02_filter_runs.py ✓
- benchmark_strategies.py moved to scripts/; rewritten to test discovery yield on unlabelled entries ✓
- S2_API_KEY obtained 2026-07-27; rate fixed to 1.1s gap (1 req/s limit) ✓
- 01_fetch_runs.py watchdog/auto-restart → deferred (workaround: tmux + manual SIGKILL+restart)
- 03_fetch_meta.py + 04_filter_meta.py split (2026-07-28): mirrors 01/02 fetch/filter pattern.
  03 is pure fetcher (no TSV output); 04 is pure analysis (no network calls). ✓
- 03_fetch_meta.py → split into 03a_fetch_xml.py (BioProject/BioSample XML) +
  03b_fetch_literature.py (PMID + PMC methods). find_cache.json bootstraps 03b on first run. ✓
- 04_filter_meta.py → renamed 04_filter_kw.py; output changed from bioproject_meta.tsv
  (one row per BioProject) to biosample_kw.tsv (one row per BioSample). BioSample is the
  fundamental unit of analysis; BioProject metadata repeated as columns. ✓
- LLM classification moved to standalone 05_llm_classify.py (OpenAI gpt-4o-mini);
  output bioproject_llm.tsv joins biosample_kw.tsv on BioProject. Reproducible, cacheable. ✓
- BioSample HTTP 414 fix (2026-07-30): efetch POST with 300 IDs per batch; stubs re-fetched. ✓

## Target journal

New Phytologist (IF ~10) or ISME Journal (IF ~11). Frame as biology discovery
(novel host-pathogen interactions + phylogenetic host susceptibility structure),
not as a methodology paper.
