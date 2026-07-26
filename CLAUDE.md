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

Legacy scripts (01_sra.py → 05_meta.py) are preserved in `legacy/` for reference.
The canonical pipeline is the four scripts below.

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

01_fetch.py   Fetch SRA run IDs + NCBI STAT taxonomy for MAL or HAL.
              Merges legacy 01_sra.py + 02_stat.py into a single data-gathering step.
              No retention gate — gate moved to 02_filter.py.
              MAL: batches 205 seed pathogen taxids; HAL: 180 seed host taxids.
              [Organism:exp] lets SRA expand to all descendant strains.
              Parallel fetch: SRA_MAX_WORKERS=8 (RunInfo) + STAT_MAX_WORKERS=32 (STAT).
              Separate rate locks: _ri_lock (SRA) + _stat_lock (STAT).
              CACHE_MIN_PCT=0.1 filter applied at STAT write time.
              Output: output/01_fetch/data/{mode}_runs.json + stat_cache.jsonl
              Falls back to output/02_stat/ cache if output/01_fetch/ absent.

02_filter.py  Three phases in one script — gate → validate → crypt.
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
              Stat cache fallback: output/01_fetch/ → output/02_stat/data/ → output/02_stat/
              Runs fallback: output/01_fetch/ → output/01_sra/data/ → output/01_sra/
              KINGDOM_THRESHOLDS in this file — edit then --skip-validate to re-run crypt.
              Unified output: output/02_filter/data/crypt.tsv with `mode` column.
              BioSample dedup applied across both modes in the combined output.

03_find.py    Fetch BioProject metadata for co-infection BioProjects.
              Reads output/02_filter/data/crypt.tsv (no --mode flag).
              Adds `modes` column (mal / hal / mal+hal) so cross-mode BioProjects visible.
              PMIDs collected from seven sources; earliest pub = depositing paper:
                1. BioProject XML <Publication> elements (submitter-supplied, zero
                   extra API calls — parsed from the efetch already made)
                2. elink bioproject→pubmed
                3. PMC full-text search for accession string
                4. PubMed text search for accession string
                (sources 5–7 only run when 1–4 find nothing)
                5. Europe PMC full-text search (preprints, supplementary-method mentions)
                6. ENA XML API for PRJEB accessions (submitter-supplied PUBMED_ID)
                7. OpenAlex full-text search (250M+ works, different indexing from NCBI/EBI)
              Cache logic: entries with PMIDs → cache hit (no re-fetch); entries without
              PMIDs → retry sources 4+5 on each run; still nothing → return cached blank.
              study_design: coinf_experiment / field_survey / unclear (keyword scan of
              title + description + abstract). coinf_experiment → manual review required.
              --hc flag restricts to same_genus_secondary=False.
              Output: output/03_find/data/bioproject_meta.tsv
```

### Running the pipeline

```bash
# Step 0: build reference DB (system python3 required)
python3 00_build.py

# Step 1: fetch SRA runs + STAT (one mode at a time — shared stat cache)
python 01_fetch.py --mode mal
python 01_fetch.py --mode hal

# Step 2: filter + detect (both modes by default, unified output)
python 02_filter.py                     # validate + gate + crypt, both modes
python 02_filter.py --skip-validate     # after manual threshold review
python 02_filter.py --mode mal          # single mode

# Step 3: BioProject metadata
python 03_find.py
python 03_find.py --hc    # restrict to same_genus_secondary=False
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
`PhibaseDB` (in 02_filter.py) derives `pathogen_taxids`, `host_taxids`,
`host_to_pathogens` at load time from the kingdom dicts — they are NOT stored in JSON.

## STAT approach

- **Parallel fetch** with global timestamp rate-lock (not per-worker sleep):
  `_rate_lock` + `_rate_last` ensures ≤ RATE requests/s total across all threads
- **STAT endpoint** (`trace.ncbi.nlm.nih.gov`) is separate from Entrez — tested to
  handle ~2–10 req/s (highly variable by time of day); RATE=15, MAX_WORKERS=32
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

## STAT detection thresholds (KINGDOM_THRESHOLDS in 02_filter.py)

| Kingdom  | Detection threshold |
|----------|---------------------|
| Fungi    | ≥ 1.0%              |
| Viruses  | ≥ 10.0%             |
| Bacteria | ≥ 5.0%              |
| Oomycota | ≥ 0.5%              |
| Nematoda | ≥ 1.0%              |

MAL in-planta gate: `MAL_MIN_HOST_PCT = 1.0` (Viridiplantae reads ≥ 1%).
HAL pathogen gate:  `HAL_MIN_PATHOGEN_PCT = 1.0` (any PHI-base pathogen/virus ≥ 1%).

## Output schema (02_filter.py → `crypt.tsv`)

| Column | MAL | HAL |
|---|---|---|
| `Run` | SRR accession | SRR accession |
| `BioSample` | NCBI BioSample accession | same |
| `host` | top Viridiplantae leaf in STAT | library organism (ScientificName) |
| `host_pct` | Viridiplantae % in STAT | Viridiplantae % in STAT |
| `host_species_pct` | STAT leaf species % | — |
| `primary_pathogen` | library organism (ScientificName) | top PHI-base/ICTV pathogen in STAT |
| `primary_taxid` | resolved taxid (or blank) | resolved taxid |
| `primary_pct` | — (library dominates) | % of analysed reads |
| `secondary_pathogens` | other PHI-base/ICTV pathogens in STAT | remaining pathogens |
| `secondary_kingdoms` | kingdoms of secondaries (e.g. Viruses; Fungi) | same |
| `n_secondary` | count | count |
| `interaction_status` | known / novel_host_range / novel_combination / unresolved | same |
| `co_infection_flag` | single / multi_species / multi_kingdom | same |
| `same_genus_secondary` | True if any secondary shares genus with primary | same |
| `biosample_n_runs` | confirmed runs sharing this BioSample | same |
| `biosample_representative` | True for one run per BioSample (highest analyzed) | same |
| `fungi_pct` etc. | kingdom % from STAT | same |
| `analyzed` | total reads analysed by STAT | same |

Plus RunInfo passthrough: `BioProject`, `SRAStudy`, `Platform`, `ScientificName`.

`co_infection_flag` priority: `multi_kingdom` > `multi_species` > `single`.
`interaction_status`: replaces old `known_interaction` bool.
  `novel_host_range` is the primary signal of interest — pathogen known to PHI-base
  but not recorded on this host species. `unresolved` means taxid lookup failed;
  cannot assess novelty. For plant viruses, `known` is returned for any PHI-base
  plant host (ICTV viruses are mapped to all 180 host seeds in DB build).
`same_genus_secondary`: True when the first word of any secondary name matches the
  first word of the primary. Use `same_genus_secondary == False` for highest-confidence
  co-detections. STAT's LCA k-mer design means species-level signals are genuinely
  diagnostic, but same-genus pairs have smaller unique-k-mer sets.
`biosample_representative`: filter to True to get one row per biological sample —
  use this for paper-level co-infection rate calculations.

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

## Entrez credentials
- API key: `NCBI_API_KEY` env var (set in `~/.bashrc`) → 10 req/s (9 used for Entrez)
- Without key: 2.5 req/s
- STAT endpoint rate: empirically ~2–10 req/s (variable); RATE=15 + 32 workers configured
- User-Agent: `crypt/01_sra (leon.lenzo@curtin.edu.au)` etc.

## Output / cache layout

Each step directory has `data/` (outputs, caches) and `logs/` subdirectories.
Scripts create them automatically on first run.

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
├── 01_fetch/
│   ├── data/
│   │   ├── {mode}_runs.json       RunInfo rows keyed by Run accession
│   │   ├── {mode}_uids.json       fetched UID set (for resumability)
│   │   ├── stat_cache.jsonl       STAT responses — append-only, shared MAL+HAL
│   │   │                          format: accession<TAB>json_data per line
│   │   └── stat_cache_index.txt   accessions only — read at startup to resume
│   └── logs/
│       ├── {mode}.log
│       └── {mode}_summary.txt
│
├── 02_stat/                       stub — 3 symlinks into 01_fetch/ for HAL process
│   ├── stat_cache.jsonl        -> ../01_fetch/data/stat_cache.jsonl
│   ├── stat_cache_index.txt    -> ../01_fetch/data/stat_cache_index.txt
│   └── hal.log                 -> ../01_fetch/logs/hal.log
│   (remove after HAL finishes)
│
├── 02_filter/
│   ├── data/
│   │   ├── crypt.tsv              unified MAL+HAL co-infection table (mode column)
│   │   ├── {mode}_confirmed.json  gate pass results (for inspection)
│   │   ├── {mode}_kingdom_dist.tsv per-run kingdom pcts (R-ready)
│   │   └── {mode}_species_dist.tsv per-detection species pcts (R-ready)
│   └── logs/
│       ├── filter.log
│       ├── filter_summary.txt
│       └── {mode}_validate_summary.txt
│
├── 03_find/
│   ├── data/
│   │   ├── bioproject_meta.tsv    one row per BioProject; modes col = mal/hal/mal+hal
│   │   └── meta_cache.json        Entrez+EBI cache (v4; resumable)
│   └── logs/
│       ├── meta.log
│       └── meta_summary.txt
│
├── legacy/                        old pipeline outputs (03_validate, 04_crypt, 05_meta)
│
└── figure/
    ├── crypt_host_tree.nwk
    └── crypt_host_tree_meta.tsv
```

```
figure/
├── crypt_host_tree.py    Build nwk + meta from crypt.tsv (system python3)
└── crypt_host_tree.R     Fan tree figure — green/orange bars per host species
                          Run from crypt/: Rscript figure/crypt_host_tree.R
```

### After HAL stat fetch completes

```bash
# Remove the 02_stat symlink stub (HAL process will have exited)
rm -r output/02_stat/

# Re-run filter + find on full MAL+HAL data
python 02_filter.py          # both modes, unified crypt.tsv
python 03_find.py            # BioProject metadata with all 6 sources
```

## Actual run results (2026-07)

| Mode | UIDs | Runs fetched | Gate pass | Gate pass rate |
|------|------|-------------|-----------|----------------|
| MAL  | 46,540 | 48,418    | 6,852     | 14.2%          |
| HAL  | —      | 559,328   | in progress (partial: 1,392 / 63k screened = 2.2%) | — |

**MAL co-infection summary** (full dataset):
- 6,852 classified → 1,148 co-infected (16.8%)
- same_genus_secondary=True: 788 (68.6% of co-infected) — k-mer bleed risk
- same_genus_secondary=False (high-confidence): **360 runs / 54 BioProjects**
- Top HC secondaries: Pepino mosaic virus, Zymoseptoria tritici, Alternaria alternata,
  Botrytis cinerea, Varicosavirus lactucae, Parastagonospora nodorum

**HAL co-infection summary** (partial — 63k/559k screened):
- 793 classified → 98 co-infected (12.4%)
- same_genus_secondary=False (HC): **32 runs / 21 BioProjects**
- Top HC secondaries: Peanut stripe virus, Bipolaris zeicola, Potyvirus phaseovulgaris,
  Pectobacterium brasiliense, Capillovirus mali, Grapevine berry inner necrosis virus
- Known interaction rate much higher than MAL (79.6% vs 37.0%) — host is precisely
  identified as library organism in HAL, improving PHI-base matching

## Performance notes

- **01_fetch SRA RunInfo fetch**: parallel POST requests (ThreadPoolExecutor, MAX_WORKERS=8).
  GET requests with 500 UIDs cause HTTP 414; POST avoids this. ~20 req/s achieved.
- **01_fetch STAT fetch**: STAT endpoint (trace.ncbi.nlm.nih.gov) is highly variable —
  observed ~2 req/s overnight to ~10 req/s during daytime.
  MAL (~48k): ~1.5 hrs. HAL (~559k): estimated 3–4 days depending on server load.
- **Shared stat_cache**: do NOT run MAL and HAL 01_fetch simultaneously — last writer
  wins on cache saves. Chain HAL after MAL:
  `until grep -q '01_fetch MAL summary' output/01_fetch/logs/mal.log; do sleep 30; done`

## Notes

- `specific_hits()` leaf detection in 02_filter.py is the core non-trivial algorithm
- `_genus(name)` returns the first word (lowercase) — used to set `same_genus_secondary`
- MAL secondary detection excludes the primary/library organism by comparing
  PHI-base seed taxids (so strains of the same species are also excluded)
- `txid33090[Host]` returns 0 results in SRA — [Host] is a free-text field.
  In-planta confirmation is done via STAT Viridiplantae reads (host_pct).
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
- `interaction_status == "novel_host_range"` is the key signal; `coinf_experiment`
  in 03_find study_design flags intentional designs for manual exclusion
- 03_find `--hc` restricts to same_genus_secondary=False; default includes all co-infected
- 03_find PMID search uses 7 sources: BioProject XML `<Publication>`, elink bioproject→pubmed,
  PMC full-text esearch, PubMed text esearch, Europe PMC, ENA XML API, OpenAlex.
  Primary paper = earliest pub date. Improved MAL coverage from 8/54 (elink only) to 64/118.
- **`_pmcids_to_pmids` bug (fixed in cache v4)**: original used `elink dbfrom=pmc db=pubmed`
  which by default returns references CITED BY the PMC articles (pmc_refs_pubmed link),
  not the articles themselves. This caused tool papers like BLAST (PMID 2231712, 1990) to
  appear as "earliest publication" because they are widely cited. Fixed by using
  `esummary db=pmc` which returns article metadata including the article's own PMID
  via the `articleids` field. Cache version bumped to `_v: 4` in 03_find.py.
- 01_fetch.py STAT fetch hangs overnight: all 32 ThreadPoolExecutor threads block on network I/O;
  SIGINT does not exit cleanly — requires SIGKILL. Happened twice during HAL fetch
  (~13 min and ~70 min idle). Watchdog/auto-restart flagged for future rebuild.

## Deferred analysis ideas (after HAL pipeline complete)

1. **Co-infection guilds** (in progress — `figure/mal_guilds.py` + `figure/mal_guilds.R`)
   Unipartite co-occurrence network: nodes = pathogens, edges = shared confirmed run,
   weight = run count. Already shows wheat rust cluster, tomato virus pair (PepMV↔ToBRFV),
   Botrytis↔Sclerotinia necrotrophic guild. Will improve with HAL data.

2. **Host susceptibility landscape**
   Co-infection rate per host species on `crypt_host_tree` figure. Are some hosts
   systematically "dirtier"? Is susceptibility phylogenetically conserved?

3. **Gate failure characterisation**
   What are the 85.8% of MAL runs failing Viridiplantae ≥1%? Stream stat_cache.jsonl
   for failed accessions — kingdom composition may reveal pure-culture contamination
   vs. low-host-content field samples.

4. **Viral threshold sensitivity**
   Quick re-run of 04_crypt.py at 5% vs 10% virus threshold to quantify impact.
   Method validation / numbers for the paper.

## Planned architecture review (next session priority)

Review script structure before adding more analysis. Key questions:
- Should 01_sra + 02_stat merge into a single fetch pipeline?
- Should MAL + HAL outputs merge into one TSV with a `mode` column?
- Where should analysis scripts live (new `06_analysis/` step, or `figure/`)?
- 02_stat.py watchdog/auto-restart design (for rebuild as 02_fetch.py)
