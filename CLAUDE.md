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

All four scripts are **fully standalone** — no shared package imports, stdlib + json only
(except ete3 and openpyxl in 00_build.py which requires system python3).

```
00_build.py   Build PHI-base + ICTV reference DB.
              MUST run with system python3 (ete3/sqlite3 incompatibility
              with miniconda).
              Auto-downloads PHI-base CSV from GitHub if absent.
              Auto-downloads ICTV VMR Excel from ictv.global if absent.
              Plant viruses: Host source in {plants, plants (S),
              invertebrates, plants} → resolved via ete3 to NCBI taxids.

01_sra.py     Fetch SRA run IDs for MAL or HAL.
              MAL: batches 205 seed pathogen taxids, 50/query,
                   txid{t}[Organism:exp] AND RNA-Seq[Strategy]
              HAL: batches 180 seed host taxids, 50/query,
                   txid{t}[Organism:exp] AND RNA-Seq[Strategy]
              [Organism:exp] lets SRA expand to all descendant strains.
              RunInfo fetched via HTTP POST (avoids 414 on large UID batches).
              Parallel fetch: ThreadPoolExecutor(MAX_WORKERS=8) + shared
              rate lock. Saves RunInfo + UID cache for resumability.

02_stat.py    Fetch NCBI STAT for all runs (parallel, rate-limited).
              STAT endpoint: trace.ncbi.nlm.nih.gov (different from Entrez).
              Empirically caps at ~10 req/s; configured MAX_WORKERS=32,
              RATE=15.0 to saturate without being artificially limited.
              Do NOT run MAL and HAL simultaneously — they share
              stat_cache.json and the last writer wins on saves.
              Applies retention gate:
                MAL: keep if Viridiplantae reads ≥ MAL_MIN_HOST_PCT (1%)
                HAL: keep if any PHI-base pathogen/virus ≥ HAL_MIN_PATHOGEN_PCT (1%)

03_crypt.py   Identify cryptic co-infections in confirmed mixed libraries.
              Uses specific_hits() leaf-detection algorithm to find the most
              specific organism name in the STAT taxonomy tree under each
              kingdom node, then cross-references against PHI-base + ICTV.
              Produces final TSV with co-infection classification.
```

## PHI-base + ICTV reference database (`output/00_build/phibase_db.json`)

Built by `00_build.py` using ete3 NCBITaxa (local SQLite at `~/.etetoolkit/taxa.sqlite`).

**PHI-base** (fungi, bacteria, oomycetes, nematodes):
- 205 pathogen seed species → ~7,985 expanded taxids (strains, subspecies, f. sp.)
- 180 host seed species → ~1,869 expanded taxids (cultivars, varieties)
- Interaction map: PHI-base host × pathogen pairs (species-level, bidirectional)

**ICTV VMR** (plant viruses):
- ~2,630 plant virus species (Host source: plants / invertebrates, plants)
- Resolved to NCBI taxids via `ncbi.get_name_translator()` — unresolved names skipped
- Expanded to include all descendant strains via `_expand_taxids()`
- Viruses mapped to all 180 PHI-base plant host seeds in `pathogen_to_hosts`
  → `known_interaction(virus, any_plant_host)` returns True
- ICTV species names added to `name_to_taxid` as supplement to NCBI names
- Auto-downloaded from `https://ictv.global/vmr/current` (redirects to latest release)

**Combined DB stats** (approximate, post-rebuild):
- ~11,000+ name entries — scientific names + synonyms + ICTV names, lowercase
- Seed maps — trace any expanded taxid → PHI-base/ICTV species seed
- `virus_taxids` key identifies which pathogen taxids are viruses

JSON structure (keys):
```
pathogen_taxids          list of all expanded pathogen taxids (PHI-base + ICTV viruses)
virus_taxids             subset of pathogen_taxids that are plant viruses (ICTV)
host_taxids              list of all expanded host taxids
contaminant_taxids       known sequencing contaminants (phiX174, E.coli K12, etc.)
name_to_taxid            {lowercase_name: taxid} — for resolving STAT org strings
taxid_to_name            {str(taxid): canonical_name}
pathogen_taxid_to_seed   {str(expanded_taxid): seed_taxid}
host_taxid_to_seed       {str(expanded_taxid): seed_taxid}
pathogen_to_hosts        {str(seed_pathogen_taxid): [seed_host_taxids]}
                         virus entries map to all 180 PHI-base plant host seeds
host_to_pathogens        {str(seed_host_taxid): [seed_pathogen_taxids]}
meta                     build date, counts, scope, vmr_path
```

## STAT approach

- **Parallel fetch** with global timestamp rate-lock (not per-worker sleep):
  `_rate_lock` + `_rate_last` ensures ≤ RATE requests/s total across all threads
- **STAT endpoint** (`trace.ncbi.nlm.nih.gov`) is separate from Entrez — tested to
  handle ~10 req/s; configured at RATE=15, MAX_WORKERS=32 to stay below ceiling
- **One-pass species resolution**: parse all kingdoms in a single pass per run
- **`specific_hits(table, node, analyzed)`** — leaf-level species detection:
  1. Find all entries with `total_count ≤ node_count` (i.e., nested under this kingdom)
  2. Detect "leaf" counts: a count is a leaf if no other count falls in `[LEAF_FRAC*c, c)`
  3. Prefer 2-word binomials over clade names or strain strings
  4. Returns `[(name, pct), ...]` sorted by pct descending
- **Mycoviruses excluded** from pathogen detection (keywords: mycovirus, mitovirus,
  hypovirus, etc. — they infect the fungus, not the plant)
- STAT shared cache at `output/02_stat/stat_cache.json` — keyed by Run accession

## STAT detection thresholds (hardcoded in 02_stat.py and 03_crypt.py)

| Kingdom  | Detection threshold |
|----------|---------------------|
| Fungi    | ≥ 1.0%              |
| Viruses  | ≥ 5.0%              |
| Bacteria | ≥ 5.0%              |
| Oomycota | ≥ 0.5%              |

MAL in-planta gate: `MAL_MIN_HOST_PCT = 1.0` (Viridiplantae reads ≥ 1%).
HAL pathogen gate:  `HAL_MIN_PATHOGEN_PCT = 1.0` (any PHI-base pathogen/virus ≥ 1%).

## Output schema (03_crypt.py → `{mode}_crypt.tsv`)

| Column | MAL | HAL |
|---|---|---|
| `Run` | SRR accession | SRR accession |
| `host` | top Viridiplantae leaf in STAT | library organism (ScientificName) |
| `host_pct` | Viridiplantae % in STAT | Viridiplantae % in STAT |
| `host_species_pct` | STAT leaf species % | — |
| `primary_pathogen` | library organism (ScientificName) | top PHI-base/ICTV pathogen in STAT |
| `primary_taxid` | resolved taxid (or blank) | resolved taxid |
| `primary_pct` | — (library dominates) | % of analysed reads |
| `secondary_pathogens` | other PHI-base/ICTV pathogens in STAT | remaining pathogens |
| `secondary_kingdoms` | kingdoms of secondaries (e.g. Viruses; Fungi) | same |
| `n_secondary` | count | count |
| `known_interaction` | PHI-base/ICTV confirms pathogen×host | same |
| `co_infection_flag` | single / multi_species / multi_kingdom | same |
| `fungi_pct` etc. | kingdom % from STAT | same |
| `analyzed` | total reads analysed by STAT | same |

Plus RunInfo passthrough: `BioProject`, `SRAStudy`, `Platform`, `ScientificName`.

`co_infection_flag` priority: `multi_kingdom` > `multi_species` > `single`.
`known_interaction` traces both pathogen and host to their PHI-base seed taxids
before checking the interaction map. For plant viruses, returns True for any
PHI-base plant host (mapped to all 180 host seeds in DB build).

## Entrez credentials
- API key: `NCBI_API_KEY` env var (set in `~/.bashrc`) → 10 req/s (9 used for Entrez)
- Without key: 2.5 req/s
- STAT endpoint rate: empirically ~10 req/s max; RATE=15 + 32 workers configured
- User-Agent: `crypt/01_sra (leon.lenzo@curtin.edu.au)` etc.

## Running

```bash
# Step 0: build PHI-base + ICTV reference DB
# MUST use system python3 (ete3 + sqlite3 incompatibility with miniconda)
python3 00_build.py               # auto-downloads PHI-base + ICTV VMR if absent
python3 00_build.py --fetch       # force re-download of both sources

# Steps 01-03: standard python (miniconda fine)

# Pre-check: count query hits before committing to a full fetch
python 01_sra.py --mode mal --count
python 01_sra.py --mode hal --count

# Fetch SRA run IDs (resumable; parallel POST-based RunInfo fetch)
python 01_sra.py --mode mal
python 01_sra.py --mode hal

# Fetch STAT + apply retention gate (resumable)
# Run MAL first, then HAL — they share stat_cache.json
python 02_stat.py --mode mal
python 02_stat.py --mode hal

# Identify cryptic co-infections + write TSV
python 03_crypt.py --mode mal
python 03_crypt.py --mode hal
```

## Output / cache layout

```
output/
├── 00_build/
│   ├── phi-base_current.csv     PHI-base CSV (fetched from GitHub)
│   ├── ictv_vmr.xlsx            ICTV VMR Excel (fetched from ictv.global/vmr/current)
│   ├── phibase_db.json          reference database (consumed by 01-03)
│   ├── build.log                full build log
│   └── _summary.txt             seed counts, taxid totals, name entries
│
├── 01_sra/
│   ├── {mode}_runs.json         RunInfo rows keyed by Run accession
│   ├── {mode}_uids.json         fetched UID set (for resumability)
│   ├── {mode}.log
│   └── {mode}_summary.txt       UIDs found, runs fetched
│
├── 02_stat/
│   ├── stat_cache.json          STAT responses — shared across MAL and HAL
│   ├── {mode}_confirmed.json    runs passing retention gate (with kingdom %)
│   ├── {mode}.log
│   └── {mode}_summary.txt       STAT coverage %, gate pass rate
│
└── 03_crypt/
    ├── {mode}_crypt.tsv         final co-infection candidate table
    ├── {mode}.log
    └── {mode}_summary.txt       flag breakdown, known/novel interactions, top secondaries
```

## Actual run results (2026-07)

| Mode | UIDs | Runs fetched | STAT fetches |
|------|------|-------------|--------------|
| MAL  | 45,312 | 47,171    | ~47k         |
| HAL  | 523,881 | 559,328  | ~559k (overnight) |

Both counts are upper bounds — runs matching taxa from multiple batches are counted
more than once. True deduplicated counts will be lower after RunInfo fetch.

esearch hit counts (upper bound, batched, may overcount cross-batch duplicates):
- MAL: ~45,721 (205 pathogen seeds, 5 batches of 50)
- HAL: ~539,909 (180 host seeds, 4 batches of 50)

## Performance notes

- **01_sra RunInfo fetch**: parallel POST requests (ThreadPoolExecutor, MAX_WORKERS=8).
  GET requests with 500 UIDs cause HTTP 414; POST avoids this. ~20 req/s achieved.
- **02_stat STAT fetch**: STAT endpoint (trace.ncbi.nlm.nih.gov) caps at ~10 req/s
  regardless of worker count. Configured MAX_WORKERS=32, RATE=15 to saturate naturally.
  Each STAT request takes ~1.8s; need ~18 workers to saturate 10 req/s.
  MAL (~47k): ~1–1.5 hrs. HAL (~559k): ~15–17 hrs.
- **Shared stat_cache**: do NOT run MAL and HAL 02_stat simultaneously — last writer
  wins on cache saves, causing redundant fetches and partial cache loss.
  Solution: chain HAL to start after MAL via tmux watcher:
  `until grep -q '02_stat MAL summary' output/02_stat/mal.log; do sleep 30; done`

## Notes

- `specific_hits()` leaf detection in 03_crypt.py is the core non-trivial algorithm
- MAL secondary detection excludes the primary/library organism by comparing
  PHI-base seed taxids (so strains of the same species are also excluded)
- `txid33090[Host]` returns 0 results in SRA — [Host] is a free-text field.
  In-planta confirmation is done via STAT Viridiplantae reads (host_pct).
- 00_build.py NCBITaxa SQLite schema: `species(taxid, parent, spname, common, rank, track)`
  and `synonym(taxid, spname)` — no `name_class` column (ete3-specific schema)
- ICTV VMR URL `https://ictv.global/vmr/current` redirects to the versioned file
  (e.g. VMR_MSL41.v1.20260721.xlsx) — the redirect target is stable per release
- Plant viruses in ICTV include insect-vectored species (geminiviruses, tospoviruses)
  via `Host source = "invertebrates, plants"` — included intentionally as they infect
  plants in the field even though the insect is the vector
