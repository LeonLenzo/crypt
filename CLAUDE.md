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

01_sra.py     Fetch SRA run IDs for MAL or HAL.
              MAL: batches 205 seed pathogen taxids, 50/query,
                   txid{t}[Organism:exp] AND RNA-Seq[Strategy]
              HAL: batches 180 seed host taxids, 50/query,
                   txid{t}[Organism:exp] AND RNA-Seq[Strategy]
              [Organism:exp] lets SRA expand to all descendant strains.
              RunInfo fetched via HTTP POST (avoids 414 on large UID batches).
              Parallel fetch: ThreadPoolExecutor(MAX_WORKERS=8) + shared
              rate lock. Saves RunInfo + UID cache for resumability.
              --count --mode mal shows per-kingdom hit breakdown.

02_stat.py    Fetch NCBI STAT for all runs (parallel, rate-limited).
              STAT endpoint: trace.ncbi.nlm.nih.gov (different from Entrez).
              Empirically caps at ~2–10 req/s (server load dependent);
              configured MAX_WORKERS=32, RATE=15.0 to saturate naturally.
              Do NOT run MAL and HAL simultaneously — they share
              stat_cache.json and the last writer wins on saves.
              Applies retention gate:
                MAL: keep if Viridiplantae reads ≥ MAL_MIN_HOST_PCT (1%)
                HAL: keep if any PHI-base pathogen/virus ≥ HAL_MIN_PATHOGEN_PCT (1%)

03_crypt.py   Identify cryptic co-infections in confirmed mixed libraries.
              Uses specific_hits() leaf-detection algorithm to find the most
              specific organism name in the STAT taxonomy tree under each
              kingdom node, then cross-references against PHI-base + ICTV.
              _genus() helper flags same-genus primary/secondary pairs
              (same_genus_secondary column) to mark k-mer bleed risk.
              Produces final TSV with co-infection classification.

04_meta.py    Fetch BioProject metadata for co-infection BioProjects.
              Reads {mode}_crypt.tsv from 03_crypt; by default filters to
              high-confidence runs (same_genus_secondary=False). Fetches
              via Entrez: BioProject title + description, linked PubMed
              articles (elink bioproject→pubmed), publication titles.
              Output: {mode}_bioproject_meta.tsv — one row per BioProject,
              sorted by run count. Resumable via {mode}_meta_cache.json.
              --all flag includes same-genus secondaries as well.
              Key curation step: distinguishes genuine cryptic co-infections
              from intentional co-infection study designs.
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
- Expanded to include all descendant strains via `_expand_taxids()`
- Viruses mapped to all 180 PHI-base plant host seeds in `pathogen_to_hosts`
  → `known_interaction(virus, any_plant_host)` returns True
- ICTV species names added to `name_to_taxid` as supplement to NCBI names
- Auto-downloaded from `https://ictv.global/vmr/current` (redirects to latest release)

**Combined DB stats** (post-rebuild 2026-07):
- 14,585 pathogen taxids (7,981 PHI-base + 6,604 ICTV viruses)
- 19,714 name entries — scientific names + synonyms + ICTV names, lowercase

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
`PhibaseDB` (in 03_crypt.py) derives `pathogen_taxids`, `host_taxids`,
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
- STAT shared cache at `output/02_stat/stat_cache.json` — keyed by Run accession

## STAT detection thresholds (hardcoded in 02_stat.py and 03_crypt.py)

| Kingdom  | Detection threshold |
|----------|---------------------|
| Fungi    | ≥ 1.0%              |
| Viruses  | ≥ 5.0%              |
| Bacteria | ≥ 5.0%              |
| Oomycota | ≥ 0.5%              |
| Nematoda | ≥ 1.0%              |

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
| `same_genus_secondary` | True if any secondary shares genus with primary | same |
| `fungi_pct` etc. | kingdom % from STAT | same |
| `analyzed` | total reads analysed by STAT | same |

Plus RunInfo passthrough: `BioProject`, `SRAStudy`, `Platform`, `ScientificName`.

`co_infection_flag` priority: `multi_kingdom` > `multi_species` > `single`.
`known_interaction` traces both pathogen and host to their PHI-base seed taxids
before checking the interaction map. For plant viruses, returns True for any
PHI-base plant host (mapped to all 180 host seeds in DB build).
`same_genus_secondary`: True when the first word of any secondary name matches the
first word of the primary — flags k-mer cross-mapping risk between sibling species.
Use `same_genus_secondary == False` to isolate high-confidence co-detections.

## k-mer bleed risk

STAT uses 31-mers. Sibling species within a genus share conserved k-mers (especially
rust fungi like Puccinia), so detecting a same-genus secondary may reflect k-mer
cross-mapping rather than true co-infection. Analysis of MAL data showed:
- 65% of secondary detections are same-genus — concentrated in Puccinia rusts
- Same-genus secondaries cluster near the detection threshold (<1%), consistent with bleed
- Cross-kingdom detections are biologically immune to this problem

**High-confidence co-infections**: filter to `same_genus_secondary == False` (diff-genus
secondary). These account for ~31% of co-infected runs but are the most credible signal.

## Entrez credentials
- API key: `NCBI_API_KEY` env var (set in `~/.bashrc`) → 10 req/s (9 used for Entrez)
- Without key: 2.5 req/s
- STAT endpoint rate: empirically ~2–10 req/s (variable); RATE=15 + 32 workers configured
- User-Agent: `crypt/01_sra (leon.lenzo@curtin.edu.au)` etc.

## Running

```bash
# Step 0: build PHI-base + ICTV reference DB
# MUST use system python3 (ete3 + sqlite3 incompatibility with miniconda)
python3 00_build.py               # auto-downloads PHI-base + ICTV VMR if absent
python3 00_build.py --fetch       # force re-download of both sources

# Steps 01-04: standard python (miniconda fine)

# Pre-check: count query hits before committing to a full fetch
python 01_sra.py --mode mal --count   # also shows per-kingdom breakdown
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

# Fetch BioProject metadata for high-confidence co-infection BioProjects
python 04_meta.py --mode mal
python 04_meta.py --mode hal
# --all flag includes same-genus secondaries (larger set, more noise)
```

### Extracting partial HAL results while 02_stat HAL is still running

While `02_stat.py --mode hal` runs in tmux, you can snapshot confirmed runs
from the partial stat_cache and run 03_crypt on them without interrupting the fetch:

```bash
# Read stat_cache (read-only), apply HAL gate, write hal_confirmed.json
python3 - << 'EOF'
import json, time
from pathlib import Path
# ... (see session notes — extract HAL confirmed from partial cache)
EOF
python 03_crypt.py --mode hal    # produces hal_crypt.tsv from partial set
# 02_stat.py will overwrite hal_confirmed.json when it finishes; re-run 03_crypt then
```

## Output / cache layout

```
output/
├── 00_build/
│   ├── phi-base_current.csv       PHI-base CSV (fetched from GitHub)
│   ├── ictv_vmr.xlsx              ICTV VMR Excel (fetched from ictv.global/vmr/current)
│   ├── phibase_db.json            reference database (consumed by 01-03)
│   ├── build.log
│   └── _summary.txt               seed counts, taxid totals, name entries
│
├── 01_sra/
│   ├── {mode}_runs.json           RunInfo rows keyed by Run accession
│   ├── {mode}_uids.json           fetched UID set (for resumability)
│   ├── {mode}.log
│   └── {mode}_summary.txt
│
├── 02_stat/
│   ├── stat_cache.json            STAT responses — shared across MAL and HAL
│   ├── {mode}_confirmed.json      runs passing retention gate (with kingdom %)
│   ├── {mode}.log
│   └── {mode}_summary.txt
│
├── 03_crypt/
│   ├── {mode}_crypt.tsv           final co-infection candidate table
│   ├── {mode}_summary.txt
│   └── {mode}.log
│
├── 04_meta/
│   ├── {mode}_bioproject_meta.tsv BioProject titles + PMIDs for HC co-infection projects
│   ├── {mode}_meta_cache.json     Entrez response cache (resumability)
│   ├── {mode}.log
│   └── {mode}_summary.txt
│
└── figure/
    ├── crypt_host_tree.nwk        NCBI cladogram of confirmed plant host species
    └── crypt_host_tree_meta.tsv   tip metadata: label, species, n_single, n_multi,
                                   n_confirmed, family
```

```
figure/
├── crypt_host_tree.py    Build nwk + meta from mal_crypt.tsv (system python3)
└── crypt_host_tree.R     Fan tree figure — green/orange bars per host species
                          (orange = co-infection, green = single pathogen)
                          Run from crypt/: Rscript figure/crypt_host_tree.R
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

- **01_sra RunInfo fetch**: parallel POST requests (ThreadPoolExecutor, MAX_WORKERS=8).
  GET requests with 500 UIDs cause HTTP 414; POST avoids this. ~20 req/s achieved.
- **02_stat STAT fetch**: STAT endpoint (trace.ncbi.nlm.nih.gov) is highly variable —
  observed ~2 req/s overnight to ~10 req/s during daytime.
  MAL (~48k): ~1.5 hrs. HAL (~559k): estimated 3–4 days depending on server load.
- **Shared stat_cache**: do NOT run MAL and HAL 02_stat simultaneously — last writer
  wins on cache saves. Chain HAL after MAL:
  `until grep -q '02_stat MAL summary' output/02_stat/mal.log; do sleep 30; done`

## Notes

- `specific_hits()` leaf detection in 03_crypt.py is the core non-trivial algorithm
- `_genus(name)` returns the first word (lowercase) — used to set `same_genus_secondary`
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
- `novel_interaction` currently cannot distinguish genuine novel co-infections from
  intentional co-infection study designs — BioProject metadata mining
  (→ mal_bioproject_meta.tsv) is the next step to resolve this
