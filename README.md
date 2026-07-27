# crypt

Screens NCBI SRA plant RNA-seq runs for cryptic co-infections using NCBI STAT pre-computed k-mer taxonomy — no read re-alignment required.

## Background

Plant disease studies deposited in SRA are designed around a single target pathogen. Field-collected samples, however, routinely harbour additional co-infecting organisms that go undetected and unreported under single-target study designs. We hypothesise that a substantial fraction of publicly available plant RNA-seq data contains secondary pathogen signal sufficient for detection via k-mer taxonomy, representing an untapped resource for co-infection epidemiology.

NCBI STAT provides pre-computed 32-mer taxonomy profiles for all SRA runs. By cross-referencing STAT outputs against PHI-base (plant–pathogen interactions) and the ICTV plant virus master species list, secondary pathogens can be identified in runs where they were not the study target — without downloading or re-aligning raw reads.

## Methods

### Screening modes

Two complementary strategies screen different segments of the public RNA-seq archive:

| Mode | Library organism | Retention gate | Signal |
|------|-----------------|----------------|--------|
| Microbe-as-library (MAL) | PHI-base plant pathogen | Viridiplantae ≥ 1% of STAT reads | Secondary pathogens co-infecting the host |
| Host-as-library (HAL) | PHI-base plant host | Any PHI-base pathogen/virus ≥ 1% of STAT reads | Pathogens present in host transcriptomes |

MAL targets pathogen-focused sequencing (population genomics, disease surveys); HAL targets host-focused sequencing (resistance transcriptomics, field transcriptomics). The two modes are complementary and largely non-overlapping.

### Reference databases

PHI-base provided plant–pathogen interaction records for 205 fungal, bacterial, oomycete, and nematode seed species, expanded via NCBI taxonomy to all descendant strains and formae speciales. The ICTV Virus Metadata Resource supplied ~2,630 plant virus species resolved to NCBI taxids. Both sources were merged into a local reference database (`00_build.py`) using `ete3 NCBITaxa`. Taxid expansion used `intermediate_nodes=True` to capture species nodes that have named strains as children in the NCBI hierarchy — required for rust f. sp. taxa and many ICTV virus entries.

### SRA query and STAT fetch

SRA accessions were retrieved via Entrez using `txid{taxid}[Organism:exp]`, which expands to all descendant strains without requiring free-text `[Host]` field matching. MAL queried 205 pathogen seed taxids; HAL queried 180 host seed taxids. STAT k-mer profiles were fetched in parallel from `trace.ncbi.nlm.nih.gov` and written to a shared append-only cache (`stat_cache.jsonl`), enabling resumable fetches and cache reuse across modes.

### Co-infection classification

For each retained run, `specific_hits()` identifies leaf-level species detections by finding k-mer counts not nested under any higher-specificity count — returning species-diagnostic signal rather than genus-level aggregates. Detected species are cross-referenced against the reference database and classified as:

- `known` — interaction recorded in PHI-base
- `novel_host_range` — pathogen known to PHI-base but not on this host species
- `novel_combination` — both organisms known but interaction not recorded
- `unresolved` — taxid lookup failed; novelty cannot be assessed

A `same_genus_secondary` flag marks cases where a secondary shares a genus with the primary. STAT's LCA k-mer design promotes shared k-mers to the genus node at database build time, preventing inter-species bleed within a genus by design. Same-genus secondaries are lower confidence not because of cross-mapping, but because closely related species retain fewer unique diagnostic k-mers after merging and cluster near the detection threshold. Cross-kingdom co-detections are biologically unambiguous. Runs from the same BioSample are deduplicated; `biosample_representative` marks the single highest-coverage run per sample.

### BioProject metadata

BioProject titles, submission dates, and linked publications were retrieved for all BioProjects in `crypt.tsv` — including single-infection projects, which serve as the denominator for field co-infection prevalence estimates. PMIDs were resolved via a six-strategy short-circuit search (BioProject XML → PMC full-text → Europe PMC → ENA XML → Semantic Scholar → PubMed), stopping at the first hit. Study design was inferred from title and abstract keywords: `coinf_experiment` (intentional mixed-infection — exclude from novel interaction counts), `field_survey`, or `unclear`.

## Results

607,746 SRA runs screened across MAL and HAL. After retention gates, 14,335 runs (11,853 BioSample-representative) across 1,797 BioProjects were confirmed.

| Mode | Runs screened | Confirmed | Gate pass |
|------|--------------|-----------|-----------|
| MAL  | 48,418       | 6,852     | 14.2%     |
| HAL  | 559,328      | ~7,483    | ~1.3%     |

Of 11,853 representative BioSamples: 10,078 (85%) single pathogen; 1,196 (10%) known co-infection; 579 (5%) novel host range — pathogen confirmed in PHI-base but not previously recorded on the detected host. 436 BioProjects contained at least one co-infected run.

Key output files:

| File | Contents |
|------|----------|
| `output/02_filter/data/crypt.tsv` | One row per confirmed run; columns: `mode`, `host`, `primary_pathogen`, `secondary_pathogens`, `co_infection_flag`, `interaction_status`, `same_genus_secondary`, `biosample_representative` |
| `output/03_find/data/bioproject_meta.tsv` | One row per BioProject; columns: `n_coinf`, `n_single`, `coinf_rate`, `study_design`, `primary_pmid` |

Filter to `biosample_representative == "True"` for sample-level statistics; additionally `same_genus_secondary == "False"` for highest-confidence co-detections.

## References

- **PHI-base**: [phi-base.org](https://phi-base.org)
- **ICTV VMR**: [ictv.global/vmr](https://ictv.global/vmr/current)
- **NCBI STAT**: Katz et al. (2021) *J Bioinform Comput Biol* — [PMC8450716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8450716/)

---

## Usage

### Dependencies

```bash
# Python — 00_build.py only; steps 01–03 are stdlib
# Must use system python3 (ete3/sqlite3 ABI conflict with miniconda)
pip install -r requirements.txt
```

```bash
# R — figure scripts (ggtree is Bioconductor; handled in install.R)
Rscript install.R
```

Optional environment variables:

```bash
export NCBI_API_KEY=...   # 10 req/s Entrez; 2.5 req/s without
export S2_API_KEY=...     # Semantic Scholar in 03_find.py; skipped without
```

### Quick start

```bash
# Step 0 — build reference DB (system python3 required)
python3 00_build.py

# Step 1 — fetch SRA runs + STAT profiles
# Run MAL before HAL (shared stat cache); use tmux for long runs
python 01_fetch.py --mode mal
python 01_fetch.py --mode hal

# Step 2 — apply retention gate, calibrate thresholds, classify co-infections
python 02_filter.py                   # both modes, unified output
python 02_filter.py --skip-validate   # after manual threshold review

# Step 3 — BioProject metadata + publication links
python 03_find.py
python 03_find.py --hc                # restrict to same_genus_secondary=False
```

All steps are resumable via their respective caches.

---

*Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)*
