# crypt

Screens NCBI SRA plant RNA-seq runs for cryptic co-infections using NCBI STAT pre-computed k-mer taxonomy — no read re-alignment required.

> **Branches:** `master` — STAT-based pipeline (current results). `kraken` — in development; replaces STAT with Kraken2 screening of all ~593k step-01 runs on Setonix HPC; eukaryotic pathogens only (fungi, oomycetes).

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

### BioProject metadata and study design

BioProject XML and BioSample XML attributes were retrieved for all accessions in `runs.tsv` (`03a_fetch_xml.py`). PMIDs were resolved via a six-strategy short-circuit search (`03b_fetch_literature.py`): BioProject XML → PMC full-text → Europe PMC → ENA XML → Semantic Scholar → PubMed, stopping at the first hit. Study design was inferred from BioSample metadata and title/abstract keywords (`04_filter_kw.py`): treatment axis (`coinf_experiment`, `abiotic_stress`, `host_study`, `single`, `unclear`) and setting axis (`field`, `lab`, `unclear`). A per-BioProject LLM classification step (`05_llm_classify.py`, GPT-4o-mini) applied a finer-grained treatment vocabulary including a `surveillance` category not captured by keywords.

## Results

608,368 SRA runs screened across MAL and HAL. After retention gates and LibrarySource pre-filtering (GENOMIC/METAGENOMIC mislabelled runs excluded), 13,323 runs were confirmed across 1,754 BioProjects.

| Mode | Runs fetched | Non-RNA excluded | Gate pass | Gate pass rate |
|------|-------------|-----------------|-----------|----------------|
| MAL  | 48,418      | 1,419           | 6,191     | 13.4%          |
| HAL  | 559,950     | 7,343           | 7,118     | 1.3%           |

Of 11,117 BioSample-representative runs: 63.8% co-infected; 55.2% high-confidence (same_genus_secondary = False); 590 with `novel_host_range` status (pathogen confirmed in PHI-base but not previously recorded on the detected host). 401 BioProjects contained at least one co-infected run.

LLM study design classification (1,754 BioProjects): 747 single-pathogen, 403 host biology, 293 surveillance, 197 abiotic stress, 100 co-infection experiment. Field-classified BioProjects show a 2× higher co-infection detection rate than lab-classified BioProjects (21% vs 10% overall; 12% vs 2% high-confidence).

Key output files:

| File | Contents |
|------|----------|
| `output/02_filter_runs/data/runs.tsv` | One row per confirmed run; `mode`, `library_organism`, `stat_pathogens`, `co_infection_flag`, `interaction_status`, `same_genus_secondary`, `biosample_representative` |
| `output/04_filter_kw/data/biosample_kw.tsv` | One row per BioSample (11,117); BioProject metadata joined; `treatment`, `study_setting`, `named_host`, `primary_pmid` |
| `output/05_llm_classify/data/bioproject_llm.tsv` | One row per BioProject (1,754); `llm_treatment`, `llm_study_setting`, `llm_named_pathogen`, `llm_rationale` |

Filter to `biosample_representative == "True"` for sample-level statistics; additionally `same_genus_secondary == "False"` for highest-confidence co-detections.

## References

- **PHI-base**: [phi-base.org](https://phi-base.org)
- **ICTV VMR**: [ictv.global/vmr](https://ictv.global/vmr/current)
- **NCBI STAT**: Katz et al. (2021) *J Bioinform Comput Biol* — [PMC8450716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8450716/)

---

## Usage

### Dependencies

```bash
# Python — 00_build.py only; steps 01–05 are stdlib + requests
# Must use system python3 for 00_build.py (ete3/sqlite3 ABI conflict with miniconda)
pip install -r requirements.txt
```

```bash
# R — figure scripts (ggtree is Bioconductor; handled in install.R)
Rscript install.R
```

Optional environment variables:

```bash
export NCBI_API_KEY=...    # 10 req/s Entrez; 2.5 req/s without
export S2_API_KEY=...      # Semantic Scholar in 03b_fetch_literature.py; skipped without
export OPENAI_API_KEY=...  # Required for 05_llm_classify.py (~$2/full run, gpt-4o-mini)
```

### Quick start

```bash
# Step 0 — build reference DB (system python3 required)
python3 00_build.py

# Step 1 — fetch SRA runs + STAT profiles
# Run MAL before HAL (shared stat cache); use tmux for long runs
python 01_fetch_runs.py --mode mal
python 01_fetch_runs.py --mode hal

# Step 2 — apply retention gate, calibrate thresholds, classify co-infections
python 02_filter_runs.py                   # both modes, unified output
python 02_filter_runs.py --skip-validate   # after manual threshold review

# Step 3a — BioProject XML + BioSample XML attributes
python 03a_fetch_xml.py

# Step 3b — BioProject PMIDs + PMC methods sections
python 03b_fetch_literature.py

# Step 4 — keyword study design classification → biosample_kw.tsv
python 04_filter_kw.py

# Step 5 — LLM BioProject classification → bioproject_llm.tsv (requires OPENAI_API_KEY)
python 05_llm_classify.py
```

All steps are resumable via their respective caches.

---

*Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)*
