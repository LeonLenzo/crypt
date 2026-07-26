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

# Fetch BioProject metadata and link publications (7 sources)
python 03_find.py
python 03_find.py --hc            # restrict to same_genus_secondary=False
```

Steps 01 and 02 are resumable. Long-running STAT fetches should be run in tmux.

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

`output/03_find/data/bioproject_meta.tsv` — one row per co-infection BioProject, including `study_design` (`coinf_experiment` / `field_survey` / `unclear`) inferred from BioProject title, description, and linked publication abstracts. `coinf_experiment` projects require manual methods review before inclusion in novel interaction counts.

## Preliminary results (MAL complete, HAL in progress — 2026-07)

**MAL** (48,418 runs screened):
- 6,852 confirmed in planta (14.4%); 902 co-infected across 118 BioProjects
- Top secondaries: *Puccinia graminis*, *Zymoseptoria tritici*, *Pepino mosaic virus*, *Alternaria alternata*, *Botrytis cinerea*
- Interaction status: ~39% `known`, ~18% `novel_host_range`, ~27% `unresolved`

**HAL** (559,328 runs; STAT fetch in progress):
- Partial screen (~2.2% gate pass rate) reveals complementary biology to MAL: legume viruses, apple and grapevine pathogens, corn and beet pathogens

**BioProject metadata**: 64 / 118 MAL BioProjects linked to publications via 7-source PMID search; remaining ~46% are likely pre-publication or data-only submissions.

## Interpreting `same_genus_secondary`

STAT uses an LCA-based k-mer merging strategy: k-mers shared between sibling species are promoted to the genus node at database build time and are not assigned to either species. Species-level detections therefore represent genuinely species-diagnostic signal — inter-species k-mer bleed within a genus is largely prevented by design. Same-genus secondaries are lower confidence not because of cross-mapping, but because closely related species retain fewer unique diagnostic k-mers after LCA merging, pushing their signal closer to the detection threshold. Cross-kingdom co-detections are the highest-confidence signal and biologically unambiguous.

## Reference databases

- **PHI-base**: [phi-base.org](https://phi-base.org)
- **ICTV VMR**: [ictv.global/vmr](https://ictv.global/vmr/current)
- **NCBI STAT**: Katz et al. (2021) *J Bioinform Comput Biol* — [PMC8450716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8450716/)

## Author

Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)
