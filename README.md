# crypt — Mining cryptic co-infections from publicly deposited plant RNA-seq

## The problem

When a wheat field shows yellow rust symptoms, agronomists know *Puccinia striiformis* is involved — but they also know that plants in the field are rarely fighting a single enemy. Co-infections are common in agricultural settings, influencing disease severity, epidemic dynamics, and host immune responses. Despite this, most published plant disease studies are designed around a single target pathogen, and that is what gets reported.

Over the past decade, plant pathologists have deposited hundreds of thousands of RNA-seq experiments in NCBI SRA. These runs were sequenced to study pathogen gene expression, monitor field populations, or characterise host resistance — but they captured everything that was in the sample. Reads from an unreported co-infecting *Zymoseptoria tritici*, *Pepino mosaic virus*, or *Botrytis cinerea* are sitting in publicly available data, unannounced.

**crypt** mines that signal without re-aligning a single read.

## How it works

[NCBI STAT](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-examples/stat/) (Sequence Taxonomy Analysis Tool) provides a pre-computed k-mer taxonomy fingerprint for every run deposited in SRA. Rather than aligning reads to reference genomes — which requires knowing what you are looking for — STAT decomposes sequencing reads into 32-nucleotide words and maps them against a comprehensive reference database spanning the full tree of life. The result is a table of organisms with their proportional read contributions, available for any of the ~30 million SRA runs, without downloading raw data.

We query SRA for runs associated with known plant pathogens or plant hosts, retrieve their STAT profiles, and cross-reference against two curated databases:

- **PHI-base** — experimentally confirmed plant–pathogen interactions (fungi, bacteria, oomycetes, nematodes)
- **ICTV VMR** — the full ICTV plant virus master species list (~2,600 plant virus species)

Any organism detected in STAT above kingdom-specific abundance thresholds, confirmed in PHI-base or ICTV but not reported as the study target, is a candidate cryptic co-infection.

## Two screening modes

The pipeline runs two complementary screens that query SRA from opposite directions:

| Mode | Full name | Library organism | Retention gate | What we detect |
|------|-----------|-----------------|----------------|----------------|
| **MAL** | Microbe-as-library | PHI-base plant pathogen | Viridiplantae ≥ 1% of STAT reads (confirms run is in planta) | Secondary pathogens co-infecting the same host plant |
| **HAL** | Host-as-library | PHI-base plant host species | Any PHI-base pathogen ≥ 1% of STAT reads | Primary and secondary pathogens present in host transcriptomes |

MAL targets studies that sequenced a pathogen directly — disease surveys, population genomics, pathotype characterisation. The Viridiplantae gate filters out pure-culture and axenic samples, keeping only runs where plant tissue is clearly present.

HAL targets studies that sequenced the host plant — transcriptomics of infected tissue, resistance gene expression, field transcriptomics. These runs already contain plant reads, so the gate instead confirms at least one PHI-base pathogen is detectable.

Together the modes cover the two dominant designs in plant disease RNA-seq and are largely non-overlapping.

## Confidence in detections

Not all secondary detections carry equal weight. The pipeline flags two tiers:

**High-confidence** (`same_genus_secondary = False`): the secondary pathogen belongs to a different genus from the primary. These are the most interpretable — a wheat rust sample (*Puccinia* spp.) that also shows signal for *Zymoseptoria tritici* or *Pepino mosaic virus* is very unlikely to be a database artefact.

**Lower-confidence** (`same_genus_secondary = True`): primary and secondary share a genus (e.g. *Puccinia graminis* + *Puccinia striiformis*). STAT uses an LCA-based k-mer merging strategy at database build time: k-mers shared between sibling species are promoted to the genus node and not assigned to either species individually. Species-level detections therefore represent genuinely diagnostic k-mers, but closely related species have smaller unique-k-mer sets and sit closer to the detection threshold. These are real detections but warrant more careful biological interpretation.

Co-infections are further classified by `interaction_status`:
- `known` — PHI-base confirms this pathogen × host combination
- `novel_host_range` — pathogen is in PHI-base, but not recorded on this host species (**primary signal of interest**)
- `novel_combination` — pathogen is in PHI-base with no recorded hosts at all
- `unresolved` — taxid could not be resolved; novelty cannot be assessed

## Dependencies

```
python >= 3.11          standard library only for steps 01–03 (miniconda fine)
system python3          required for 00_build.py (ete3 + sqlite3 ABI conflict)
ete3                    NCBI taxonomy resolution (00_build.py only)
openpyxl                ICTV VMR Excel parsing (00_build.py only)

R packages (figure scripts only):
  ggraph, igraph, ggplot2, dplyr, tibble, ggtree, ape, cowplot
```

NCBI API key strongly recommended — set `NCBI_API_KEY` in `~/.bashrc`. Without it, Entrez rate-limits to 2.5 req/s; with it, 10 req/s. The STAT fetch (step 01) hits a separate endpoint and is not Entrez-rate-limited.

## Quick start

```bash
# Step 0: build PHI-base + ICTV reference DB
# Must use system python3 (ete3/miniconda sqlite3 conflict)
python3 00_build.py               # auto-downloads PHI-base CSV and ICTV VMR

# Steps 1–3: standard python (miniconda fine)

# Step 1: fetch SRA run IDs + STAT taxonomy
# Run MAL first, then HAL — they share a stat cache
# Use tmux — these run for hours; _Tee writes log, do NOT add | tee
python 01_fetch.py --mode mal
python 01_fetch.py --mode hal

# Step 2: apply retention gate, calibrate thresholds, classify co-infections
python 02_filter.py               # both modes, unified output
python 02_filter.py --mode mal    # single mode

# Step 3: fetch BioProject metadata and identify associated publications
python 03_find.py
python 03_find.py --hc            # restrict to same_genus_secondary=False
```

Steps 01 and 02 are fully resumable — interrupt and restart freely.

## Output

`output/02_filter/data/crypt.tsv` — one row per confirmed run:

| Column | Description |
|--------|-------------|
| `mode` | `mal` or `hal` |
| `host` | Plant host (STAT-detected leaf species for MAL; library organism for HAL) |
| `primary_pathogen` | Library organism (MAL) or top STAT-detected pathogen (HAL) |
| `secondary_pathogens` | Additional PHI-base/ICTV organisms detected; semicolon-separated |
| `co_infection_flag` | `single` / `multi_species` / `multi_kingdom` |
| `same_genus_secondary` | `True` if any secondary shares genus with primary |
| `interaction_status` | `known` / `novel_host_range` / `novel_combination` / `unresolved` |
| `biosample_representative` | `True` for one run per biological sample (use for sample-level stats) |

Filter to `co_infection_flag != "single"` for co-infected runs.  
Filter additionally to `same_genus_secondary == "False"` for high-confidence detections.  
Filter to `biosample_representative == "True"` for biological-sample-level statistics.

`output/03_find/data/bioproject_meta.tsv` — one row per co-infection BioProject:

| Column | Description |
|--------|-------------|
| `study_design` | `coinf_experiment` / `field_survey` / `unclear` (keyword inference) |
| `primary_pmid` | Earliest linked publication (depositing paper heuristic) |
| `modes` | `mal` / `hal` / `mal+hal` |

`coinf_experiment` BioProjects require manual methods-section review to confirm whether the co-infection was intentional (experimental design) or incidental (cryptic). Intentional designs should be excluded from novel interaction counts.

## Preliminary results (MAL complete, HAL in progress — 2026-07)

**MAL** (48,418 runs screened):
- 6,852 confirmed in planta (14.4%)
- 902 co-infected runs across 118 BioProjects
- **High-confidence** (diff-genus secondary): subset of above, across field surveillance and pathogenomics projects
- Top secondaries: *Puccinia graminis*, *Zymoseptoria tritici*, *Pepino mosaic virus*, *Alternaria alternata*, *Botrytis cinerea*
- Interaction status: ~39% `known`, ~18% `novel_host_range`, ~27% `unresolved`

**HAL** (559,328 runs queried; STAT fetch in progress):
- Partial results show ~2.2% gate pass rate (expected — host transcriptomes contain less pathogen signal than pathogen-focused libraries)
- Biology is complementary to MAL: legume viruses, apple/grapevine pathogens, corn and beet pathogens — largely non-overlapping with MAL

**BioProject metadata**: 64 / 118 MAL BioProjects (54%) linked to publications. Remaining 46% are likely unpublished or pre-publication datasets.

## Reference databases

- **PHI-base**: [phi-base.org](https://phi-base.org) — Urban et al., plant–pathogen interaction database
- **ICTV VMR**: [ictv.global/vmr](https://ictv.global/vmr/current) — ICTV virus master species list
- **NCBI STAT**: Katz et al. (2021), *J Bioinform Comput Biol* — [PMC8450716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8450716/)

## Author

Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)
