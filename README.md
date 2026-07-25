# crypt

Screens NCBI SRA RNA-seq runs for cryptic plant co-infections using NCBI STAT pre-computed k-mer taxonomy — no read alignment required.

## Background

Field plant samples submitted to SRA for a single pathogen study frequently contain co-infecting organisms that go unreported. This pipeline detects those secondary pathogens by cross-referencing STAT k-mer taxonomy profiles against PHI-base (plant–pathogen interactions) and ICTV (plant viruses).

## Two screening modes

| Mode | Library organism | What we detect | Gate |
|------|-----------------|----------------|------|
| **MAL** (microbe-as-library) | PHI-base plant pathogen | Secondary pathogens co-infecting the host plant | Viridiplantae ≥ 1% of STAT reads confirms the run is in planta |
| **HAL** (host-as-library) | PHI-base plant host | Primary + secondary pathogens in the host transcriptome | Any PHI-base pathogen ≥ 1% of STAT reads |

## Dependencies

```
python >= 3.11          (miniconda fine for steps 01–03)
system python3          (required for 00_build.py — ete3 + sqlite3 conflict)
ete3                    (taxonomy resolution in 00_build.py)
openpyxl                (ICTV VMR parsing in 00_build.py)

R packages (figure only):
  ggtree, ape, ggplot2, dplyr, tibble, cowplot
```

NCBI API key recommended (set `NCBI_API_KEY` in `~/.bashrc`); without it rate limits to 2.5 req/s.

## Quick start

```bash
# 1. Build reference DB (system python3 required)
python3 00_build.py

# 2. Fetch SRA run IDs
python 01_sra.py --mode mal
python 01_sra.py --mode hal

# 3. Fetch STAT taxonomy + apply gate (run MAL first — shared cache)
python 02_stat.py --mode mal
python 02_stat.py --mode hal

# 4. Classify co-infections
python 03_crypt.py --mode mal
python 03_crypt.py --mode hal

# 5. Fetch BioProject metadata for high-confidence co-infection BioProjects
python 04_meta.py --mode mal
python 04_meta.py --mode hal
```

Steps 01–03 are resumable. Long-running STAT fetches should be run in tmux.

## Output

`output/03_crypt/{mode}_crypt.tsv` — one row per confirmed run, columns include:

| Column | Description |
|--------|-------------|
| `host` | Plant host species (STAT-detected for MAL; library organism for HAL) |
| `primary_pathogen` | Library organism (MAL) or top detected pathogen (HAL) |
| `secondary_pathogens` | Additional PHI-base/ICTV pathogens detected in STAT |
| `co_infection_flag` | `single` / `multi_species` / `multi_kingdom` |
| `same_genus_secondary` | `True` if secondary shares genus with primary (k-mer bleed risk) |
| `known_interaction` | PHI-base confirms this pathogen × host pair |

Filter to `co_infection_flag != "single"` for co-infected runs.
Filter additionally to `same_genus_secondary == "False"` for high-confidence detections.

## Results (2026-07)

**MAL** (complete):
- 48,418 runs screened → 6,852 confirmed (14.2%)
- 1,148 co-infected (16.8% of confirmed)
- **360 high-confidence** co-infections (diff-genus secondary) across **54 BioProjects**
- Top secondaries: *Pepino mosaic virus*, *Zymoseptoria tritici*, *Alternaria alternata*

**HAL** (partial — 63k / 559k runs screened):
- 2.2% gate pass rate
- 98 co-infected (12.4% of confirmed); 32 high-confidence across 21 BioProjects
- Top secondaries: *Peanut stripe virus*, *Bipolaris zeicola*, *Pectobacterium brasiliense*

## Key caveat — k-mer bleed

STAT uses 31-mers. Closely related species (same genus) share conserved k-mers, so a same-genus secondary detection may reflect cross-mapping rather than true co-infection. ~65–69% of secondary detections are same-genus in both modes. Cross-kingdom detections are the most reliable signal.

## Figure

```bash
# Build tree data (system python3)
python3 figure/crypt_host_tree.py

# Render fan tree (from crypt/)
Rscript figure/crypt_host_tree.R
```

Produces a fan phylogeny of confirmed plant host species with bars showing gate-passed run counts (orange = co-infection, green = single pathogen detected).

## Reference databases

- **PHI-base**: [phi-base.org](https://phi-base.org) — plant–pathogen interaction database
- **ICTV VMR**: [ictv.global](https://ictv.global/vmr/current) — plant virus master species list
- **NCBI STAT**: pre-computed k-mer taxonomy for all SRA runs

## Author

Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)
