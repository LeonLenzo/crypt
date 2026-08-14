# kraken — Orthogonal Kraken2 Species-level Validation

## Rationale

The STAT-based screening pipeline relies on pre-computed k-mer profiles built from the entire SRA at NCBI scale. While this enables whole-corpus screening, it introduces a known reliability problem: k-mers shared across closely related species are promoted to the lowest common ancestor (LCA) node in the STAT database, preventing species-level resolution for taxonomically dense clades. A kallisto pilot on 45 MAL runs dominated by *Puccinia striiformis* f. sp. *tritici* (PST, wheat yellow rust) confirmed this directly — STAT assigned zero eukaryotic percentage to all 15 PST-dominated runs, while both kallisto and Kraken2 detected PST at 65–68% of classified reads. STAT resolves PST reads to the Dikarya node (kingdom-level), not species.

This raises a general concern: STAT detections are reliable where pathogen k-mer sets are sufficiently unique at SRA scale, but blind spots exist for organisms that share abundant k-mers with close relatives in the database. Kraken2 with a purpose-built, curated reference addresses this by controlling the database composition explicitly.

## Database design

The Kraken2 database contains only PHI-base eukaryotic pathogen CDS sequences (fungi and oomycetes). Host sequences are deliberately excluded — because the DB contains only pathogen sequences, `pct_classified` in Kraken2 output directly represents pathogen burden without requiring normalisation against host signal. Host reads simply go unclassified. Including host CDS would require the same Viridiplantae gate logic used in STAT MAL, reintroducing the ambiguity the orthogonal approach is meant to resolve.

### Genus fill-in (no masking)

Rather than masking to create species-diagnostic k-mers, the DB uses **genus-level completeness** to achieve specificity. For each PHI-base genus detected in the pilot, the best-annotated assembly for every non-PHI-base species within that genus is also included. This genus fill-in means that k-mers shared between species within a genus resolve — via Kraken2's LCA algorithm — to the genus node rather than being spuriously assigned to the wrong species. Species-level detections therefore represent reads for which no within-genus congener claims the k-mer, providing a natural specificity filter without manual masking.

### Read acquisition and subsampling

Reads are downloaded from ENA FTP in advance (`download_control.slurm`) and stored locally on Setonix scratch, then classified from local files (`classify_control.slurm`). Each run is subsampled to 500,000 reads. This cap is sufficient for robust species-level detection at all thresholds used: at the lowest threshold (0.5% for Fungi), the expected signal is ≥ 2,500 reads — well within reliable detection range for Kraken2 with confidence=0.15 and min-hit-groups=3.

### Coverage and gaps

CDS-based sequences are used throughout rather than whole-genome sequences: RNA-seq reads derive from spliced transcripts and align poorly to genomic sequence across introns, and CDS sequences are shorter and more specific. Seeds lacking NCBI annotation (no gene model) are excluded — without CDS coordinates there is no transcript-compatible sequence to add.

## Control validation

Validation uses two complementary components:

### Stratified SRA controls (`control/sample.py`)

A stratified set of 2,473 real SRA runs was drawn from the full corpus to characterise STAT–Kraken2 agreement across the parameter space. Runs are stratified into strata A–I by STAT eukaryotic signal level and study design (LLM classification):

| Stratum | Description | Target |
|---------|-------------|--------|
| A | Non-plant true negative (external; non-plant host) | 200 |
| B | Plant host-only — HAL gate-fail, 0% fungi/oomycete | 500 |
| B2 | Pathogen-only — MAL gate-fail, 0% host reads | 200 |
| C | Single pathogen, lab, low burden (STAT 0.5–2%) | 200 |
| D | Single pathogen, lab, medium burden (2–10%) | 200 |
| E | Single pathogen, lab, high burden (>10%) | 200 |
| F | Single pathogen, field | 400 |
| G | Intentional co-infection experiment | all |
| H | HC field co-infected, diff-genus | 400 |
| I | Same-genus secondary (LCA collapse test) | all |

These runs have STAT-inferred ground truth (not known composition); the set is designed to answer: where STAT and Kraken2 agree, how consistent is the quantitative signal? Where they disagree, which is more reliable?

Reads are pre-downloaded to Setonix scratch via `download_control.slurm`, then classified via `classify_control.slurm`. Run from Setonix:

```bash
# Step 1 — download reads
sbatch kraken/slurm/download_control.slurm
# Step 2 — classify (after download completes)
sbatch kraken/slurm/classify_control.slurm
```

### In silico controls (`control/insilico.py`, stratum J)

To obtain controls with **known composition** (not STAT-inferred ground truth), synthetic mixed-infection libraries are constructed by mixing reads from pure-culture SRA runs into a plant host background at specified ratios (0.1%, 0.5%, 1%, 5%, 10%, 25%). For obligate biotrophs without pure-culture RNA-seq (e.g., PST), reads are simulated from the reference genome using ART Illumina prior to mixing. These in silico controls directly calibrate KINGDOM_THRESHOLDS and validate same-genus specificity. Kristina Cihatova is collaborating on the panel design and ART-simulated read generation for the rust fungi panel.

Production classification of the full corpus follows after control results validate the DB and parameter choices.

## Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Reads per run | 500,000 | Sufficient for species-level detection; limits ENA FTP bandwidth per run |
| Kraken2 confidence | 0.15 | Standard setting; reduces spurious species-level calls at root/genus level |
| Min hit groups | 3 | Requires minimers from ≥ 3 independent positions; reduces single-fragment hits |
| Workers | 8 | Respects ENA concurrent connection guidelines |
| Threads per worker | 4 | Workers × threads matched to available cores on Setonix |

Classification is performed on Setonix HPC (Pawsey Supercomputing Centre) via SLURM. Reads are downloaded from ENA FTP to Setonix scratch in a separate step, then classified from local files. Results are written to an append-only cache (`classify/data/kraken_cache.jsonl`) enabling resumable runs.

## Pilot results

A 45-run pilot (2026-08-05) followed by a 32-run high-confidence validation set (2026-08-07) confirmed:

- STAT correctly detects the majority of HAL-mode eukaryotic co-infections where the secondary organism has a sufficiently distinct k-mer profile
- The PST blind spot is confined to rust fungi (Pucciniales); non-rust fungal detections in the pilot set showed broad STAT–Kraken2 agreement
- The pathogen-only DB (no host sequences) produces cleaner `pct_classified` values than earlier host-inclusive DB pilots, where host background introduced noise at low pathogen burden thresholds

## Limitations

**ENA FTP bandwidth.** Streaming 500,000 reads per run from ENA for ~10,000 runs is network-intensive and subject to ENA connection rate limits. The classify step requires Setonix compute access and takes on the order of days for the full corpus.

**CDS-only database misses splicing junctions.** Kraken2 k-mers span 35 bp; reads that bridge exon–exon junctions (absent from CDS sequences) will not match database k-mers, slightly reducing sensitivity. In practice this effect is small for highly expressed genes, where most reads fall within exons.

**Database completeness.** The 101 annotated seed species represent a subset of the 205 PHI-base seed taxids; seeds lacking annotation are not in the database. Co-infections involving unannotated or poorly assembled pathogens are invisible to Kraken2 but detectable by STAT (which uses the full SRA k-mer set including environmental and low-coverage organisms).

**Rust blind spot not fully characterised.** The STAT PST blind spot was confirmed by the pilot, but the extent to which other rust clades (e.g., *Melampsora*, *Hemileia*) share this problem has not been systematically tested. The control validation set includes rust-dominated strata specifically to characterise this.

## Key output files

| File | Contents |
|------|----------|
| `kraken/output/classify/data/kraken_cache.jsonl` | Append-only classification cache; one JSON object per run |
| `kraken/output/classify/data/reports/` | Raw Kraken2 report.txt files (pilot runs; full taxonomy tree) |
| `kraken/control/output/data/run_ids.txt` | 2,473-run stratified control validation set |
| `kraken/control/output/data/control_runs.tsv` | Control set with stratum assignments and STAT annotations |
