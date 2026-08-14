# kraken — Orthogonal Kraken2 Species-level Validation

## Rationale

The STAT-based screening pipeline relies on pre-computed k-mer profiles built from the entire SRA at NCBI scale. While this enables whole-corpus screening, it introduces a known reliability problem: k-mers shared across closely related species are promoted to the lowest common ancestor (LCA) node in the STAT database, preventing species-level resolution for taxonomically dense clades. A kallisto pilot on 45 MAL runs dominated by *Puccinia striiformis* f. sp. *tritici* (PST, wheat yellow rust) confirmed this directly — STAT assigned zero eukaryotic percentage to all 15 PST-dominated runs, while both kallisto and Kraken2 detected PST at 65–68% of classified reads. STAT resolves PST reads to the Dikarya node (kingdom-level), not species.

This raises a general concern: STAT detections are reliable where pathogen k-mer sets are sufficiently unique at SRA scale, but blind spots exist for organisms that share abundant k-mers with close relatives in the database. Kraken2 with a purpose-built, curated reference addresses this by controlling the database composition explicitly.

## Database design

The Kraken2 database contains only PHI-base eukaryotic pathogen CDS sequences (fungi and oomycetes). Host sequences are deliberately excluded. This design choice reflects a key difference from the STAT approach: because the DB contains only pathogen sequences, `pct_classified` in Kraken2 output directly represents pathogen burden without requiring normalisation against host signal. Including host CDS would require the same Viridiplantae gate logic used in STAT MAL, reintroducing the ambiguity that the orthogonal approach is meant to resolve.

### Masking

Before database construction, each pathogen's CDS sequences are masked against k-mers shared with any other pathogen in the reference set using BBDuk (k=35, mincount=2). This species-diagnostic masking ensures that every k-mer retained in the database is unique among the included pathogens. Without masking, organisms sharing taxonomic order (e.g., *Melampsora* appearing in PST-dominated runs due to shared Pucciniales k-mers) produce false-positive detections. The masking step is pathogen-vs-pathogen only — no host sequences are used in masking.

### Host removal

Because the database contains only pathogen sequences, host reads do not match any database entry and are reported as `unclassified`. `pct_classified` in the Kraken2 output therefore directly represents pathogen burden in the library — it is equivalent to running host-removal followed by pathogen classification, without the computational overhead. This is the key architectural difference from approaches that build combined host+pathogen databases: those require an explicit host-removal step and normalisation against classified host counts. The pathogen-only approach is simpler, faster, and makes the pathogen fraction directly interpretable.

### Read subsampling

Each run is subsampled to 500,000 reads from the ENA FTP stream. This cap is sufficient for robust species-level detection at all thresholds used: at the lowest threshold (0.5% for Fungi), the expected signal is ≥ 2,500 reads — well within reliable detection range for Kraken2 with confidence=0.15 and min-hit-groups=3. Reads arrive in the original ENA submission order (FASTQ deposit is not sorted by quality or coverage), so the subsample is representative of the full library composition. The cap limits ENA concurrent bandwidth load and reduces per-run runtime without sacrificing sensitivity at the thresholds of interest.

### Coverage and gaps

| Component | Source | Annotated seeds | Database |
|-----------|--------|-----------------|----------|
| Fungi | PHI-base plant entries | 101 species with CDS annotation | fungi + oomycetes only |
| Oomycota | PHI-base plant entries | included above | |

Seeds lacking NCBI annotation (no gene model, therefore no `cds_from_genomic.fna`) are excluded — without CDS coordinates there is no transcript-compatible sequence to add. *Nicotiana benthamiana*, the most widely used model host, has no annotated assembly in NCBI Datasets as of 2026-08 despite published chromosome-level assemblies (Bally et al. 2022); this is a known gap not addressable without manual assembly integration.

CDS-based sequences are used throughout rather than whole-genome sequences, for two reasons: (1) RNA-seq reads derive from spliced transcripts and align poorly to genomic sequence across introns; (2) CDS sequences are shorter and more specific, reducing database size and false-positive rates.

## Control validation set

Before running Kraken2 classification across the full 10,995-run corpus, a stratified control validation set (2,473 runs) was constructed from `runs.tsv` to characterise the relationship between STAT detections and Kraken2 detections across the parameter space. Runs were stratified into strata A–I defined by STAT eukaryotic signal level and biosample representation, with target sizes ensuring adequate coverage of both STAT-positive (co-infected) and STAT-negative (single) biosample categories. The control set is designed to answer: where STAT and Kraken2 agree, how consistent is the quantitative signal? Where they disagree, which is more reliable?

The 2,473-run control set is built (`kraken/control/output/data/run_ids.txt`) and ready to submit on Setonix:

```bash
git pull && sbatch kraken/slurm/kraken_classify_control.slurm
```

Production classification of the full corpus follows after control results validate the DB and parameter choices.

## Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Reads per run | 500,000 | Sufficient for species-level detection; limits ENA FTP bandwidth per run |
| Kraken2 confidence | 0.15 | Standard setting; reduces spurious species-level calls at root/genus level |
| Min hit groups | 3 | Requires minimers from ≥ 3 independent positions; reduces single-fragment hits |
| Workers | 8 | Respects ENA concurrent connection guidelines |
| Threads per worker | 4 | Workers × threads matched to available cores on Setonix |

Classification is performed on Setonix HPC (Pawsey Supercomputing Centre) via SLURM, streaming reads directly from ENA FTP without local storage. Results are written to an append-only cache (`classify/data/kraken_cache.jsonl`) enabling resumable runs.

## Pilot results

A 45-run pilot (2026-08-05) followed by a 32-run high-confidence validation set (2026-08-07) confirmed:

- STAT correctly detects the majority of HAL-mode eukaryotic co-infections where the secondary organism has a sufficiently distinct k-mer profile
- The PST blind spot is confined to rust fungi (Pucciniales); non-rust fungal detections in the pilot set showed broad STAT–Kraken2 agreement
- The pathogen-only DB (no host sequences) produces cleaner `pct_classified` values than the earlier masked+host DB pilot, where host background introduced noise at low pathogen burden thresholds

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
