# stat — STAT k-mer Screening

## Rationale

NCBI STAT provides pre-computed 32-mer taxonomic profiles for every run deposited in the SRA. Because these profiles are generated at NCBI infrastructure scale and cached server-side, they can be queried via a REST API without downloading or re-aligning raw reads. For a corpus of ~600,000 plant-associated RNA-seq runs, this is the only practical approach to whole-corpus screening — a Kraken2 equivalent would require petabytes of read data and months of compute time. The STAT approach accepts some loss of sensitivity in exchange for scale.

## Methods

### Screening modes

Two complementary query strategies target different segments of the public archive:

| Mode | Abbrev | Library organism queried | STAT retention gate |
|------|--------|--------------------------|---------------------|
| Microbe-as-library | MAL | PHI-base plant pathogen | Viridiplantae reads ≥ 1% (confirms in-planta origin) |
| Host-as-library | HAL | PHI-base plant host | ≥ 1 eukaryotic PHI-base pathogen ≥ 1% |

MAL targets pathogen-focused sequencing experiments (population genomics, disease surveys, transcriptome studies of isolated pathogens still in contact with host tissue). HAL targets host transcriptomics (resistance studies, field phenotyping, abiotic stress experiments) where the pathogen was not the study target. The modes are complementary and non-overlapping in their SRA accession sets.

Notably, `txid{n}[Host]` does not work as an SRA eSearch field — it is a free-text field, not taxonomy-linked. In-planta confirmation is therefore achieved via the STAT `Viridiplantae` percentage rather than metadata filtering.

### Reference database

`stat/stat_build.py` constructs a local reference database (`phibase_db.json`) from PHI-base (plant entries, eukaryotic pathogens only):

- **PHI-base** fungi + oomycetes: 205 seed pathogen taxids expanded to all descendant strains, formae speciales, and named variants via `ete3 NCBITaxa`. The expansion uses `intermediate_nodes=True`, which is essential — the default returns only leaf taxa, silently dropping species nodes that have named strains as children (affecting rust f. sp. taxa and many other plant pathogens).

The database contains 665 pathogen taxids (Fungi + Oomycota scope), a 21,352-entry name lookup table, and a 305,355-name Viridiplantae allowlist used to identify host species from STAT profiles without relying on `[Host]` metadata. Bacteria, nematodes, and viruses are excluded from the reference — detection is eukaryotic only (see Eukaryotic scope below).

### Co-infection detection

For each run passing the retention gate, `specific_hits()` in `stat/stat_filter.py` identifies leaf-level species detections by finding k-mer counts not nested under any more-specific count — returning species-diagnostic signal rather than genus-level aggregates. This is the core non-trivial algorithm; genus-level counts in STAT reflect LCA promotion of shared k-mers and are not informative for co-infection detection.

Detected species are classified by their interaction status relative to the study host:

| `interaction_status` | Meaning |
|----------------------|---------|
| `known` | Interaction recorded in PHI-base for this host–pathogen pair |
| `novel_host_range` | Pathogen known to PHI-base but not recorded on this host species |
| `novel_combination` | Both organisms known; combination not in PHI-base |
| `unresolved` | Taxid lookup failed; novelty cannot be assessed |

A `same_genus_secondary` flag marks cases where a detected secondary shares a genus with the primary organism. These are lower-confidence not because of cross-mapping — STAT's LCA design promotes shared k-mers to the genus node, preventing species-level bleed within a genus by construction — but because closely related species retain fewer unique diagnostic k-mers, pushing detections toward the threshold.

### Eukaryotic scope

Bacterial and viral co-detections are excluded. PolyA+ library selection systematically depletes prokaryotic mRNA (typically < 0.5% of reads in plant transcriptome data), making bacterial STAT percentages unreliable as co-infection indicators. Viral detections were also excluded: plant virus k-mer profiles lack the species-level discrimination needed to distinguish closely related strains at the signal-to-noise ratios present at low read fractions. The final detection scope covers eukaryotic pathogens only, with per-kingdom thresholds calibrated against pilot data:

| Kingdom | Detection threshold |
|---------|---------------------|
| Fungi | ≥ 0.5% STAT reads |
| Oomycota | ≥ 0.5% STAT reads |
| Nematoda | ≥ 1.0% STAT reads |

A LibrarySource pre-filter (`TRANSCRIPTOMIC`, `TRANSCRIPTOMIC SINGLE CELL`, `METATRANSCRIPTOMIC`, `VIRAL RNA`) removes GENOMIC and METAGENOMIC runs that appear in RNA-Seq strategy results due to submitter labelling errors (1,365 MAL and 6,573 HAL runs removed).

## Results

| Mode | Runs fetched | Non-RNA excluded | With STAT | Gate pass | Gate pass rate |
|------|-------------|-----------------|-----------|-----------|----------------|
| MAL  | 48,418      | 1,365           | 47,048    | 6,191     | 13.2%          |
| HAL  | 559,950     | 6,573           | 552,884   | 5,772     | 1.0%           |

The 607,197-entry STAT cache is the computational foundation of all downstream analyses. After deduplication to one run per BioSample (`biosample_representative == True`), 9,002 biological samples are retained across 1,285 BioProjects. Of these, 1,099 (12.2%) show at least one eukaryotic co-infection and 340 (3.8%) are high-confidence detections (secondary pathogen of a different genus). A total of 1,480 biosample-representative runs carry a `novel_host_range` classification — secondary pathogens present on host species not previously recorded in PHI-base.

The much lower HAL gate pass rate (1.0% vs 13.2% MAL) reflects the difference in study design: most host-focused libraries were collected under controlled conditions with a single inoculated pathogen, which reduces the incidence of field-acquired co-infections.

## Key output files

| File | Contents |
|------|----------|
| `stat/output/stat_build/data/phibase_db.json` | Reference DB: taxid maps, name allowlists, kingdom assignments |
| `stat/output/stat_fetch/data/stat_cache.jsonl` | Unified cache: RunInfo + STAT per run, append-only (607,197 entries) |
| `stat/output/stat_filter/data/runs.tsv` | 10,995 confirmed runs; one row per run with co-infection classification |
