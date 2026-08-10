# crypt — Cryptic Co-infection Mining from Public Plant RNA-seq

*Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)*

## Background

Plant disease studies deposited in the NCBI Sequence Read Archive (SRA) are designed around a single target pathogen. Field-collected samples, however, routinely harbour additional co-infecting organisms whose signal is present in the sequencing data but goes undetected and unreported under single-target study designs. We hypothesise that a substantial fraction of the public plant RNA-seq archive contains secondary pathogen signal sufficient for taxonomic detection — representing an untapped resource for co-infection epidemiology at landscape scale.

This project mines public SRA data for evidence of unreported co-infection using NCBI STAT pre-computed k-mer taxonomy profiles, requiring no raw read download or re-alignment. Two complementary query strategies — one targeting pathogen-focused libraries, one targeting host-focused libraries — together screen the plant-associated RNA-seq corpus for eukaryotic secondary pathogens (fungi, oomycetes, nematodes).

## Modules

The pipeline is organised into three sequential modules, each with its own rationale, methods, findings, and limitations documented in the respective module README.

| Module | Purpose | README |
|--------|---------|--------|
| **[stat/](stat/)** | STAT k-mer screening of 608,368 SRA runs; primary co-infection detection | [stat/README.md](stat/README.md) |
| **[metadata/](metadata/)** | BioProject/BioSample enrichment, literature linkage, and LLM study design classification | [metadata/README.md](metadata/README.md) |
| **[kraken/](kraken/)** | Orthogonal Kraken2 species-level validation of STAT detections | [kraken/README.md](kraken/README.md) |

## Headline results

Screening 608,368 SRA runs (46,315 MAL + 546,816 HAL after non-RNA exclusions), STAT k-mer profiling yielded 10,995 confirmed runs across 1,285 BioProjects. After deduplication to one run per biological sample, 9,002 biosample-representative runs were retained. Of these, 1,099 (12.2%) showed evidence of at least one eukaryotic co-infection; 340 (3.8%) were high-confidence detections where the secondary pathogen belongs to a different genus from the primary. A total of 1,480 biosample-representative runs showed novel host range interactions — secondary pathogens detected on host species not previously recorded in PHI-base for that pathogen.

LLM-based study design classification of 1,285 BioProjects revealed a marked setting effect: field-collected samples showed a 22.7% co-infection rate versus 11.4% in controlled laboratory conditions, consistent with the ecological complexity of field environments. The majority of BioProjects were single-pathogen studies (895/1,285), reinforcing that co-infection is the rule rather than the exception under field conditions yet remains systematically understudied.

## Scope

**Host scope:** Viridiplantae (plant hosts) only, anchored by PHI-base plant–pathogen interaction records and the ICTV plant virus master species list.

**Pathogen scope:** Eukaryotic pathogens only (fungi, oomycetes, nematodes). Bacterial co-detections are excluded — PolyA+ library selection systematically depletes bacterial mRNA, making bacterial STAT percentages unreliable indicators of co-infection. Viral detections were also excluded as STAT's k-mer approach lacks the specificity to discriminate closely related plant virus strains at the thresholds used here.

## References

- **PHI-base:** Urban et al. (2020) *Nucleic Acids Res* — [phi-base.org](https://phi-base.org)
- **ICTV VMR:** [ictv.global/vmr](https://ictv.global/vmr/current)
- **NCBI STAT:** Katz et al. (2021) *J Bioinform Comput Biol* — [PMC8450716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8450716/)
- **Kraken2:** Wood et al. (2019) *Genome Biol* — [PMC6883579](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6883579/)
