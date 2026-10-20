# crypt — Cryptic Co-infection Mining from Public Plant RNA-seq

*Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)*

## Background

Plant disease studies deposited in the NCBI Sequence Read Archive (SRA) are designed around a single target pathogen. Field-collected samples, however, routinely harbour additional co-infecting organisms whose signal is present in the sequencing data but goes undetected and unreported under single-target study designs. We hypothesise that a substantial fraction of the public plant RNA-seq archive contains secondary pathogen signal sufficient for taxonomic detection — representing an untapped resource for co-infection epidemiology at landscape scale.

This project mines public SRA data for evidence of unreported co-infection using NCBI STAT pre-computed k-mer taxonomy profiles, requiring no raw read download or re-alignment. Two complementary query strategies — one targeting pathogen-focused libraries, one targeting host-focused libraries — together screen the plant-associated RNA-seq corpus for eukaryotic secondary pathogens (fungi, oomycetes, nematodes).

## Modules

The pipeline has **two live modules**, each with its own rationale, methods, findings and
limitations documented in its README. `01_literature` builds the cohort from the papers;
`02_classify` builds the reference database, fetches the reads, classifies them and calls
detections.

Four earlier modules were retired on 2026-10-20 and moved to `superseded/`, which is on disk and
out of the repository. [SUPERSEDED.md](SUPERSEDED.md) records what each did, what replaced it,
and which of its outputs the live code still reads. The short version: the pipeline originally
screened SRA for pathogen signal and then read the papers (`01_stat`, `02_literature`,
`03_kraken`, `04_kallisto`). It now reads the papers first, because a signal-gated denominator
cannot support a prevalence estimate: the STAT gate's pass rate is host-dependent (wheat
3.42%, soybean 0.55%).

Module and step directories are numbered in workflow order, so the tree reads as the pipeline: a
numbered directory is a step, an unnumbered one (`slurm/`, `utilities/`, `report/`, `data/`,
`logs/`) is support.

| Module | Evidence | Purpose | README |
|--------|----------|---------|--------|
| **[01_literature/](01_literature/)** | the papers, as the frame | construct a cohort from published field RNA-seq studies rather than from detector hits | [01_literature/README.md](01_literature/README.md) |
| **[02_classify/](02_classify/)** | the reads, by k-mer | db_v4, read fetching, Kraken2 classification, detection calls | [02_classify/README.md](02_classify/README.md) |

```
01_literature/  seeds  data  curate  retrieve  report  release  studies  controls
01_literature/ prompts  data
02_classify/    01_refs  02_build  03_cohort  04_classify  05_detect
                seeds/  slurm/  report/  utilities/
superseded/     01_stat  02_literature  03_kraken  04_kallisto    RETIRED, untracked
```

## Headline results

**Status, 2026-10-20.** The results below were produced against **db_v3** on a STAT-gated set of
runs. db_v4 is built and validated, and 1,184 runs have been reclassified against it, but that
is 15% of the 7,688-run literature cohort and those 1,184 are exactly the old STAT-gated set, so
**no rate on this page is a prevalence estimate.** The db_v4 reclassification of the full cohort
is the active front. See `CLAUDE.md`, "What is NOT quotable yet", for the four selections that
sit in the current denominator.

Screening 608,368 SRA runs (46,315 MAL + 546,816 HAL after non-RNA exclusions), STAT k-mer profiling yielded 10,995 confirmed runs across 1,285 BioProjects. After deduplication to one run per biological sample, 9,002 biosample-representative runs were retained. Of these, 1,099 (12.2%) showed evidence of at least one eukaryotic co-infection; 340 (3.8%) were high-confidence detections where the secondary pathogen belongs to a different genus from the primary. A total of 1,480 biosample-representative runs showed novel host range interactions — secondary pathogens detected on host species not previously recorded in PHI-base for that pathogen.

LLM-based study design classification (full-text-gated: 732/1,285 BioProjects with retrievable manuscript text, 6,467 BioSamples) revealed a marked setting effect: field-collected samples showed an 11.1% biotic-only cryptic co-infection rate versus 4.0% in greenhouse and 7.8% in other controlled conditions — roughly 2.8x higher in the field than in the greenhouse, consistent with the ecological complexity of unconstrained wild pathogen exposure. The majority of classified BioSamples were single-pathogen-focus studies, reinforcing that co-infection is the rule rather than the exception under field conditions yet remains systematically understudied.

**This setting contrast is provisional.** It rests on `llm_study_setting == "field"` identifying genuinely field-collected samples, and hand-checking shows it does not reliably do so: a *Botrytis cinerea* inoculation atlas (PRJNA1217477, 211 samples), *in vitro* work on artificial surfaces (PRJNA526829, 60) and controlled compatible/incompatible inoculations (PRJNA328045, 39) are all classified `field`. These three are now registered in [`_exclusions.py`](_exclusions.py) and dropped at the analysis stage, which takes the field cohort from 3,282 BioSamples to 2,972. A title screen flags roughly 13% of the cohort as plausibly non-field, so the registry is a record of what has been checked, not a complete list. The same conflation runs through the other provenance fields: `geo_loc_name` records the submitting institution as often as the collection site ("USA: California, Davis" is UC Davis, not a paddock), and `collection_date` ranges back to 1925 for culture-collection isolates. The per-sample provenance verification that separates these axes is designed but not built, so the field rate and the 2,643-sample denominator should not be quoted as settled. See [SUPERSEDED.md](SUPERSEDED.md) and [01_literature/README.md](01_literature/README.md): the per-sample provenance verification named here WAS subsequently built, in `01_literature`, which reached host and tissue 100%, location 99.8% and year 97.2% on the 7,688-run cohort.

Multi-strategy literature resolution linked 58.0% of BioProjects (746/1,286) to a primary publication (PMID or DOI), with full-text methods sections retrieved for 57.0% via PMC/Unpaywall/manual PDF fallback. The remainder are principally data-only submissions and unpublished surveillance datasets.

The STAT co-infection rates above are the screen's own estimate and are superseded by the Kraken2 pass below. STAT resolves k-mers shared between close relatives to their lowest common ancestor, which costs species-level resolution in exactly the taxonomically dense clades co-infection work depends on: a kallisto pilot confirmed STAT reports 0% eukaryotic signal for all 15 runs dominated by *Puccinia striiformis* f. sp. *tritici*, while Kraken2 and kallisto both place it at 65 to 68%. The screening funnel above stands; the co-infection rate derived from it does not.

**Kraken2 species-level detection.** The target is all 2,719 field/aerial BioSamples from the LLM-classified set, deliberately not narrowed to already-flagged co-infections, since the point is catching what STAT missed. The current database (`db_v3`) holds 2,115 BUSCO-screened assemblies (2,020 pathogen, 95 host) spanning 1,043 distinct taxa, against db_v2's 326. Breadth was the fix for a specific failure: with only 326 taxa most reads had exactly one plausible match and received a confident species call whether or not the species was present, which is why soybean rust was "detected" in every wheat sample. Kraken2's LCA can only decline to resolve a read when two or more database taxa compete for its k-mers, so the rebuild added competitors rather than additional nameable pathogens. Measured on three wheat runs against identical parameters, db_v3 cut the soybean-rust artefact by 65 to 72%, brought host reads from 0 to 2.3-3.0% (there was no *Triticum* in db_v2 at all), reduced unclassified reads by 11 to 17 percentage points, and left the true target unchanged. 3,220 of 3,225 runs are classified; the 5 failures are known-failed downloads.

BUSCO completeness screening is retained as a QC record rather than a filter. Of the 2,063 candidate assemblies scanned, 1,986 pass, 7 fall below the completeness bar and 70 have no CDS at all; median completeness among the 1,993 scored is 98.4%. All 1,986 passing assemblies are selected, and a per-taxid fallback adds 34 more that failed but are their taxon's only representative, so the screen cannot remove a taxon. It prunes surplus assemblies within a taxon and nothing else, and the attrition that actually matters is missing annotation, ten times more common than failing completeness.

**Detection filtering and co-occurrence.** A raw read count cannot separate a real infection from a classifier pile-up, so each (host x species) detection is scored on the shape of its k-mer accumulation curve: a real infection accumulates new distinct k-mers as sequencing deepens, an artefact stacks reads onto the same fragment. Two thresholds then gate what counts as a co-infection, and they answer different questions. `score >= 0.7` is the artefact filter, keeping ~94% of known-real detections and removing 100% of the known off-host artefacts. `coverage >= 1%` is a presence floor, meaning the organism's DNA is in the sample. The resulting rates are 47% of wheat samples co-infected and 39% of maize. Co-occurrence edges are tested against a prevalence null (Fisher's exact, BH-adjusted p < 0.05) rather than counted raw, because two common pathogens share a sample by chance; only wheat, maize and barley carry three or more significant edges. These are **db_v3 numbers and are superseded as a measurement**, though the reasoning stands. db_v4 replaced the fitted `score` with a depth-normalised saturation statistic, and coverage is not comparable across database versions because its denominator is each build's own inspect table. The current derivation, including why the 1% floor now has a basis it did not have under db_v3, is in [02_classify/05_detect/README.md](02_classify/05_detect/README.md).

**What is not yet resolved.** DNA presence is not proof of active co-infection, and coverage cannot separate the two: spore drift, surface inoculum and mixed samples all deposit real organism DNA. The clearest case is the residual soybean rust on wheat, which read-level alignment confirms is genuinely *Phakopsora*-like sequence rather than misclassified stripe rust (of 248,504 reads, 94% map to *Phakopsora* and 18 to *P. striiformis*), yet *P. pachyrhizi* cannot infect wheat. This limitation is inherent to SRA-based co-infection inference and is disclosed rather than thresholded away. Resolving it is the job of competitive read assignment, which asks which genes a detection's reads hit and whether a call survives when congeners and the host compete for the same reads. A first version (`superseded/04_kallisto`) pseudo-aligned a stratified set of 174 runs against per-host kallisto indices with no failures, but never produced a full `calls.tsv`, so no detection has been confirmed or rejected by competition yet. It is being rebuilt as `02_classify/04_align` against db_v4's reference set, and its target set is now defined rather than guessed: the 1,097 low-coverage high-saturation calls, plus the congener clusters.

A second, sharper case for it emerged on the first db_v4 run. One run calls fourteen *Alternaria* species, every one at coverage 5.6 to 7.3% and saturation 0.070 to 0.080. Fourteen genuinely co-infecting congeners would differ wildly in abundance; one organism's reads spread across fourteen near-identical references looks exactly like this, and the uniformity is the diagnostic. No per-detection threshold can catch it, because each call is individually well-supported. It is a property of the database's congener density, which breadth deliberately increased.

## Scope

**Host scope:** Viridiplantae (plant hosts) only, anchored by PHI-base plant–pathogen interaction records and the ICTV plant virus master species list.

**Pathogen scope: fungi and oomycetes only.** Bacteria and viruses are excluded for the same underlying reason before any classifier is involved: polyA+ selection, which most of this cohort used, removes most bacterial and many viral reads, so their percentages do not measure what was in the tissue. On top of that, STAT's k-mers cannot discriminate closely related plant virus strains at these thresholds, and viral genomes are too small and fragmented for the Kraken2 pass. Nematodes are excluded because annotated assemblies are absent for most PHI-base seeds, so there is nothing to build a reference from; the Nematoda entry left in the STAT gate is inert, admitting 2 of 10,995 retained runs and none on nematode signal alone.

## References

- **PHI-base:** Urban et al. (2020) *Nucleic Acids Res* — [phi-base.org](https://phi-base.org)
- **ICTV VMR:** [ictv.global/vmr](https://ictv.global/vmr/current)
- **NCBI STAT:** Katz et al. (2021) *J Bioinform Comput Biol* — [PMC8450716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8450716/)
- **Kraken2:** Wood et al. (2019) *Genome Biol* — [PMC6883579](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6883579/)
