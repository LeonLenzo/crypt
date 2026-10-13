# Reference scope: what goes into the database, and what is missing

This step selects the assemblies the Kraken2 database is built from and downloads their CDS.
`kraken_search.py` does pathogens, `kraken_hosts.py` does hosts. Both are CDS-only: these are
RNA-seq libraries, so an assembly with no gene model cannot be used, and `02_build`'s
`kraken_build.py` skips any that lacks one.

`scope_residents.py` is not part of the build. It measures what the current scope leaves out.

## The gap this step has, stated plainly

db_v3 was seeded from PHI-base and expanded to congeners, so it is a plant **pathogen**
database. But the module's own design rule is that breadth exists to give the classifier
**competitors**, not to cover more pathogens: Kraken2's LCA can only decline to resolve a read
when two or more database taxa share its k-mers.

A read from an organism with **no database representation at all** does not go unclassified.
It goes to the nearest taxon that is present.

That failure mode is invisible to every diagnostic in `05_filter`. The reads are real, the
organism is real, and only the assignment is wrong, so the accumulation curve rises, coverage
grows, and the score passes. No threshold can catch it, because nothing about the signal is
artefactual. This is a different problem from cross-mapping between database taxa, which is
what `score` was built for.

## Is core_nt the answer? No, and not because of size

core_nt is the obvious alternative, and the one Colvin & Chopra 2026
(doi:10.64898/2026.05.21.727028) used on two of our own BioProjects. Comparing the 2025-10-15
build's inspect table against db_v3's (`data/db_comparison.tsv`):

| | db_v3 | core_nt |
|---|---:|---:|
| total minimizers | 4.80 B | 59.40 B |
| species-rank nodes | 944 | 1,116,463 |
| **Fungi** | **4.05 B** (84.3% of db) | **2.93 B** (4.9% of db) |
| fungal species nodes | 877 | 124,718 |
| Oomycota | 185 M | 93 M |
| Viridiplantae | 568 M (CDS) | 5.66 B (genomic) |

**core_nt is 12.4x larger overall and holds 0.72x db_v3's fungal k-mer space.** Its fungal
breadth is nomenclatural, not genomic:

| per-species fungal ceiling | db_v3 | core_nt |
|---|---:|---:|
| median | 3,948,105 | **51** |
| species above 1 M minimizers | 768 (87.6%) | **634 (0.5%)** |

Half of core_nt's fungal species are one barcode sequence. It has **fewer** genome-scale fungal
species than db_v3, and only 221 of them are shared.

Two consequences:

- **Our coverage criterion is not portable to it.** At a median ceiling of 51 minimizers,
  "fraction of k-mer space seen" is not a stable quantity. This is worth stating generously
  about Colvin & Chopra rather than as a criticism: their read-count floor was forced by their
  database, not chosen over ours. It also independently justifies the fraction criterion, which
  previously rested on the *Phakopsora* argument alone.
- **The barcode entries are an active liability.** A conserved rDNA or ITS read matches hundreds
  of them and is promoted to genus or family by the LCA; Bracken then redistributes it back down
  to species on abundance priors. That is how core_nt plus Bracken manufactures species calls,
  arriving at the *Phakopsora* failure by a different route.

And the comparison above understates it, because db_v3 is CDS and core_nt is genomic. For
filamentous ascomycetes roughly half the genome is coding, so core_nt's transcript-compatible
fungal space is nearer 1.3 to 1.6 B against db_v3's 4.05 B, all of which is reachable by
RNA-seq reads.

Eight of the 32 pathogens `05_filter` detects on cereal hosts are fragments-only in core_nt
(`data/corent_vs_dbv3_species.tsv`), including *Pyrenophora teres* at 21,000 minimizers against
5.43 M, and *Fusarium culmorum*, *F. avenaceum* and *F. sporotrichioides* all near-absent while
*F. graminearum* sits at 3.6 M. On core_nt a mixed *Fusarium* infection has nowhere to go but
*F. graminearum*.

**Decision: core_nt is not adopted, and the 246 GB index is not downloaded.** Only its inspect
table is kept, in `_refdata/databases/kraken2_core_nt_20251015/`.

## What core_nt did tell us

It holds 413 genome-scale fungal species db_v3 lacks (`data/corent_only_genome_scale.tsv`).
Triaging those by class suggested they were mostly clinical, yeast and wood-decay fungi with
only ~25 relevant Dothideomycetes. **That triage was wrong**, and the abundance data below is
why: class-level reasoning is the wrong lens for this question.

## What is actually on a cereal leaf

Colvin & Chopra profiled 2,194 field samples and report per-sample fungal abundance. Their
tables are in `01_literature/studies/doi_10.64898_2026.05.21.727028/`, and
`scope_residents.py --residents` ranks them against db_v3 (`data/residents_colvin.tsv`).

**32% of maize and 39% of sorghum species-level fungal abundance sits on taxa absent from
db_v3.** The abundant absentees:

| species | maize | sorghum | prevalence |
|---|---:|---:|---:|
| *Plenodomus biglobosus* | 2.15% | **6.47%** | 96 / 100% |
| uncultured *Glomus* | 5.40% | 3.13% | 95 / 100% |
| *Cladosporium herbarum* | — | 2.34% | 89% |
| *Tilletiopsis washingtonensis* | — | 1.86% | 100% |
| *Filobasidium floriforme* | 1.30% | 1.58% | 82 / 100% |
| *Boeremia exigua* | 0.69% | 1.55% | 83 / 100% |
| *Cladosporium halotolerans* | 1.11% | 1.18% | 96 / 100% |
| *Kwoniella europaea* | 0.95% | — | 85% |
| *Dioszegia hungarica* | — | 0.76% | 100% |

*Plenodomus biglobosus* is sorghum's second most abundant species at 100% prevalence, and was
independently the largest single entry in the core_nt-only list at 21.6 M minimizers. Two
unrelated probes landed on the same organism. *Filobasidium* and *Kwoniella* are
Tremellomycetes, *Tilletiopsis* is Exobasidiomycetes and *Glomus* is Glomeromycetes, the three
classes the class-level triage dismissed.

Of the nine families Colvin names as consistently abundant across ITS and metagenomic surveys
as well as their own data, db_v3 holds six (`data/colvin_families.tsv`). It holds **zero**
Cladosporiaceae, Glomeraceae and Debaryomycetaceae. Aspergillaceae, which the core_nt triage
implied was a gap, is in fact covered at 155 species and 623 M minimizers.

### Colvin's table is not ground truth, and the artefact is visible in it

In **sorghum**, at 100% prevalence: *Parastagonospora nodorum* 2.86%, *Alternaria dauci* 1.43%,
*Plenodomus lingam* 1.24%, *Ascochyta rabiei* 1.22%, *Pyrenophora tritici-repentis* 1.15%. A
Nebraska sorghum leaf does not carry wheat, carrot, canola and chickpea pathogens at those
levels. Those are Pleosporales reads from unrepresented residents scattering across whichever
Pleosporales their database held. (*Malassezia restricta*, a human skin yeast, at 0.78% in maize
is the same story.) **Read it at genus and family level as "these lineages are present", never
as species truth.**

It also predicts what to expect from us: when the 1,120 literature-first sorghum runs go through
db_v3, Pleosporales pathogen calls should be inflated, because sorghum's dominant residents are
the ones we are missing. Sorghum currently has **zero** detections only because the kraken pass
ran on the 02_literature cohort.

## Wheat and rice, from the papers rather than from a guess

My first resident list was a guess from general phyllosphere literature, and it was wrong in
both directions: it missed the dominant taxa and padded with genera no cereal survey names.
Replaced 2026-10-13 by reading four papers in full. `RESIDENT_GENERA` in `scope_residents.py`
carries every entry with its source and evidence.

**Fischer 2025** (doi:10.3389/fpls.2025.1597807) is the authority for wheat: a meta-analysis of
33 studies, 190 locations, 924 species in 464 genera. Of its **13 core genera** (present across
all five wheat-growing regions, ranked by average abundance, lower is more abundant), **four are
absent from db_v3**:

| rank | genus | db_v3 |
|---:|---|---|
| 2.4 | *Alternaria* | 32 nodes |
| **4.75** | ***Cladosporium*** | **absent** |
| 5.65 | *Fusarium* | 113 |
| 7.52 | *Penicillium* | 48 |
| 7.55 | *Epicoccum* | 2 |
| 8.26 | *Aspergillus* | 154 |
| 8.36 | *Bipolaris* | 13 |
| **8.9** | ***Chaetomium*** | **absent** |
| **9.06** | ***Nigrospora*** | **absent** |
| 9.23 | *Stemphylium* | 2 |
| 9.33 | *Trichoderma* | 48 |
| **9.6** | ***Candida*** | **absent** |
| 9.66 | *Pyrenophora* | 9 |

*Cladosporium* is second only to *Alternaria* across the global wheat literature. And decisively:
**seven of the top 20 genera by average abundance rank are yeasts, and all seven are absent**
(*Sporobolomyces*, *Cryptococcus*, *Dioszegia*, *Filobasidium*, *Udeniomyces*, *Rhodotorula*,
*Vishniacozyma*).

**Allen 2026** (doi:10.1111/mec.70485, wheat flag leaves, ITS, 68 samples, two North Carolina
sites) corroborates independently: 13 of the 22 genera in its top-20 taxa figure are absent,
including *Cladosporium* (dominant at one site), *Cystofilobasidium*, *Aureobasidium*, *Mrakia*,
*Leucosporidium*, *Bullera*, *Selenophoma*, *Erratomyces*.

**Yin 2023** (doi:10.1128/spectrum.04371-22, rice leaves, ITS, 30 composite samples, Hainan):
**rice's single most abundant genus is *Papiliotrema* at 14.5%, and it is absent.** Then
unclassified Pleosporales 14.4% and *Symmetrospora* 5.7%. Among the 29 core fungal ASVs:
*Nigrospora* 13.2% (absent), *Pyrenophora* 13.1% (present), *Papiliotrema* 12.5%,
*Phaeosphaeria* 4.5%. Rice is thinner than wheat: it has no meta-analysis, and Dastogeer 2022
gives only random-forest discriminators, not abundances.

***Nigrospora* is the single highest-value target**: core in wheat (rank 9.06) AND rice (13.2%),
18 assemblies, zero gene models.

### Two naming artefacts and one taxonomic trap

- ***Magnaporthe* is not a gap.** It is the old name for *Pyricularia*, which db_v3 holds
  (*P. oryzae*, 109 cereal detections), and has zero assemblies under the old name.
- ***Phaeosphaeria* needs checking.** *P. nodorum* is a synonym of *Parastagonospora nodorum*,
  which db_v3 holds, so Yin's 4.5% may be the same organism under an old name.
- ***Cryptococcus*, *Candida* and *Rhodotorula* are wastebasket genera that were SPLIT.**
  *Vishniacozyma*, *Papiliotrema* and *Filobasidium* were all carved out of *Cryptococcus*, so
  "*Cryptococcus*" in an older ITS survey of the kind Fischer pools may mean exactly the genera
  listed separately above. Adding modern *Cryptococcus* sensu stricto, which is largely clinical,
  is NOT the same as covering what those surveys called *Cryptococcus*. Resolve against current
  NCBI taxonomy, not the paper's name.

## Annotation is the binding constraint, and the cost is 508 assemblies

CDS-only means an assembly is not enough. Of the 50 genera named across the four papers plus
Colvin, 15 are already in db_v3 and the other 35 split three ways
(`data/resident_annotation.tsv`):

| route | genera | assemblies | named species |
|---|---:|---:|---:|
| **add now** (annotation exists) | 22 | **465** | 107 |
| **miniprot** (assemblies, no gene models) | 8 | **43** | 25 |
| no assembly anywhere | 5 | 0 | 0 |
| **db_v4 addition** | **30** | **508** | **132** |

Against db_v3's 2,063 candidate assemblies and 877 fungal species that is **+25% assemblies and
+15% fungal species**. Expect ~95% to survive selection, by the db_v3 precedent (of 2,063
candidates, 1,986 passed BUSCO and 70 had no CDS).

**Why the count is assemblies and not a per-species quota:** db_v3's policy is a FLOOR WITH NO
CAP (`kraken_search.py`, `DEFAULT_FLOOR = 4`). Capping a deep species was measured to discard
real k-mer content, up to 62% of a species complex's observed content, so every annotated
assembly of a selected species is taken. That is safe for bias, because Kraken2 stores each
distinct minimizer once and redundant assemblies add no weight, but it is not free in build
time, disk or RAM. **If the total is too large, the lever is the genus list, never a cap.**

The weight is lopsided: five genera supply 354 of the 465 add-now assemblies.

| genus | assemblies | species | why it is big |
|---|---:|---:|---|
| *Aureobasidium* | 135 | 8 | pullulan production |
| *Cryptococcus* | 96 | 15 | clinical |
| *Candida* | 44 | 6 | clinical |
| *Rhodotorula* | 42 | 13 | lipid biotech |
| *Malassezia* | 37 | 16 | clinical (and a useful contamination marker: Colvin's maize 0.78%) |

Dropping the three clinical genera would save 177 assemblies, but *Candida* is a wheat core genus
and *Cryptococcus* a top-20 wheat yeast, so that trade has a real cost. The 8 miniprot genera
cost only 43 assemblies, so round 2 is cheap in build terms and expensive only in gene-prediction
time.

Round 2 is not optional: *Cystofilobasidium* and *Vishniacozyma* are wheat core and *Nigrospora*
is core in both wheat and rice. The 5 dead ends need a stated limitation, though only
*Erratomyces* and *Selenophoma* are genuine wheat-leaf losses; *Tomentella* and *Tylospora* are
ectomycorrhizal and probably soil carryover, like Colvin's uncultured *Glomus* at 5.4% of maize.

**Partial coverage is still protective.** A *Vishniacozyma* read does not need *Vishniacozyma*
in the database; it needs something close enough to promote to Filobasidiales rather than be
confidently mis-called as a Pleosporales pathogen. *Filobasidium floriforme* and *Papiliotrema
laurentii* are both annotated Filobasidiales, so the order gets a home even where the genus does
not. The goal is a nearer home, not a correct one.

## What this test can and cannot do

It can only move the co-infection rate **down**. If unrepresented residents are inflating
pathogen calls, adding them reduces detections; there is no version where the rate rises. The
current wheat 47% is therefore an upper bound until this is settled. That is a reason to do it
ourselves rather than be asked.

## A trap in `_classified.tsv`

3,382 of its 24,467 rows (13.8%) carry `score = nan`, all `cls = unclassified`: `detect.py`
could not fit a curve for that host x species cloud. **Test the floors positively (`>=`), never
by negating (`not <`).** NaN compares False against everything, so `if score < FLOOR: continue`
admits those rows instead of skipping them. It inflated the cereal species count here from 32 to
44 before it was caught, pulling in *Bipolaris oryzae*, *Claviceps purpurea* and *Aspergillus
flavus* on coverage alone with no score at all; 393 of the NaN rows clear the 1% coverage floor.
`05_filter/figures/prep_network.py` tests positively and is correct, so the published
co-infection rates are unaffected.

## Layout

```
kraken_search.py      select pathogen assemblies + download CDS (the build's input)
kraken_hosts.py       resolve + download host CDS
scope_residents.py    measure what the scope misses; writes the tables below
data/
  db_comparison.tsv              db_v3 vs core_nt, structural
  corent_vs_dbv3_species.tsv     per-species ceilings for the 32 cereal detections
  corent_only_genome_scale.tsv   413 genome-scale fungal species db_v3 lacks
  residents_colvin.tsv           Colvin taxa ranked by abundance, with db_v3 membership
  colvin_families.tsv            the nine dominant families against db_v3
  resident_annotation.tsv        NCBI assembly/annotation counts per resident genus
```

Run from the repository root:

```
python 03_kraken/01_search/scope_residents.py --compare --residents
python 03_kraken/01_search/scope_residents.py --annotation      # network, ~2 min
```

Inputs it needs that are not in the repository: `_refdata/databases/kraken2_core_nt_20251015/
inspect.txt` (see its `PROVENANCE.md`) and Colvin's supplements in `01_literature/studies/`
(gitignored; URLs in that directory's `SOURCE.txt`).
