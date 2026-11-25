# 03_align: which genes the reads come from, and whether they are pathogenicity genes

Sketched 2026-11-02, revised 2026-11-17 and 2026-11-25. Nothing here is built. This file is the design and
its reasoning; **the deliverable is `release/04_align`**, which does not exist yet. It
supersedes the `02_classify/04_align` named in older documents, which could never have had that
path because `02_classify/04_` is `04_classify`. Align is its own module for the reason
CLAUDE.md gives: a different method answering a different question, and nesting it under the
classifier inverted the dependency.

**Revised 2026-11-17 by the methods audit**, on three points: the cohort is now wheat, the
reference set moves from db_v4 to the db0/db1/db2 ladder, and the host-response question moves
into this module as a third job rather than sitting unowned in `TODO-RESULTS.md`. Everything
about the quantifier and the statistic is unchanged. The measured-scope table near the end is
flagged as stale and says why.

**FOUR DECISIONS TAKEN 2026-11-25**, each recorded in full in the section it belongs to and
summarised here because they were what blocked a build. Nothing is built yet.

1. **The host reference is Chinese Spring and nothing else**, Ensembl Plants
   `triticum_aestivum_refseqv2` (GCA_018294505.1, IWGSC RefSeq v2.1), already on disk from
   `02_reference`. The align index must NOT inherit db2's host sink, which is 19 wheat
   cultivars. See "The host reference".
2. **Question 3 survives as a within-site blocked contrast, and it is the question the chapter
   wants.** Site becomes the blocking factor rather than a result and geography is a nuisance
   variable rather than an axis. 29 informative blocks, 260 runs, median 8 runs per block,
   measured rather than assumed by `blocks.py` in this directory. See "The design, and the
   constraint that drives it". Closes open question 0.
3. **The host is quantified at gene level with homoeologues collapsed, and the collapse is
   derived from our own alignments**, not from a published triad table. See "What the index
   already gives it". Closes open question 6.
4. **Gene set annotation lives in this module.** holobase is not part of this chapter, Leon's
   call. See "Gene sets". Closes open question 2.

## Why this module is now the spine, not an appendix

Leon's call, 2026-11-02. The prevalence story is **demoted from the claim to the instrument**.
The cohort cannot support a rate: the unselected arm rests on 24 site-years, and the confound
is three-way rather than two, because the disease-selected data are global surveys (median one
run per site-year) and the unselected data are single-field trials (median seven). Those are
different kinds of study being contrasted, which no model repairs. A heavily caveated
percentage is also not a claim anyone needs.

So `02_classify` becomes the SCREEN: it says which organisms are present and worth
characterising, which is what it is good at and is fully defensible. This module carries the
finding. The reasoning is that transcription is the one thing RNA buys that amplicon surveys
and metagenomes cannot, and for an obligate biotroph it is close to an activity assay.

It also fixes the thesis arc. "Co-infections occur at rate X" leads nowhere; "unreported
co-infections are transcribing pathogenicity programmes in field samples" leads straight into
Chapter 2's in planta mechanism.

**Chapter 2 is this module's ground truth, not a crossover risk.** `phd/02-coinfection-rnaseq`
is dual RNA-seq of wheat with *P. nodorum* and *P. tritici-repentis*. Its in-vitro half is
COMPLETE and already runs this exact stack: `txi_*.rds` is tximport over pseudoalignment, then
DESeq2, then `GO_enrichment_nodorum.csv` and `GO_enrichment_repentis.csv`. The in-planta design
is 72 libraries with mocks, single infections, n = 3 and two timepoints. So the no-contrast
statistic below can be run on libraries where the species present are known and the correct
answer has already been computed with replicates. That is a validation Chapter 1 cannot
generate for itself.

## Prior art, Undermind 2026-11-02

Search `functional enrichment per taxon mixed fungal transcriptomes`, 134 papers. The approach
is not new in kind and must not be written as if it were.

**The direct precedent is Morales-Cruz et al., closed-reference metatranscriptomics of the
grapevine trunk-disease complex.** Field wood samples, several co-occurring fungi in one
library, each taxon's putative virulence transcripts quantified including cell wall degradation
and secondary metabolism. It also optimised its multi-species mapping by simulation, which is
the cross-mapping problem below, handled rather than ignored. Indexed in the workspace as the
2017 bioRxiv preprint; find the journal version before citing it. Massonnet et al. on the same
system add a warning that bears on our reference set: a single isolate misses variable
secondary-metabolite clusters, so per-species reference breadth is not optional.

What that work did NOT do is turn the profile into a calibrated test. It is descriptive
quantification of virulence transcripts, not a statistic with a null. **That is the gap this
module occupies**, together with three narrower ones: foliar cereal tissue rather than wood,
congeneric resolution (the grapewine complex is not six *Puccinia* species), and a stated
detection floor rather than a reported abundance.

Adjacent and worth reading: Franzosa et al. for species-level functional profiling of
metatranscriptomes, which is how the bacterial world solved the same shape of problem;
Gutiérrez & Hablützel for taxonomically resolved pathway profiling in eukaryotic biomes;
Marcelino et al. for mixed-community fungal metatranscriptomics as a proof of concept; and Xu
et al., the closest ecological analogue, functional dynamics of a phyllosphere mycobiome under
rust infection.

Tooling exists and should not be rebuilt: **singscore** for within-sample rank scoring,
**FungiFun3** and **FunFEA** for fungal gene set enrichment, **GOAT** for robust set
enrichment, and **Crossmapper** for estimating cross-mapping under a chosen reference set and
read length, which is exactly the measurement the trap below needs.

### What the literature says we may conclude

This is the part that constrains the writing, not the code.

- Fungal RNA supports transcriptional **activity**, not pathogenic causation. Urquhart & Idnurm
  show directly that high expression need not imply a gene is required for disease.
- **CAZyme programmes overlap between parasitism and saprotrophy** (Olson et al.), and some
  pathogens switch to saprotrophic growth on plant residue (Gay et al.). So cell wall
  degrading enrichment on its own does not separate infection from decay, which is what the
  claim ceiling below already says.
- **The obligate biotroph asymmetry has support.** Rust expression in infected tissue is read
  as active host-dependent growth precisely because rusts cannot grow on dead tissue (Tao
  et al.; Bruce et al. on the *Puccinia triticina* effective secretome). Still evidence of
  activity, not proof a given gene drives disease.
- A within-sample score is **relative prominence inside that sample**. Not induction, not
  differential expression, and not evidence that the expressed fraction is unexpectedly high.

## The three questions, in order

`03_classify` establishes that an organism's DNA is present above three floors. It cannot say
which genes the reads came from, which makes three different things impossible:

1. **Is this organism actually here, or are its reads a congener's?** The old `04_kallisto`
   module was built to answer exactly this and its design is carried forward below.
2. **Is this organism expressing a pathogenic lifestyle?** New, and the reason the module is
   being revived now rather than when it was first specified. Leon 2026-11-02: we accept that
   causation of disease is not available from archived field samples, but "this fungus is
   expressing genes associated with pathogenicity" is a weaker claim the data can carry.
3. **Does the HOST respond differently when more than one pathogen is present?** Leon's third
   question for the paper, moved here on 2026-11-17 from `research/TODO-RESULTS.md` section
   1b, where it had no owning module. It belongs here and not in its own module because it
   needs the same artefact as the other two and no other: a per-gene quantification against an
   index that already holds host CDS. Nothing else in the pipeline produces one.

The order is not negotiable. A functional test run on a cross-mapping artefact characterises
the function of an organism that is not there, and a host-response contrast built on
mis-assigned co-infection calls contrasts the wrong samples. **Confirm presence, then
characterise.**

**Question 3 is structurally different from questions 1 and 2 and must not inherit their
design.** The whole argument below for a within-sample, no-contrast statistic exists because
there is no uninfected control for a fungus found opportunistically in somebody else's
experiment. For the host that is not true: every run has a host, singly and multiply infected
runs both exist in the same deposits, and the contrast is between runs rather than within one.
So question 3 CAN use a conventional between-sample model, and should, with the design
constraints in "The host-response design" below. Carrying the unique-reads rule across to it is
a separate judgement, also below, and the answer is probably not.

## What carries forward from 04_kallisto, unchanged

Its README holds the design; `SUPERSEDED.md` holds what replaced what. Three things survive
intact because they were right and were measured rather than assumed:

- **One index per host**, holding host CDS plus every species Kraken2 saw in that host at any
  level, plus the measured k-mer neighbours of the flagged species. Host CDS are in the index
  so host reads are absorbed rather than forced onto a pathogen. A genuinely present species
  missing from the index has its reads pushed onto whatever relative is present, which inflates
  the very thing the module tests, so breadth here is cheap insurance.
- **Tagged transcript headers**, `{species_slug}|{accession}|{transcript_id}`. This is what
  makes the functional layer possible at all: per-target abundances aggregate to species while
  keeping a gene-level identifier. The scheme was designed for breadth and happens to be
  exactly what a gene-set test needs.
- **`02_build/data/neighbours.tsv`**, 398 pairs at >=1% containment over 45 species, computed
  from cached FracMinHash sketches at k=31 rather than from taxonomy. Containment, not
  phylogenetic rank, predicts cross-assignment: at congeneric divergence synonymous
  substitutions break a 31-mer every few bases, so rank is the wrong unit. It already produced
  a usable negative result, that the 470 low-level *Ascochyta* detections on wheat share less
  k-mer space with *P. nodorum* (0.13%) than with *Alternaria alternata* (0.16%) and are
  generic Dothideomycete background rather than bleed.

What does NOT carry forward is the reference set. **It must be whichever database tier the
detections being confirmed were called against**, which as of 2026-11-17 is the db0/db1/db2
ladder from `release/02_reference`, not db_v4. Coverage and containment both divide by the
build's own inspect table, so mixing tiers measures two things against different denominators
and nothing is comparable. Pin the tier in the module's own output and state it beside every
number. The 398-pair `neighbours.tsv` above was computed on db_v4's species set and is a
MEASUREMENT CARRIED ACROSS, not a denominator; it needs recomputation on the shipped tier
before it is used as one.

**The recomputation already has a tool, 2026-11-25, and it is not a sketch pipeline.**
`release/02_reference/figures/crossover/compute_matrix.py` computes exact pairwise minimizer
sharing from the shipped CDS at the build's own `l=31 / k=35`, runs on Setonix, and has run:
254 wheat-row taxa, block-diagonal by genus at mean 0.015 within a genus against 0.0001
across. That is what `neighbours.tsv` measures, on the current build and at the right k, so
the route is to re-run that script over the taxa this module needs rather than to rebuild
FracMinHash sketching. **Its present scope is not enough.** Of the 90 non-host taxa `db2`
calls, 39 are in the matrix and 51 are not, and the 51 are the *Alternaria* congeners and the
resident core, which is exactly where the presence question bites.

## The quantifier is an ALIGNER, and there is no EM

Decided 2026-11-02 after working out that the EM was the source of the bias below, not an
incidental detail. The old module used kallisto; this one uses **bowtie2** against the combined
tagged CDS.

**Why not kallisto or salmon.** Pseudoalignment judges a read by exact k-mer compatibility,
which is binary. An aligner scores mismatches, which is graded, and that distinction is
load-bearing here: conserved PROTEIN is not conserved NUCLEOTIDE, synonymous sites diverge, and
`neighbours.tsv` exists precisely because synonymous substitutions break a 31-mer every few
bases at congeneric divergence. k-mer methods discard that signal; an aligner reads it. Salmon's
selective alignment plus decoys would be the middle road and **is not available**: Setonix
modules carry bowtie2, bwa, bwa-mem2, hisat2, STAR and samtools, with no salmon and no
minimap2.

**Why bowtie2 specifically.** The references are CDS, so there are no introns to splice across
and a splice-aware aligner (STAR, hisat2) is overhead plus, for STAR, a large RAM footprint.
bwa-mem is a reasonable substitute. Host CDS stays in the reference as a read sink, which is the
role salmon's decoys would otherwise have played.

**Report multiple alignments and classify, do not filter on MAPQ.** bowtie2 needs `-k` so that
ambiguity is visible at all, and in `-k` mode its MAPQ no longer carries the usual meaning, so
a MAPQ filter would silently do something other than what it looks like. Classify instead by
inspecting the species tags of the reported alignments.

### The statistic: species-unique reads, the same rule on both sides

For each read, take the species tags of its alignments at or near the best score. If that set
holds exactly one species, the read is species-unique and counts for that species and gene.
Otherwise it is ambiguous and is **discarded**. Both the gene set and the background are then
built from unique reads only, so the deflation that the EM applied to one side of the ratio
cannot arise: there is no apportionment step.

This costs power, because every conserved gene leaves the background. Power is recoverable by
scoping to deep libraries. A statistic that cannot manufacture its own result is worth more than
a powerful one that can.

**`near the best score` is a real parameter and must be exposed.** A read two mismatches from
species A and three from B should probably count for A; a hard cutoff is arbitrary. Require the
best alignment to beat the runner-up by a stated margin in alignment score, and report the
statistic's sensitivity to it. Do not bury the margin as a constant.

### Positional evenness, which the EM could not give at all

The BAM says WHERE in the transcript reads fall, so a gene genuinely expressed (reads spread
along the CDS) separates from a conserved domain pulling traffic (reads piled in one window).
That is a direct control on the failure mode below, and it is the same logic `05_detect` already
applies one level up as "reads piled, not spread". Compute spread per (species x gene) and
require it before a gene joins either set.

### Presence, before function

The three facts the old module defined still stand, read off alignments rather than an EM:

- **species-unique reads.** Real presence is far above zero; pure cross-mapping is near it.
- **gene breadth**, distinct genes carrying unique reads. Real presence spreads over many
  genes; cross-mapping concentrates on a few conserved ones.
- **competitive survival**, counts with neighbours in the reference against alone. A collapse
  when competitors are present marks the detection as cross-mapping.

There is precedent for this in the project already: the *Phakopsora* adjudication recorded in
CLAUDE.md realigned 248,504 reads, found 94% mapping and only 18 going to *P. striiformis*.
Read-level alignment is how this project has settled cross-mapping arguments before, and
*Melampsora* at ~0.9% is still flagged as needing the same check.

## Function: the part that is new

### The test, and why not differential expression

There is no contrast available and there never will be: no replicates, no uninfected control,
no time course, because every library is an opportunistic field sample from somebody else's
study. So the test is **within one organism, within one run**, which has the incidental benefit
of cancelling library size, host fraction and the organism's absolute abundance. Those three
confounds have obstructed every cross-sample comparison in this project.

### Two statistics, and the recommendation changed after the literature search

Both statistics below run on **species-unique read counts**, per the quantifier section: the
same uniqueness rule applies to the gene set and to the background, and nothing is apportioned.

**Primary, within-sample rank scoring.** Rank the transcripts assigned to one organism in one
run by abundance and score the gene set against that ranking (singscore, Foroutan et al., which
was built for the no-contrast case; FungiFun3 or FunFEA supply the fungal set definitions, GOAT
the enrichment statistic). Published, widely used, and it makes a bounded claim: the set is
prominent within this sample.

**Secondary and DESCRIPTIVE ONLY, the genomic repertoire.** The expressed fraction of an
organism's effector or CAZyme repertoire against the expressed fraction of its gene content
overall, as a 2x2 per (run x species x gene set). I proposed this as the primary test on
2026-11-02 before searching, on the reasoning that the organism's own genome is an internal
null needing nothing borrowed from another sample. **The literature does not support treating
it as a calibrated test**: a within-sample score is relative prominence and is not evidence
that a genome-wide expressed fraction is unexpectedly high, and detection of a transcript as
"expressed" is itself depth-dependent, so the fraction moves with library depth rather than
with biology. It stays as a descriptive statistic and a sensitivity check. It is not the
headline number.

### The trap this design exists to avoid, and the residue that survives it

Worth stating in full, because it is the reason for every choice above and it will be the first
thing a reviewer probes.

Effectors are among the fastest-evolving, most lineage-specific sequences in a fungal genome,
so an effector read usually matches exactly one reference transcript. Housekeeping and
ribosomal transcripts are conserved by selection, so their reads match several. Under a
competitive EM the two halves of the ratio are then treated differently by the same algorithm:
the effector numerator is unique, full weight, intact, while the background denominator is
ambiguous and gets apportioned away to competitors. Deflate the denominator and leave the
numerator alone and the ratio rises for reasons that have nothing to do with the organism.

**And the bias is worst exactly where this study lives.** The EM gives ambiguous reads to
whichever species has more unique support. An abundant declared *P. striiformis* alongside a
low-level undeclared *P. triticina* takes nearly all the shared housekeeping; the triticina
keeps its own unique effector reads because nothing competes for them. Its apparent
transcriptome becomes mostly effectors. That species is the undeclared secondary co-infection,
so the bias concentrates on the target population and grows as the evidence weakens. It is not
a nuisance at the margin.

**Unique-reads-only removes the mechanism**, because nothing is apportioned. What survives:

- **The irreducible core.** Where two references are identical over a read length, no aligner
  helps; that is information-theoretic. Those reads are discarded from BOTH sides, which is the
  honest handling, and alignment lets uniqueness be computed PER BASE rather than per
  transcript, which is strictly finer than an equivalence class.
- **Differential power, not differential bias.** A species with few diagnostic positions
  yields few unique reads on both sides, so its test is underpowered rather than inflated.
  Report the unique-read count beside every result so an underpowered test is visible.
- **Selection into the test.** A species reached `05_detect` by covering db_v4 minimizer space,
  possibly via its diagnostic regions, so entry is already correlated with uniqueness. The
  off-host null is what bounds this.

**Three controls, and they go in the FIRST pass, not the second.** The magnitude of what
remains is unknown and must be measured rather than argued, but the direction is predictable in
advance and points at the desired result, which is the one situation where the control cannot
be deferred. Otherwise "we found enrichment" and "the method produces enrichment" are
indistinguishable afterwards.

1. **The off-host null.** Run the whole test on a rust in a host it cannot infect. Any
   enrichment is the method's floor. Cheapest and most direct.
2. **The axenic controls.** A single-pathogen library has no competitor. Quantify it with and
   without congeners in the reference; the difference is the bias, measured on ground truth.
3. **Chapter 2.** Known two-species co-infection, mocks, replicates, and the correct answer
   already computed by DESeq2 and GO enrichment.

Crossmapper would estimate cross-mapping for a chosen reference set and read length, as
Morales-Cruz et al. did by simulation, but it is not in Setonix modules either and controls 1
and 2 measure the same thing on real libraries.

### Gene sets

All computable from the CDS already on disk, and annotation is a one-off per reference reused
across every run: secreted CAZymes and the plant cell wall degrading subset (dbCAN3), candidate
effectors (SignalP then EffectorP, with ApoplastP for localisation), secondary metabolite
clusters (antiSMASH, fungal mode), peptidases (MEROPS), and PHI-base homologs, against the
PHI-base release `02_reference` already pins.

**DECIDED 2026-11-25: annotation lives in THIS module**, under `seeds/`, with the version and
licence recorded per source the way `02_reference` records its own pins. holobase is not part
of this chapter, Leon's call. The cost is accepted and stated: a second consumer would have to
import from here or repeat the work, and moving a finished annotation set into `_refdata`
later is a copy, while waiting on holobase now blocks the module.

### The prediction is made in advance, and it can fail

Biotrophs are CAZyme-depleted and necrotrophs CAZyme-enriched is a published signature, and the
trophic axis for that split already exists in the detection figure. So the test has an a priori
direction rather than being a search for whatever comes out. If obligate biotrophs show
necrotroph-like cell wall degrading enrichment, the method is wrong and we learn that.

## The host-response design

Moved here 2026-11-17 from `TODO-RESULTS.md` section 1b, with its reasoning intact and one
consequence of the wheat decision added. The question: do samples carrying more than one
pathogen differ transcriptionally from singly infected ones?

### What the index already gives it, and the one thing it does not

Host CDS is in the index for the other two questions, as a read sink so host reads are absorbed
rather than forced onto a pathogen. Quantifying it per gene rather than discarding it is the
whole of the extra work on the alignment side, so question 3 costs one more count table and no
second index, no second alignment pass and no second reference decision. That is why it is not
its own module.

The one thing it does not give is the uniqueness rule. **Do not carry species-unique-reads
across to the host, unexamined.** That rule exists to stop one fungus's reads being apportioned
onto a congener, and bread wheat is hexaploid: the A, B and D homoeologues of a gene are far
closer to each other than two congeneric fungi are, so a large and NON-RANDOM fraction of wheat
reads align equally well to two or three homoeologues of their own gene. Discarding every one
of those does not remove a bias, it removes the conserved core of the transcriptome, and it
removes more of it from exactly the gene families that retain all three copies. The Kraken2
version of this is already a recorded trap in CLAUDE.md, where 27.7% of the *Triticum* clade's
minimizers sat at the genus node.

**DECIDED 2026-11-25: gene level, homoeologues collapsed, and the collapse comes from our own
alignments.** Gene level is free, because the Ensembl CDS headers carry `gene:` and transcripts
aggregate to a locus with no extra table. The collapse is the part that has to be built, and
there were two routes. The published 1:1:1 triads (Ramirez-Gonzalez et al. 2018) are defined on
v1.1 gene IDs, cover only the 1:1:1 subset so roughly half the gene set is left unassigned, and
need a v2.1-to-v1.1 ID map before they can be used at all. **The route taken is to derive the
groups from the BAM:** put an edge between two host genes whenever a read aligns equally well
to both, and take connected components as the quantification units. That is defined by our own
read length rather than by somebody else's ortholog call, it covers every gene including the
dyads and the tandem duplicates, and it needs no external table.

**Its failure mode is transitive merging through large gene families**, where one promiscuous
read chains two unrelated loci into one component. So the edge carries a minimum weight, the
share of a gene's reads it accounts for, and the component size distribution is reported rather
than assumed. Check the components against the published triads afterwards, as a validation and
never as an input.

### The host reference

**DECIDED 2026-11-25: Chinese Spring, one cultivar**, Ensembl Plants
`triticum_aestivum_refseqv2`, GCA_018294505.1, IWGSC RefSeq v2.1. It is already on disk and
recorded in `02_reference/data/raw/host_cds_manifest.tsv` as
`Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.cds.all.fa.gz`, fetched with its md5, alongside
the v1.0 annotation of the same cultivar.

**The align index must not inherit db2's host sink, and this is the one place the two modules
deliberately disagree.** `db2` holds 19 wheat cultivars because a classifier wants every host
read absorbed and does not care which cultivar absorbed it. Quantification cares: 19
near-identical gene sets make nearly every host read ambiguous between cultivars before
homoeology is even considered, and that is a decision the uniqueness rule should never have to
make. One coherent gene model per locus is what a count table needs.

The cost is that the cohort is global cultivars and Chinese Spring is not any of them, so genes
absent from Chinese Spring are invisible and divergent alleles align with mismatches. bowtie2
tolerates the second. The first costs coverage of the dispensable genome, and **it does not bias
the contrast**, because it is the same reference on both sides of single versus multi. State it
as a limitation on which genes can be seen, not on which samples can be compared.

### The design, and the constraint that drives it

**Never compare across sites.** The big wheat projects are surveys, not replicated experiments,
and site is confounded with everything.

The wheat cohort's structure, from `release/01_cohort`:

| project | runs | what it is |
|---|---:|---|
| PRJEB39201 | 550 | the global survey. Many sites, few runs each |
| PRJEB31334 | 153 | multi-site survey |
| PRJEB15280 | 138 | multi-site survey |
| PRJEB84390 | 100 | multi-site survey |
| PRJNA718488 | 30 | INOCULATED, Zurich, one site. The manipulated contrast |

**The wheat decision cost this question its primary design, and a measurement has now replaced
the argument about it.** The design recorded on 2026-10-23 made **PRJNA746402** primary: maize,
one site in Pontevedra, one design, 59 informative runs split 14 single against 45 multi. It is
maize, so it is out, and no single wheat project replaces it. **29 blocks do**, measured on the
`db2` calls on 2026-11-25:

| | |
|---|---:|
| field *T. aestivum* runs | 1,379 |
| exactly one wheat pathogen called | 732 |
| two or more | 547 |
| none | 100 |
| runs carrying real coordinates | 1,062 |
| project x site blocks | 506 |
| blocks holding >= 2 single AND >= 2 multi | **29** |
| runs inside those blocks | **260** |
| median runs per informative block | **8** |

Computed by `research/03_align/blocks.py`, which prints this table and the grid sweep below it.

`single` is exactly one call with `is_wheat_pathogen = yes`, `multi` is two or more, and the 100
runs with none are reported but are not in the contrast: a run with no call carries a resident
community rather than nothing, which is a third state and not a zero. The resident taxa are
deliberately not in the definition.

PRJEB39201 carries 17 of the 29 blocks (149 runs), then PRJEB31334 6 (45), PRJEB36485 2 (17),
PRJEB15280 2 (12) and PRJEB84390 1 (11). **PRJNA718488 is the 29th and stays OUT of the
estimate**: it is inoculated at one site, so it is the manipulated check and answers a different
question. The block key is BioProject plus coordinate, and it is stable across grid precision:
rounding to 100 m, 1 km and 10 km gives 29, 29 and 28 informative blocks, so the inventory is
not an artefact of where the rounding falls.

**A median of 8 runs per block, rather than the 4 the four-cereal note assumed, is what makes
this a design rather than a consolation.** The analysis is a within-block contrast over 29
blocks, reporting within-block effects only and never a pooled one, which is also what makes
geography a nuisance variable here rather than a result. Only 1,062 of the 1,379 runs carry
real coordinates, so a third of the cohort cannot enter a blocked analysis at all; that is a
limitation to state, not a gap to fill with country centroids, which would invent blocks.

**Condition on severity, designed in rather than checked after.** A run carrying two pathogens
probably carries more disease, so match or covary on the dominant pathogen's coverage.
Otherwise the result reads as amount of infection rather than diversity of it, which is a
different claim and a less interesting one.

**Causality is not recoverable from this design either way** and the Discussion has to say so.
A host difference between singly and multiply infected field samples is an association. Nothing
in an archived survey establishes which came first.

## What the module can claim, and what it cannot

**Can:** "species X's expressed transcriptome is enriched for secreted plant cell wall degrading
enzymes relative to its own genomic repertoire, in N of M runs where it was detected."

**Cannot:** causation of disease, and for a necrotroph not even the distinction between
infection and colonisation of dying tissue, because a necrotroph on senescing tissue and a
saprotroph express much the same enzymes.

**The sharp exception is obligate biotrophy.** A rust or a powdery mildew cannot grow on dead
tissue, so expression of a biotroph-specific programme is evidence of a living interaction
rather than of decay. This is the strongest claim in the module and it is available precisely
where the cohort is richest.

## Measured scope, 2026-11-02 — STALE, AND KEPT FOR ITS SHAPE NOT ITS NUMBERS

**Read this section for what it establishes and not for any figure in it.** Every number below
comes from `research/02_classify/data/gold/detections.tsv`, which is db_v4 over the FOUR-CEREAL
cohort. Both halves of that have since changed: the chapter is wheat only as of 2026-11-17, and
the database is the db0/db1/db2 ladder. Maize and sorghum carried most of the secondary-detection
volume, so the event counts will fall substantially and the species list will shift toward the
wheat pathogens and the wheat resident core.

**What survives the recount is the finding that depth is not the binding constraint**, because
that is a statement about the depth distribution of secondary detections rather than about how
many there are. Recompute the table against the shipped tier on the wheat cohort before scoping
anything from it.

Depth is not the binding constraint, which was not obvious in advance. Of 3,387 classified
runs, 1,869 carry two or more fungal detections. Secondary (rank 2 and below) detections have
a median of 21,856 reads on CDS, 83.5% above 5,000 and 51.8% above 20,000. A ~12,000-gene
fungal transcriptome at 20,000 reads cannot support a per-gene test and does not need to; the
statistic aggregates a set of a few hundred genes.

Candidate events at >=20,000 reads:

| subset | events | species | runs |
|---|---|---|---|
| all secondary detections | 4,483 | 98 | 1,637 |
| pathogen = yes | 2,044 | 48 | 1,176 |
| obligate biotroph genera | 928 | 6 | 854 |
| pathogen, non-biotroph | 1,116 | 42 | 604 |

The 20,000 figure is a provisional cut, not a derived floor. It should be set the way the depth
floor was set: by rarefaction against the control arm, where the answer is known. Do not carry
it into a manuscript as if it were measured.

## The first build: job 1 only

**Build presence confirmation and stop there.** The three jobs are ordered for a reason and the
first one is also the only one that waits on nothing: it needs no gene set annotation, no
enrichment statistic and none of the host decisions above, and it is one small index rather than
one per host. It is also what the shipped module already defers to: `03_classify`'s README
names two things no threshold will settle and both are this job, and the `samples` and
`multiplicity` captions each end by handing the question here.

### Target 1: the *Alternaria* clusters

186 runs carry a flagged congener cluster and *Alternaria* is the genus in 176 of them. Across
the cohort it is 2,316 calls over 240 runs and 20 species, at a median within-run coverage CV of
0.245 against *Puccinia*'s 0.678. One run calls fourteen of them within a few percent of each
other, which is not what fourteen co-infecting organisms look like, so the uniformity is the
diagnostic and the question is whether the genus is one organism here or twenty. Nothing at the classifier's disposal can answer it, because the
genus' k-mer space is shared by construction.

The answer is the three presence facts, per species: species-unique reads, gene breadth, and
positional spread within a transcript. Plus competitive survival, the same reads quantified with
and without the congeners in the index, which is control 2 and is cheap here because the index is
small.

### Target 2: the tier-A cases, where `db2` can remove a false call but cannot make the true one

*Alternaria triticina* has **0** `db2` calls: it is a wheat pathogen that was in the reference
scope, had no gene model, and the model we predicted makes it a competitor only. *Stemphylium
lycopersici* fell from 22 `db1p` calls to 1 in the same move. So the database is now right about
the organism that is not there and still silent about the one that is.

**Alignment is where a predicted-CDS taxon can be confirmed, because the evidence is different
in kind.** Species-unique reads spread over many genes and spread evenly within them is not what
cross-mapping onto a predicted gene model produces. That is a presence claim the classifier
could not make rather than a floor being lowered. **It must be labelled wherever it is quoted:**
a positive resting on CDS we predicted ourselves is a different provenance from a positive
resting on an NCBI annotation, and the output carries the distinction as a column, not as a
sentence in a caption.

### The index, and the reads

One index, not one per host: 114 Pleosporaceae taxa are in the reference, 41 of them
NCBI-annotated and 21 with CDS predicted here, of which 71 are *Alternaria* (27 annotated, 7
predicted), plus Chinese Spring CDS as the sink. Tagged headers as above. The rest of the
family is in scope and unbuildable, which is itself part of the adjudication: a cluster that
collapses onto the one congener with a gene model is a different finding from one that spreads.

The reads are 176 runs and **917 Gbases**, mean 5.2 Gbases per run, concentrated in PRJEB39201
(96 runs), PRJEB15280 (27), PRJEB84390 (18) and PRJEB65589 (12). Every number in this section
is printed by `research/03_align/scope.py`, against the shipped tables rather than from a note.

**Throughput comes from a timed dry run on three runs of different depth, never from a guess.**
bowtie2 `-k` over 917 Gbases is this module's first real cost and nothing in the project has
measured it. Size the SLURM request from that dry run, and remember the memory request is the
CPU bill.

Order inside job 1: build the index, then the dry run and the sizing, then the 176 runs, then
the off-host null on the same index (control 1) **before** any adjudication is read. A
cross-mapping floor measured after the result is a floor fitted to the result.

## Layout, proposed

    03_align/
      01_index/     bowtie2 index over the SHIPPED tier's CDS (db2), tagged headers, Chinese
                    Spring CDS as the sink. One index per question's taxon set, not per host
      02_align/     bowtie2 -k + samtools per run; BAM is the product, not a count table
      03_assign/    read -> species-unique or ambiguous; per-gene unique counts; positional
                    spread; the three presence facts. ALSO emits the host gene-level count
                    table, which costs no extra alignment pass
      04_function/  gene set annotation, then the within-sample enrichment test
      05_host/      the single-versus-multi host contrast, site-blocked. Between-sample and
                    conventional, unlike 04_function. Reads 03_assign's host table
      controls/     the off-host null, the axenic with/without-congener pair, and the
                    Chapter 2 validation. Built alongside 04_function, never after it
      seeds/        gene set definitions and their provenance
      report/       figures and tables
      _common.py

**Job 1 is `01_index`, `02_align`, `03_assign` and `controls/` only.** `04_function`,
`05_host` and `seeds/` are the second build and are specified above but not scoped.

## Open questions before anything is built

0. **CLOSED 2026-11-25. Question 3 survives, run blocked.** 29 blocks hold both classes at a
   median of 8 runs each, which is twice what the four-cereal note assumed, so the blocked
   design is the design rather than a fallback. Within-block effects only, severity covaried,
   PRJNA718488 kept out as the manipulated check. See "The design, and the constraint that
   drives it".
1. **Which subset.** The obligate-biotroph events are where the claim is strongest and where
   alignment cost is affordable; the full secondary-detection set is the ambitious version.
   Scoping is no longer a concession now that prevalence is not the claim. **Recount on the
   wheat cohort and the shipped db tier first**: the counts in the measured-scope section above
   are four-cereal db_v4 and the biotroph set in particular was carried by rust surveys that
   are mostly wheat, so it may shrink less than the total does.
2. **CLOSED 2026-11-25. Annotation lives in this module**, under `seeds/`, with a pin and a
   licence per source. holobase is not part of this chapter.
3. **Cost, unsized.** Two separate one-offs: annotating db_v4's reference set, and aligning. The
   second is the new one, because alignment is materially slower than pseudoalignment and
   `02_classify` has no measurement for it. Size both from a dry run, never from holobase
   counts, which over-count usable references roughly threefold.
4. **The score margin** for calling a read species-unique. Needs a value and a sensitivity
   curve, not a constant.
5. **bowtie2 end-to-end or local.** Local soft-clips, which helps against divergent references
   but can let a conserved domain carry a read that end-to-end would reject. Test on the
   controls rather than choosing on principle.
6. **CLOSED 2026-11-25. Gene level, homoeologues collapsed, collapse derived from the BAM**
   rather than from the published triads, with a minimum edge weight and the component size
   distribution reported. See "What the index already gives it, and the one thing it does not".
