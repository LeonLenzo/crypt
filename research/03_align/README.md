# 03_align: which genes the reads come from, and whether they are pathogenicity genes

Sketched 2026-11-02. Nothing here is built. This file is the design and its reasoning; it
supersedes the `02_classify/04_align` named in `SUPERSEDED.md` and CLAUDE.md, which could never
have had that path because `02_classify/04_` is `04_classify`. Align is its own module for the
reason CLAUDE.md already gives: a different method (competitive EM, not k-mer counts) answering
a different question.

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

## The two questions, in order

`02_classify` establishes that an organism's DNA is present above three floors. It cannot say
which genes the reads came from, which makes two different things impossible:

1. **Is this organism actually here, or are its reads a congener's?** The old `04_kallisto`
   module was built to answer exactly this and its design is carried forward below.
2. **Is this organism expressing a pathogenic lifestyle?** New, and the reason the module is
   being revived now rather than when it was first specified. Leon 2026-11-02: we accept that
   causation of disease is not available from archived field samples, but "this fungus is
   expressing genes associated with pathogenicity" is a weaker claim the data can carry.

The order is not negotiable. A functional test run on a cross-mapping artefact characterises
the function of an organism that is not there. **Confirm presence, then characterise.**

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

What does NOT carry forward is the reference set. It must be db_v4's, or coverage and
containment are measured against different denominators and nothing is comparable.

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
clusters (antiSMASH, fungal mode), peptidases (MEROPS), and PHI-base homologs, which holobase
already holds.

**Open: where annotation belongs.** It is a property of the reference set, not of a run, and
`_refdata`/holobase is the single source of reference truth. But its only consumer today is this
module, and the repo's convention is that a thing lives with its consumer. Leon's call. Putting
it in holobase is the better long-run answer and the more expensive one now.

### The prediction is made in advance, and it can fail

Biotrophs are CAZyme-depleted and necrotrophs CAZyme-enriched is a published signature, and the
trophic axis for that split already exists in the detection figure. So the test has an a priori
direction rather than being a search for whatever comes out. If obligate biotrophs show
necrotroph-like cell wall degrading enrichment, the method is wrong and we learn that.

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

## Measured scope, 2026-11-02, from `02_classify/data/gold/detections.tsv`

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

## Layout, proposed

    03_align/
      01_index/     bowtie2 index per host over db_v4 CDS, tagged headers, host CDS as sink
      02_align/     bowtie2 -k + samtools per run; BAM is the product, not a count table
      03_assign/    read -> species-unique or ambiguous; per-gene unique counts; positional
                    spread; the three presence facts
      04_function/  gene set annotation, then the within-sample enrichment test
      controls/     the off-host null, the axenic with/without-congener pair, and the
                    Chapter 2 validation. Built alongside 04_function, never after it
      seeds/        gene set definitions and their provenance
      report/       figures and tables
      _common.py

## Open questions before anything is built

1. **Which subset.** The 928 obligate-biotroph events (854 runs, 6 species, >=20,000 reads) are
   where the claim is strongest and where alignment cost is affordable. The full 4,483-event /
   98-species set is the ambitious version. Scoping is no longer a concession now that
   prevalence is not the claim.
2. **Where annotation lives.** It is a property of the reference set and `_refdata`/holobase is
   the single source of reference truth, but its only consumer today is this module and the
   repo's convention is that a thing lives with its consumer. Holobase is the better long-run
   answer and the more expensive one now. Leon's call.
3. **Cost, unsized.** Two separate one-offs: annotating db_v4's reference set, and aligning. The
   second is the new one, because alignment is materially slower than pseudoalignment and
   `02_classify` has no measurement for it. Size both from a dry run, never from holobase
   counts, which over-count usable references roughly threefold.
4. **The score margin** for calling a read species-unique. Needs a value and a sensitivity
   curve, not a constant.
5. **bowtie2 end-to-end or local.** Local soft-clips, which helps against divergent references
   but can let a conserved domain carry a read that end-to-end would reject. Test on the
   controls rather than choosing on principle.
