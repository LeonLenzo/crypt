# crypt — Cryptic Co-infection Mining from Public Plant RNA-seq

*Leon Lenzo, Curtin University (leon.lenzo@curtin.edu.au)*

## What this is

Plant disease studies deposited in the NCBI Sequence Read Archive are designed around a single
target pathogen. Field-collected samples routinely harbour additional co-infecting organisms
whose signal is present in the sequencing data but goes unexamined and unreported under
single-target designs.

This project identifies those organisms in archived wheat RNA-seq, and asks whether they are
transcriptionally active rather than merely present. It does not ask how often they occur: the
cohort cannot support a prevalence, so **the rate is the instrument and not the claim**. The
classification stage is a screen, saying which organisms are present and worth characterising;
the alignment stage carries the finding.

## Scope

| | |
|---|---|
| **Host** | wheat, *Triticum aestivum* |
| **Also in the cohort** | related cereal hosts deposited inside cohort BioProjects, as a host contrast; lab libraries, as the axenic and inoculated controls |
| **Organisms** | fungi |
| **Reference scope** | every fungus USDA Fungus-Host records on any plant, read at genus, with an annotated assembly |
| **Material** | public RNA-seq runs. Nothing was sequenced for this chapter |

Bacteria are out because polyA+ selection depletes bacterial mRNA, so their percentages do not
measure what was in the tissue. Viruses are out because their genomes are too small and
fragmented for the k-mer pass. Oomycetes are out.

**Scope is any plant, never wheat.** The study exists to find unreported host associations, and
a reference absent from the database is invisible and unrecoverable, because its reads go
silently to the nearest taxon that is present. "Recorded on wheat" and "known pathogen" are
therefore columns on the reference table, never filters on it. A detection reads off those
columns: a species with neither is a resident, one labelled a wheat pathogen is expected, and
one labelled a pathogen but never recorded on wheat is a candidate unreported interaction.

## Method

**1. The cohort is built from the papers.** It starts from a list of DOIs and works outward
through titles, full text and BioSample records, with every curated value traced to a quoted
sentence. A cohort built instead from a detector cannot support any claim about what is
present, because the gate selects for the thing being measured and its pass rate varies by
host.

**2. The reference is built for competition.** Kraken2's LCA can only decline to resolve a read
when two or more database taxa compete for its k-mers. A narrow database gives most reads
exactly one plausible match and a confident species call whether or not the species is there,
so breadth is what makes a call mean anything. The scope rule reaches the resident community
and not only the pathogens, because reads from an unrepresented organism are not left
unclassified.

The database is an incremental ladder, so the effect of each addition is measurable on its own:
`db0` baseline, `db1` adds congeners, `db2` adds 135 gene-predicted genomes for genera with no
annotation at all. Predicted CDS is a competitor only, never a positive call.

**3. Detection is three floors, in order.** An absolute floor on distinct k-mers, then a
coverage floor relative to the taxon's own share of the database, then a shape test asking
whether distinct k-mers accumulate as a real organism's would or whether reads are piled on the
same fragments. The order matters: a relative floor alone is meaningless on a small reference,
and a shape test only applies where there are enough reads to judge spread.

**4. Alignment resolves what k-mers cannot.** Counting k-mers cannot say which genes a
detection's reads came from, which leaves open whether an organism is really present or its
reads are a congener's, and whether it is expressing a pathogenic lifestyle. Both are answered
by aligning against a combined per-host CDS index with tagged transcript headers, then counting
species-unique reads only, applying the same rule to the gene set and to the background. The
same index gives the host response to multiple infection.

## What is established

These measure the method rather than the world, which is why they stand without a denominator
argument.

- **The off-host null is near zero.** Running the detection filter on organisms in hosts they
  cannot infect gives an artefact rate of 0.032% median and 1.97% maximum, with none of 6,285
  rows called. This is what gives the coverage floor a basis: 1% sits roughly 30x the null's
  75th percentile.
- **A control library returns the right answer.** On an inoculated *Puccinia triticina*
  library, exactly two species pass all three floors in both runs, the inoculated pathogen at
  27.6% and 31.1% coverage and the host at 7.8% and 7.1%, out of 607 and 639 species rows. The
  top false positive sits 12x below the floor.
- **Gene prediction is as good as submitter annotation.** vipsania over 135 genomes in genera
  with none: BUSCO 135/135 scored, median completeness 99.4% against 98.4% for the downloaded
  annotations, and zero prediction failures. The four targets below 90% are assemblies that
  were poor before prediction touched them.
- **Depth is what limits recovery, not read length.** From the rarefaction, read length costs
  nothing beyond its contribution to total k-mers.

## What is not resolved

**DNA presence is not proof of active co-infection, and coverage cannot separate the two.**
Spore drift, surface inoculum and mixed samples all deposit real organism DNA. The clearest
case is residual *Phakopsora* on wheat: read-level alignment confirms it is genuinely
*Phakopsora*-like sequence rather than misclassified stripe rust, 94% of 248,504 reads mapping
to *Phakopsora* and 18 to *P. striiformis*, yet *P. pachyrhizi* cannot infect wheat. This is
inherent to archive-based inference and is disclosed rather than thresholded away.

**Congener density manufactures confident-looking calls.** One run calls fourteen *Alternaria*
species, every one at coverage 5.6 to 7.3% and saturation 0.070 to 0.080. Fourteen genuinely
co-infecting congeners would differ wildly in abundance; one organism's reads spread across
fourteen near-identical references looks exactly like this, and the uniformity is the
diagnostic. No per-detection threshold catches it, because each call is individually
well-supported.

Both are what the alignment stage exists to settle.

## Layout

```
release/     the reviewer deliverable. Public sources only, rebuilt forward from its own
             inputs. 01_cohort  02_reference  03_classify  04_align
research/    the working tree: history, reasoning, and designs written before they are built
tools/       dataflow.py, and docxedit/ for the manuscript
```

Each module has its own README with its rationale, its methods and its current state.

| | |
|---|---|
| `docs/pipeline_schematic.md` | the chain: what feeds what, across both trees |
| `CLAUDE.md` | decisions, traps and constraints. No status, deliberately |
| `SUPERSEDED.md` | the four retired modules, what replaced each, and which of their outputs live code still reads |

## References

- **PHI-base:** Urban et al. (2020) *Nucleic Acids Res* — [phi-base.org](https://phi-base.org)
- **USDA Fungus-Host Distributions:** Farr & Rossman, Systematic Mycology and Microbiology Laboratory
- **Kraken2:** Wood et al. (2019) *Genome Biol* — [PMC6883579](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6883579/)
- **BUSCO:** Manni et al. (2021) *Mol Biol Evol*
- **bowtie2:** Langmead & Salzberg (2012) *Nat Methods*
