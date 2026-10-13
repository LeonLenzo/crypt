# Undermind prompt: published field plant RNA-seq

Purpose: size the universe of published studies that sequenced bulk RNA from field-collected
aerial plant tissue, with no pathogen or host constraint.

The hypothesis is that polymicrobial infection is more common in the field than is reported,
and that public RNA libraries carry the evidence. "Cryptic" is defined against what the
authors said: a detection is cryptic when the study never reported that organism. So the
frame needs two things from each study, not one. It needs the samples, and it needs **what
the authors claimed about pathogens**, including claiming nothing. Without the second, a
detection is just a detection and the word cryptic cannot be applied to it.

That is also why this search returns a baseline the database-first route cannot: the rate at
which field studies report co-infection *themselves*. That is the "than we think" half of the
hypothesis, and it has never been measured here.

Run as written, record the result in `../data/`, and note the date and any wording changed.
A prompt that was edited between runs produces two corpora that cannot be compared.

---

## Prompt, v1 (not yet run)

```
I am building a comprehensive corpus of published studies that sequenced
bulk RNA from aerial plant tissue collected in the field. The aim is to
establish how often such samples contain pathogen signal that the original
authors did not report, so I need the studies themselves, not studies about
co-infection.

Find primary research papers that:
  - sequenced bulk RNA (total RNA or mRNA, not single-cell, not amplicon or
    metabarcoding, not targeted capture) from plant tissue;
  - sampled aerial tissue: leaf, stem, spike, head, panicle, flower, fruit,
    seed, bark or whole shoot. Exclude root, rhizome, tuber, nodule and
    rhizosphere samples, and exclude whole-seedling or whole-plant samples
    that include roots;
  - sampled that tissue from field, farm, orchard, vineyard, forest or other
    uncontrolled outdoor settings, including surveillance and population
    genomics surveys, landscape transects, and multi-site or multi-year
    collections. Exclude plants grown in a glasshouse, growth chamber or in
    vitro, exclude detached-leaf assays, and exclude any study that
    deliberately inoculated its plants, including co-inoculation and
    mixed-infection experiments;
  - deposited the reads in a public archive (SRA, ENA, DDBJ, GSA).

Be comprehensive. I want every study meeting those criteria, not a selection
of the best examples. Include studies regardless of what they were looking
for: pathogen surveillance, abiotic stress, drought, phenology, ripening,
population genetics, breeding trials, gene expression atlases, and anything
else. Do not restrict by host species or by pathogen, and do not prioritise
studies that mention disease.

For each study, report:
  - host species, and which tissue was sequenced;
  - country, region and years of collection;
  - number of samples sequenced;
  - library preparation (polyA selection, rRNA depletion, total RNA);
  - the archive accession;
  - whether per-sample metadata linking each accession to a collection site
    and date is given in the paper or its supplementary material;
  - which pathogens or diseases, if any, the authors state they were
    studying or report having detected, and whether the paper describes more
    than one pathogen in any single sample. If the study is not about
    pathogens at all, say so explicitly.
```

## Why each clause is there

**"bulk RNA"**, with single-cell, amplicon and capture excluded. Detection depends on
untargeted sequencing of whatever RNA was in the tissue. A metabarcoding study answers the
co-infection question directly and differently, and would be double-counting; a capture panel
cannot see an organism it was not designed for.

**The setting clause lists settings rather than saying "field"**, because "field" in a title
frequently means a field *trial* of deliberately inoculated plants. The explicit exclusion of
detached-leaf and inoculation designs is the clause doing the most work.

**"Be comprehensive" and "regardless of what they were looking for".** This is the whole
point of the new frame, and the emphasis is on completeness rather than on any one kind of
study. The existing cohort was selected by querying known PHI-base pathogens, so it is 72%
single-pathogen-focus studies by construction. The fix is not to swap that for a preference
for non-pathogen studies, which would be the same mistake inverted. Every qualifying study
counts once, surveillance and drought alike, because a denominator is only a denominator if
nothing was chosen for it.

The no-pathogen-focus examples are named not because those studies are worth more, but
because a semantic search will under-return them by default: they do not use the vocabulary
being searched on. Naming them is a correction to the retrieval, not a ranking.

**Aerial tissue only.** Root and rhizosphere material carries a large soil-derived
community, so a fungal detection there is as likely to be a saprotroph or a root endophyte as
a pathogen, and the signal-to-noise ratio for the question being asked is much worse than in
leaf or spike tissue. This also keeps the new frame consistent with the existing one: the
Kraken2 target set is field *and aerial*, and `llm_tissue` already splits the classified
cohort 5,339 aerial to 862 non-aerial. Whole-seedling and whole-plant samples are excluded
for the same reason, since they bring roots in through the back door.

**"Which pathogens, if any, the authors state they were studying."** This is the clause the
hypothesis actually turns on, and the first draft of this prompt did not have it. Cryptic is
defined relative to the authors' own claim: in a study naming pathogen X it means anything
beyond X, and in an abiotic or no-stress study it means *any* pathogen at all, since nothing
was expected. Both readings need the claim on record. Asking for it at retrieval also yields
the author-reported co-infection rate for free, which is the comparator the whole "more
common than we think" framing needs and which the database-first route cannot produce,
because it only ever saw studies that already had a detection.

**Co-inoculation experiments are excluded by name**, not just by the general inoculation
clause. A deliberate multi-pathogen design is explicitly not cryptic: a second organism there
is the experiment. Naming them separately matters because such papers use co-infection
vocabulary heavily and are exactly what a semantic search for this topic will surface first.

**The per-sample metadata question.** Asking the search to report it makes supplement quality
a property of the retrieved record rather than something discovered months later. It is also
the inclusion criterion that fixes provenance at the frame instead of patching it downstream.

**Library preparation.** 6,060 of the 6,467 currently classified samples have library prep
unstated in SRA. polyA selection versus rRNA depletion changes detection sensitivity for the
organisms being counted, so a prevalence estimate needs it, and papers normally state it in
the methods even when the submission does not.

## What to check first in the result

1. **Does it return co-infection papers instead of field studies?** The search is for studies
   that *could* contain unreported co-infection, not studies *about* co-infection. If the top
   of the return is review articles and deliberate mixed-infection experiments, the prompt is
   being read as a topic search and needs rewording before the yield means anything.
2. **Does it return PRJNA1217477?** That is the UC Davis *Botrytis* inoculation atlas, already
   registered in `../../_exclusions.py` as unusable. It is published, it has a supplement, and
   it has per-sample metadata, so it passes every mechanical filter and is rejected only by
   reading the methods. If the search returns it, the exclusion clause is not working and the
   prompt needs tightening before the yield means anything. PRJNA526829 and PRJNA328045 are
   the same test, less severe.
3. **The pathogen-focus balance.** If nearly every hit is a surveillance or disease study, the
   search is reproducing the old bias and the frame is incomplete. Expect this to some degree,
   since the studies that never mention a pathogen are the ones a text search finds hardest to
   reach. Judge the return on coverage, not on how interesting any individual study looks.
4. **Tissue.** Spot-check that root and whole-seedling studies were actually excluded rather
   than returned and left for us to filter.
5. **Accession presence.** Count how many records come back with a usable archive accession
   rather than "available on request".

## Expected failure mode

Semantic retrieval will over-return surveillance and under-return studies that never mention
a pathogen, because the latter do not use the vocabulary being searched on. That skew matters
for completeness, not because those studies are individually worth more: a frame that is
systematically missing one category cannot give an honest prevalence estimate in either
direction. It is a limit of searching on text rather than of the prompt, and it is the reason
the metadata screen over the wider SRA universe remains the complementary half of this work
rather than a discarded alternative.
