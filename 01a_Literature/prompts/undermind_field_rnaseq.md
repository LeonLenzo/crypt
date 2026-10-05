# Undermind prompt: published field plant RNA-seq

Purpose: size the universe of published studies that sequenced bulk RNA from field-collected
plant tissue, with no pathogen or host constraint. This is the sampling frame for a prevalence
estimate, so what is being searched for is the *denominator*, not the positives.

Run as written, record the result in `../data/`, and note the date and any wording changed.
A prompt that was edited between runs produces two corpora that cannot be compared.

---

## Prompt, v1 (not yet run)

```
I am building a comprehensive corpus of published studies that performed bulk
RNA sequencing on plant tissue collected from the field, to estimate how
commonly plant samples carry fungal, oomycete or nematode infections and
co-infections.

Find primary research papers that:
  - sequenced bulk RNA (total RNA or mRNA, not single-cell, not amplicon or
    metabarcoding, not targeted capture) from plant tissue;
  - sampled that tissue from field, farm, orchard, vineyard, forest or other
    uncontrolled outdoor settings, including surveillance and population
    genomics surveys, landscape transects, and multi-site or multi-year
    collections. Exclude studies whose plants were grown in a glasshouse,
    growth chamber, or in vitro, and exclude detached-leaf and deliberate
    inoculation experiments;
  - deposited the reads in a public archive (SRA, ENA, DDBJ, GSA).

Include studies regardless of what they were looking for. I specifically want
studies with no pathogen focus at all, for example abiotic stress, drought,
phenology, ripening, population genetics, breeding trials and gene expression
atlases, as well as pathogen surveillance. Do not restrict by host species or
by pathogen.

For each study, report: host species, country and years of collection, number
of samples sequenced, library preparation (polyA selection, rRNA depletion,
total RNA), the archive accession, and whether per-sample metadata linking
each accession to a collection site and date is available in the paper or its
supplementary material.
```

## Why each clause is there

**"bulk RNA"**, with single-cell, amplicon and capture excluded. Detection depends on
untargeted sequencing of whatever RNA was in the tissue. A metabarcoding study answers the
co-infection question directly and differently, and would be double-counting; a capture panel
cannot see an organism it was not designed for.

**The setting clause lists settings rather than saying "field"**, because "field" in a title
frequently means a field *trial* of deliberately inoculated plants. The explicit exclusion of
detached-leaf and inoculation designs is the clause doing the most work.

**"regardless of what they were looking for"**, with the no-pathogen-focus examples named.
This is the whole point of the new frame. The existing cohort was selected by querying known
PHI-base pathogens, so it is 72% single-pathogen-focus studies by construction. A drought
time-course in a paddock is worth more here than another rust survey, because it carries no
expectation about what should be present.

**The per-sample metadata question.** Asking the search to report it makes supplement quality
a property of the retrieved record rather than something discovered months later. It is also
the inclusion criterion that fixes provenance at the frame instead of patching it downstream.

**Library preparation.** 6,060 of the 6,467 currently classified samples have library prep
unstated in SRA. polyA selection versus rRNA depletion changes detection sensitivity for the
organisms being counted, so a prevalence estimate needs it, and papers normally state it in
the methods even when the submission does not.

## What to check first in the result

1. **Does it return PRJNA1217477?** That is the UC Davis *Botrytis* inoculation atlas, already
   registered in `../../_exclusions.py` as unusable. It is published, it has a supplement, and
   it has per-sample metadata, so it passes every mechanical filter and is rejected only by
   reading the methods. If the search returns it, the exclusion clause is not working and the
   prompt needs tightening before the yield means anything. PRJNA526829 and PRJNA328045 are
   the same test, less severe.
2. **The pathogen-focus balance.** If nearly every hit is a surveillance or disease study, the
   search is reproducing the old bias and the no-pathogen-focus half of the frame is missing.
   Expect this to happen to some degree: a drought study is valuable here precisely because it
   never mentions a pathogen, which is also what makes it hard to retrieve.
3. **Accession presence.** Count how many records come back with a usable archive accession
   rather than "available on request".

## Expected failure mode

Semantic retrieval will over-return surveillance and under-return the abiotic-stress and
phenology studies that are the most valuable part of the frame. That is a limit of searching
on text, not of the prompt, and it is the reason the metadata screen over the wider SRA
universe remains the complementary half of this work rather than a discarded alternative.
