# `disease_exposure` — the curation rule

**Written 2026-10-21, BEFORE the 41 projects were read, deliberately.** The levels below were
fixed in advance so they could not be fitted to the detection rates, which were already known
for 25 projects when this was written. That ordering problem is real and is disclosed in
Methods rather than hidden: the covariate was *conceived* after observing that six rice and
sorghum projects returned zero, so only its definition, not its motivation, is blind.

## What it records

**What the source experiment did to its plants that bears on whether foliar fungal infection
could establish and be sampled.** It is a property of the study design, readable from the
methods section, and it is assigned WITHOUT reference to what the classifier found.

It is NOT a judgement about whether a pathogen was present, and it is NOT a filter. Nothing is
removed from the cohort on the basis of this field. It exists to stratify.

## Why it is needed

The cohort is unbiased with respect to **pathogen signal**, because no step in its construction
conditions on whether a pathogen was present or declared. It is NOT unbiased with respect to
**disease exposure**, because it inherits whatever the underlying experiments did, and those
experiments were designed to answer other questions. Large public cereal RNA-seq is dominated
by abiotic-stress physiology and by breeding trials, and those two kinds of study treat disease
very differently.

## The levels, in order of increasing suppression

Assign the STRONGEST level that applies. A study that both excludes rain and applies fungicide
is `protected`, not `sheltered`.

| level | meaning | typical evidence |
|---|---|---|
| `open` | grown outdoors with no documented measure that reduces disease exposure | a field trial that mentions fertiliser and nothing else |
| `sheltered` | rain, dew or free water deliberately excluded or interrupted for any part of the growing period | rain-out shelter; closing canopy; "protecting the field from rain"; drought imposed by rain exclusion rather than by withholding irrigation alone |
| `protected` | a chemical or biological treatment applied that suppresses fungi | fungicide, seed treatment with a fungicidal active, "general plant protection" |
| `indoor` | the sequenced plants were grown entirely under cover | glasshouse, growth chamber, polytunnel, controlled environment |

**Structural rain exclusion goes in `setting`, not here.** Leon's call, 2026-10-21. `setting`
already qualifies rather than splits: `field (chamber-raised, transferred outdoors)` and
`field (inoculated in situ)` are existing values. A rain-out facility is a fact about where the
plant grew, so it belongs there as **`field (rain-excluded)`**, which keeps the project in the
cohort while recording exactly what differed. Qualifying also avoids reclassifying a study
against its own description: the Baihe authors call theirs a field experiment.

Use the `setting` qualifier only where exclusion covers the WHOLE experiment. Where it covers
one arm and the arms cannot be separated by a predicate, record `disease_exposure: sheltered`
at project level and say so in the note. PRJNA510920 is that case: rain exclusion applies to
the W-deplete half of an open IRRI field, and the SampleName codes `4s`/`4w`/`9s`/`9w` are
never defined in the paper.
| `unknown` | the methods do not say, or full text is unavailable | |

## A separate, ORTHOGONAL field: `tissue_selection`

Sampling protocol determines whether disease, if present, would be *collected*. It is recorded
separately because it is a different mechanism from suppression, and a study can do either,
both or neither.

| level | meaning | typical evidence |
|---|---|---|
| `representative` | a healthy or typical plant chosen; edge plants avoided; a named upper leaf taken | "a single representative plant... avoiding edge plants... the fourth leaf below the flag leaf" |
| `positional` | tissue taken by a fixed position or stage with no health criterion | "flag leaf at anthesis"; "days after anthesis" |
| `symptomatic` | tissue chosen because it showed symptoms | lesion sampling, disease surveys |
| `unspecified` | the methods do not say | |

## Rules of assignment

1. **Only the sequenced plants count.** Papers frequently describe several experiments. If the
   RNA-seq samples came from the field arm, a glasshouse arm elsewhere in the paper does not
   make the project `indoor`. Record which arm the quote describes in the note.
2. **Absence of evidence is `unknown`, never `open`.** A paper that says nothing about disease
   management is not evidence that none was applied. `open` requires a methods section detailed
   enough that a measure would have been mentioned.
3. **Quote the sentence.** Every assignment gets an `evidence_id` into `provenance.tsv` with
   the locus and the quote, like every other curated field.
4. **A negative scan is evidence and should be recorded as one.** "fungicide 0, pesticide 0,
   inoculation 0, greenhouse 0 across the full text" supports `open` where the methods are
   otherwise detailed. PRJNA505891 already carries exactly such a note and it is the model.
5. **Do not look at the detection rate while assigning.** If you already know it, say so in the
   note.

## Known hazard

At the time of writing, every project known to carry a suppressing feature returned a 0%
detection rate, and no high-rate project carried one. That is either the mechanism working as
plant pathology predicts (leaf wetness drives foliar infection) or an artefact of having seen
only six such projects. **It means the field and the outcome are currently confounded**, so the
stratified result must be reported with the project count per stratum, and a stratum carried by
one or two projects must not be presented as a rate.
