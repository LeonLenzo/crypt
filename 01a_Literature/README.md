# 01a_Literature — build the cohort from the papers, not from the database

## Rationale

The pipeline's existing entry point is `01_stat`: query NCBI STAT for runs carrying signal
from a known PHI-base pathogen, then find the papers behind the runs that pass. That order
makes the cohort a function of what the detector found, and it biases the result in a way no
amount of downstream filtering can undo.

**The gate selects for pathogen presence, so the denominator is not a denominator.** Across
the analysed cohort, "any pathogen detected" is 93 to 100% by construction. A co-infection
rate computed on it answers "given that we already found something, how often do we find a
second thing", which is not the question.

**The gate pass rate is host-dependent.** Measured over the full 543,639 HAL runs in
`stat_cache.jsonl`: wheat 3.42%, maize 1.07%, Arabidopsis 0.78%, soybean 0.55%. Wheat
dominates the cohort because its libraries clear the gate roughly six times more readily than
soybean's, not because more wheat is sequenced. Arabidopsis has six times more screened runs
than wheat. Every cross-host comparison inherits this.

**The field cohort is survey-dominated, and the surveys are the usable part.** Of 3,282
field-classified BioSamples in 323 BioProjects, the top 10 projects hold 45% and projects of
20 or more samples hold 63%. The largest are the rust field pathogenomics corpora (PRJEB39201
at 539 samples across 30 countries, PRJEB31334, PRJEB15280, PRJNA256347, PRJEB65589,
PRJEB36485). These are real multi-site, multi-year surveys, and they are also the studies that
published per-sample supplements, which is why supplement mining lifted location coverage to
86%. Sitting among them are an inoculation atlas, an in vitro study and a controlled
inoculation, all classified `field`, all now in `../_exclusions.py`.

This module inverts the order. Find the studies first, by what they did rather than by what
was detected in them, then pull their accessions.

## What this costs, and what it does not

It is tempting to read literature-first as trading coverage for cleanliness. The trade is
smaller than it looks.

**The analysed cohort is already literature-gated.** Only BioProjects with retrievable full
text are classified at all: 732 of the 1,287 that survived the STAT gate. Every one of the
6,467 analysed BioSamples comes from a published study. So this module does not introduce
publication bias into the analysis. That bias is already present. What changes is that the
*pathogen* gate comes out of the frame.

**The frame is "published field plant RNA-seq", not the true universe.** That should be
stated plainly in methods rather than claimed away. Unpublished surveillance deposits are
genuinely field and genuinely missing from a literature search. Sizing that gap is the job of
the complementary metadata screen over the wider SRA universe (roughly 11,300 non-Arabidopsis
BioProjects), which keys on "is this field" and never on pathogen signal. The two are a pair,
not alternatives: the literature route gives a well-characterised core, the metadata route
measures what the core missed.

**The unreliable component leaves the critical path.** Today `llm_study_setting` is
load-bearing for the headline, and it is wrong often enough to matter: it calls a UC Davis
inoculation atlas `field`. Under this frame, "is this field bulk RNA-seq" becomes an inclusion
criterion applied by hand to a few hundred candidate papers, not a per-BioProject LLM judgment
over a corpus that is mostly controlled work. The failure mode changes from invisible false
positives to a methods section someone reads.

## Design

Nothing here is built. The intended shape:

1. **Search.** Semantic literature retrieval for bulk RNA-seq of field-collected plant tissue,
   with no host or pathogen constraint. See `prompts/undermind_field_rnaseq.md` for the prompt
   and the reasoning behind each clause.
2. **Resolve paper to accession.** The inverse of what `02_literature/01_search` does, and it
   fails differently: papers deposit without a per-run mapping, give a BioProject with no
   sample-level key, or say "available on request". This attrition sets the real yield and
   should be measured on the first 20 papers before anything is committed to.
3. **Apply inclusion criteria by hand.** Field-collected, bulk RNA, public reads, and
   per-sample metadata linking each accession to a site and a date. Record the reason for each
   rejection, so the frame is reproducible rather than a judgment that happened once.
4. **Hand the accession list downstream.** The detection stack does not change. `03_kraken`
   and `04_kallisto` are the instrument, and the instrument is fine; only the frame was wrong.
   `stat_cache.jsonl` also stays useful as a free per-host fungal-signal layer over whatever
   cohort this produces.

## The test set

Any frame worth adopting rejects these three on its own terms, without being told:

| BioProject | n | What it actually is |
|---|---|---|
| PRJNA1217477 | 211 | UC Davis *Botrytis cinerea* inoculation atlas |
| PRJNA526829 | 60 | *in vitro* on artificial surfaces |
| PRJNA328045 | 39 | controlled compatible/incompatible inoculation |

PRJNA1217477 is the hard case and the one to watch. It is published, it has a supplement, and
it has per-sample metadata, so it passes every mechanical filter in the design above. What
rejects it is reading the methods and finding deliberate inoculation. If the search returns
it, the prompt's exclusion clause is not working and the yield does not mean anything yet.

Full reasons and evidence in [`../_exclusions.py`](../_exclusions.py).

## Relationship to 02_literature

`02_literature` **enriches a cohort that already exists**: it takes the BioProjects the STAT
gate produced and adds metadata, full text, LLM study-design classification and provenance.
`01a_Literature` **constructs a cohort**. They are different operations on the same kind of
material, which is why this is a separate module rather than another step inside that one.

The two share tooling in the obvious direction. Supplement mining
(`02_literature/02_text/supp_provenance.py`) and geocoding
(`02_literature/03_classify/geocode_localities.py`) apply unchanged to whatever cohort this
module produces, and the per-sample metadata requirement in step 3 is the inclusion criterion
version of the provenance problem those scripts were written to patch.

## Status

Not started. The only artefact is the prompt, which has not been run. Next step is to run it
and size the result, then measure paper-to-accession yield on a sample of 20 before committing
to the approach.

## Layout

```
01a_Literature/
  prompts/   search prompts, versioned, with the reasoning for each clause
  data/      search returns and the resolved accession lists (gitignored)
```
