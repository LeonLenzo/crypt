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

**Scope: aerial tissue only.** Root and rhizosphere material carries a large soil-derived
community, so a fungal detection there is as likely to be a saprotroph or root endophyte as a
pathogen, and the signal-to-noise ratio is much worse than in leaf or spike tissue. This also
keeps the new frame consistent with the existing one, where the Kraken2 target is field and
aerial and `llm_tissue` splits the classified cohort 5,339 aerial to 862 non-aerial.

## What the frame has to carry

The hypothesis is that polymicrobial infection is more common in the field than is reported,
and that public RNA libraries hold the evidence. **Cryptic is defined against the authors'
own claim**: in a study naming pathogen X it means anything detected beyond X, and in an
abiotic or no-stress study it means any pathogen at all, since nothing was expected.
Deliberate multi-pathogen designs are explicitly not cryptic.

That definition puts a requirement on the frame that a pure prevalence frame does not have.
Each study must come with **what its authors claimed about pathogens, including claiming
nothing**. Without that, a detection is just a detection and cannot be called cryptic. It is
the single most important thing this module collects, and it is why the search asks for the
authors' stated pathogens alongside the accessions.

It also produces a comparator the database-first route cannot. Asking every retrieved study
what it reported gives the rate at which field studies report co-infection **themselves**,
which is the "than we think" half of the hypothesis and has never been measured here. The
existing cohort cannot supply it, because the STAT gate only ever admitted studies that
already had a detection.

### A deposit with no manuscript is rejected

**Inclusion criterion, Leon's call 2026-10-06: an archive deposit with no available manuscript
does not enter the cohort, however good its metadata.** Recorded in `seeds/decisions.tsv` at
`level = policy`.

This follows from the paragraph above rather than adding to it. If cryptic is defined against
what the authors claimed, then a deposit nobody has written up has no claim to define it
against. There is no methods section to verify the setting, no statement of which pathogen was
expected, and so a detection in it cannot be classified as cryptic or not. Such a deposit is
not weak evidence, it is evidence of a different kind, and mixing it in would quietly change
what the denominator means.

Two consequences worth stating plainly, because this criterion costs something real:

- It rejects exactly the material the parked metadata field screen is designed to surface
  (`crypt-metadata-field-screen`): large unpublished surveillance deposits, often with better
  per-sample metadata than published studies. That route therefore cannot feed this cohort. It
  remains the right tool for measuring **what this cohort missed**, which is a separate job.
- "Has a manuscript" is a claim about the DATA, not about a single paper. The Sato 2024
  Arabidopsis cluster qualifies: its RNA-seq has no dedicated paper, but the field experiment
  is described in Sato et al. 2024 and the sequencing in Tomita et al. v3, so both the setting
  and the authors' silence on pathogens are documented. A deposit fails this test when nothing
  anywhere describes it.

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

1. **Search.** Semantic literature retrieval for bulk RNA-seq of field-collected aerial plant
   tissue, with no host or pathogen constraint. The aim is completeness: every qualifying
   study counts once, surveillance and drought alike, because a denominator is only a
   denominator if nothing was chosen for it. See `prompts/undermind_field_rnaseq.md` for the
   prompt and the reasoning behind each clause.
2. **Resolve paper to accession.** The inverse of what `02_literature/01_search` does, and it
   fails differently: papers deposit without a per-run mapping, give a BioProject with no
   sample-level key, or say "available on request". This attrition sets the real yield and
   should be measured on the first 20 papers before anything is committed to.
3. **Apply inclusion criteria by hand.** Field-collected, aerial tissue, bulk RNA, public
   reads, per-sample metadata linking each accession to a site and a date, and a recorded
   statement of what the authors said about pathogens. Record the reason for each
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

Two top-level ideas: **seeds/ is authored, data/ is derived.** `curate/CURATING.md` is the
runbook; read it before curating a project.

```
01a_Literature/
  _layout.py _hostgroup.py _cohort.py _resolver.py   shared, importable, not entry points
  retrieve/   triage sources accession_papers resolve_accessions   external fetching
  resolvers/  ada21 kashima sato24 epicon range28    one per awkward study
  curate/     apply_curation review scaffold import_provenance + CURATING.md
  report/     cohort host_summary offtarget_screen project_brief truncation_check
  seeds/      AUTHORED and irreplaceable, tracked: curation joins provenance decisions
              found run_scope backfill_accessions, and undermind/ (the search returns)
  data/bronze/  immutable external captures: runinfo, biosample_attrs, geo, ena (ignored)
  data/silver/  generated tables: runs bioprojects papers paper_bioproject ... (ignored)
  data/gold/    small tracked outputs: cohort host_summary offtarget, the gap analyses
  studies/    per-study drop zone: papers, supplements, resolver join tables (ignored)
  prompts/    search prompts, versioned, with the reasoning for each clause
```

The rule: if losing a file costs a refetch it belongs under `data/`; if it costs re-reading
papers it belongs in `seeds/`.
