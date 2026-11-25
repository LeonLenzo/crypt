# Getting to Results — crypt / phd/01-review

Written 2026-10-23 at the end of the Methods drafting session. Background and the reasoning
behind each item are in the memory note `crypt-session-state-2026-10-23`; this file is the
worklist only.

The manuscript is a .docx in OneDrive kept open in Word. Edit it through
`tools/docxedit/` (`word_review.ps1`, `word_reply.ps1`, `word_insert.ps1`,
`word_replace.ps1`), never by writing the file.

## The inclusion funnel this list assumes

    cohort as curated                7,688 runs  41 projects
    exposure = open                  6,282       30   (-1,406)
    + natural infection only         6,239       28   (-43)
    + grown in the field throughout  6,135       27   (-104)
    + above depth floor 0.463 Gkmer  3,580       26   (-2,555)

Retained set: maize 1,374, sorghum 1,118, wheat 1,032, 56 other, **rice 0**. Every rice run is
either disease-suppressed by design or below the depth floor. The rice pilot agrees: of 81
classified rice runs across four projects, 13 had any fungal call and 2 had a pathogen-labelled
one, those being an insect pathogen and a *Penicillium*; the rest were yeasts and saprotrophs.

## 1. Blocked on data

- [ ] **Classify the rest of the retained set. QUEUE REORDERED 2026-10-26 by pilot fungal
      share, not by size**, because a project's pilot predicts whether the whole thing is
      readable (`crypt-host-fraction-confound`). **2,346 runs are classified** as of
      2026-10-26 afternoon: PRJNA753546 landed 94/94 in 202 min with 0 failures and its 91 new
      runs classified in 8m53s (job 51608745), 186-379 species rows per run. Classifying while
      the fetch is still running costs nothing -- classify.py skips a non-empty report as
      `cached`, so the chained sweep just confirms them.

      Submitted 2026-10-26 as job 51599349 (fetch) with 51599350 (classify) chained on
      `afterany`, 445 runs and 0.64 TiB, in this order:

        PRJNA753546  maize  89 runs  51.4% fungal  100% co-infected
        PRJNA505891  maize  21        30.9%        100%
        PRJNA746402  maize 222        20.9%         74%
        PRJNA896566  maize  66        20.0%        100%
        PRJNA721299  other  32         9.0%        100%
        PRJNA564795  maize   8         8.8%        100%
        PRJNA369690  maize   3        15.4%        100%
        PRJNA516527  maize   2         5.7%        100%
        PRJEB31334   wheat   2        75.8%         42%

      ABANDONED, 1,231 runs and 3.55 TiB, every pilot below 1.2% fungal: the PRJEB83049 tail
      (280 runs), all of PRJEB67964 (730), PRJNA527782 (177), PRJNA383416 (38), and the
      handful left in PRJEB37824 and PRJNA1119650. Fetching them costs four days of the single
      ENA stream for runs that cannot answer the question.

      The previous fetch (51403212) hit its 24 h wall on 24 Oct and was not resubmitted: 330
      runs fetched, 24 HTTP 403 failures clustered late, per-stream rate halving from 0.78 to
      0.40 MB/s, 48 part-downloads left in reads_tmp.
- [x] **Controls: the pipeline is BUILT and the first control is classified** (2026-10-26).
      The primary control is a controlled single-pathogen infection, not an axenic culture;
      design in `01_literature/controls/DESIGN.md`, candidate set in `control_candidates.tsv`,
      order in `FETCH_PRIORITY.md`. Fetch with
      `02_classify/utilities/controls/fetch_controls.py` (plan comes from ENA's filereport at
      fetch time, reads land in `data/bronze/controls/`, NOT under `reads/`), classify with
      `classify.slurm` under `CONTROLS=1`.

      *P. triticina* / PRJNA613154: **exactly two species pass all three floors in both runs,
      the inoculated pathogen (coverage 27.6% / 31.1%) and the host (7.8% / 7.1%), out of 607
      and 639 species rows with a db_v4 denominator.** Top false positive is *H. vulgare* at
      0.08% coverage, 12x below the floor; *P. striiformis* sits ~2,800x below *P. triticina*.
      So at this depth saturation never becomes the deciding test -- everything dies on the
      absolute or coverage floor first. Read it as ONE point at the shallow end: these are the
      two thinnest runs in the control set (0.11 Gkmer, 49 nt single-end), and saturation is
      not scale-free.
- [ ] **Two gaps before the controls can answer the floor question.**
      1. `detect.py` drops any run absent from `run_depth.tsv`, so it cannot process a control
         at all. Needs a `control_depth.tsv` with `read_len` and `layout` per run; the ENA
         capture in `data/bronze/controls/meta/` already carries what it takes
         (`base_count / read_count`, and the fastq file count gives layout).
      2. Span the DEPTH axis, which is the sampling rule. Job 51609023 is queued `afterany` the
         cohort fetch: **PRJEB29046** (Ptr, 6 runs, the complete template -- mock, infected and
         axenic arms in one deposit, so it is the real test) and **PRJEB58154**
         (*Z. tritici*, 18 runs, highest undeclared load at 343).
- [ ] **Competitive read assignment (`04_align`)** on a now-defined target set: *Alternaria*
      (260 runs carrying >=3 species of the genus), and the maize *Cercospora* and *Bipolaris*
      clusters, where the real pathogen sits at ~42% coverage and congeners trail it.
- [ ] **Resolve *Stemphylium lycopersici***: 35 runs, 32 of them wheat, 2.1% median coverage, no
      sibling of the genus called in those runs, so it cannot be congener smear. It is currently
      on the pathogen map at 22 localities.

## 1b. Opened 2026-10-23 by the sorghum classify

- [ ] **Fungal share, not crop, explains the unselected-field contrast.** PRJEB83049 landed
      (264 runs, all ok) and is 0.289% fungal, median 0 fungi called, 4/284 co-infected,
      despite being the DEEPEST project in the cohort at 4.09 Gkmer. PRJEB67964's pilot is
      worse at 0.061%. Both are the same farm, UNL Havelock, Nebraska. Meanwhile unselected
      humid-site maize is 20-31% fungal with 5-12 pathogen genera and 91-100% co-infected,
      declaring nothing. Detail in the memory note `crypt-host-fraction-confound`.
- [ ] **The 0.463 Gkmer floor is on TOTAL depth and does not catch this.** Effective fungal
      depth (total Gkmer x fungal share) is bimodal with a 150x gap and nothing between 0.14
      and 1.48. It CANNOT become an inclusion criterion, because unlike total depth it is an
      outcome, so filtering on it conditions on the result.
- [ ] **Decide what replaces the pooled headline.** It is 38.8% genus-level over 1,617 runs,
      down from 46.7% over 1,333, purely because 284 blind runs entered the denominator.
- [ ] **Site-level climate join** to test the humidity reading. The curated coordinates exist;
      the maritime/continental split in `prep_pathogen_map.py` is a country list and puts
      South Carolina with Nebraska, so it cannot test it.
- [ ] Curate the missing location for PRJNA505891, the strongest single piece of evidence
      (21 runs, 100% co-infected, 30.9% fungal, median 10 pathogen genera) and currently
      without a curated site.

## 1b. 03_align, the host-response question

**MOVED 2026-11-17 to `research/03_align/README.md`, "The host-response design".** It is job 3
of that module, not a separate thing, and the design now lives beside the index it depends on.
Two changes were made there and are not reflected below: the wheat decision removed
**PRJNA746402**, which was the primary and is maize, so the question has a secondary design and
no primary; and the host-quantification rule is an open question, because hexaploid wheat makes
the species-unique rule the wrong instrument within the host.

What follows is the 2026-10-23 version, kept because the site-structure table is the evidence
the blocked design rests on. The four-cereal numbers in it need recomputing on the wheat cohort.

Leon's third question for this paper: do samples carrying more than one pathogen differ
transcriptionally from singly infected ones? Per-gene host quantification is not produced
anywhere in the current pipeline.

**Never compare across sites.** The big wheat projects are surveys, not replicated experiments:
PRJEB39201 is 538 runs over 240 sites, roughly two runs per site, and no site in any of the
multisite projects holds 10 runs in both infection classes.

What the site structure does allow:

    project       runs  sites  sites with both classes  runs in those  median runs/site
    PRJEB39201     538    240                       55            274               4
    PRJEB31334     153     75                       14             80             4.5
    PRJEB15280     132     75                       12             50               3
    PRJEB84390     100     32                        7             52               5
    PRJNA746402     61      1                        1             59              59
    PRJNA718488     30      1                        1             27              27
    PRJEB47693      23      8                        0              0               0

Design that follows from it:

- [ ] **Primary: PRJNA746402**, maize, Spain (Pontevedra, Salcedo), one site, one design, 59
      informative runs split 14 single against 45 multi.
- [ ] **Secondary: a site-blocked analysis of PRJEB39201**, using the 55 sites that contain both
      classes (274 runs) with site as a block. Per-site power is negligible at a median of four
      runs, so the strength is in the number of blocks: a stratified or mixed model reporting
      within-site effects only, never a pooled contrast.
- [ ] **Condition on severity, designed in rather than checked after.** A run carrying two
      pathogens probably carries more disease, so match or covary on the dominant pathogen's
      coverage. Otherwise the result reads as amount of infection rather than diversity of it.
- [ ] PRJNA718488 (Zurich, 27 runs, 17 against 10) is inoculated, so keep it as the manipulated
      contrast rather than as part of the estimate.
- [ ] Causality is not recoverable from this design either way, and the Discussion has to say so.

## 2. Decisions

Settled 2026-10-23:

- [x] **The funnel above is the inclusion rule**, and rice is reported as an EXCLUDED stratum
      with the reason given per project (rain exclusion, plant protection, or below the depth
      floor), rather than dropped silently. The exclusions rest on design read from each study's
      methods, and saying so is what distinguishes them from conditioning on the outcome.
- [x] **PRJNA718488 is a validation set, outside the denominator.** Excluded from the prevalence
      estimate by the natural-infection rule, but reported as the internal positive control:
      *P. triticina* 8/15 at 14.70% coverage inoculated against 3/15 at 1.98% fungicide-protected.
- [x] **Headline is all evidence tiers, with both bounds reported**: 46.7% genus-level, with
      41.3% (any PHI-base host evidence) and 39.2% (host-resolved only).

Still open:

- [ ] Story axis. **HELD OPEN deliberately until the maize and sorghum runs land**, because the
      physiology rate currently rests on 117 of 2,515 runs. Focus is a binary collapse of
      `sampling_selection`, and in the retained set it tracks crop almost exactly, so the choice
      is presentational rather than analytical.

        pathology   1,065 runs  12 projects   976 classified   46.2% co-infected
                    wheat 1,017, maize 27, other 18, barley 2, sorghum 1
        physiology  2,515 runs  15 projects   117 classified   59.0% co-infected
                    maize 1,347, sorghum 1,117, other 36, wheat 15

- [ ] Curate guild labels for the 36 gap species, 34 of which are simply absent from holobase.
      *Microdochium majus* is the test case: genus-level FUNGuild only, recorded by USDA on
      wheat, barley and maize, 90 calls currently counted as nothing.

## 3. Methods text

- [ ] Add the disease-exposure definition (open / sheltered / protected / undetermined, curated
      from each study's methods independently of its result) and the funnel with its reasons.
- [ ] Add the exclusions registry: three BioProjects excluded on design, with evidence.
- [ ] Restore reference-data provenance: a fixed snapshot of NCBI Taxonomy, GBIF Backbone, USDA
      Fungus-Host (CC0), PHI-base (CC-BY 4.0) and FUNGuild, with version and licence per source.
      Keep it separate from the assembly retrieval date.
- [ ] Fix "collected under uncontrolled field conditions"; move depth into Detection and
      classification so S1 reads continuously; fix "admitted though route two".
- [ ] Zotero entries still missing: Kraken2 (Wood et al. 2019), Lander and Waterman 1988,
      Roberts et al. 2004.

## 4. Figures

- [ ] Pathogen map: mark evidence tier, or restrict to host-resolved species. *Cladosporium
      cladosporioides* (37 localities) and *Epicoccum nigrum* (34) currently plot as surveyed
      pathogens next to *Z. tritici* (151).
- [ ] Colour `detect` panel (c) by evidence tier rather than guild. One line in `detect.R`;
      `detect_species.tsv` already carries the column.
- [ ] Regenerate all figure data and figures after the classify queue lands. Everything is
      current as of 2026-10-23, including the `prep_geography.py` paired-library fix.

## 5. Writable now

The wheat sections. 971 of 1,032 retained wheat runs are classified, declared-pathogen recovery
is 96.6% (1,061/1,098), the host call agrees with curation for 97.4% of runs, and the
depth-versus-rate correlation is r = -0.055 over 19 countries, i.e. no depth artefact. The crop
comparison is the part that has to wait.

## Deferred

- Benchmark against Zebra-style cover filters and KrakenUniq unique counts.
- Split `target` into `phibase` and `funguild` route names in `routes.py`, at the next genuine
  rebuild only, since re-running `01_refs` means fresh NCBI queries and a different assembly set.
- The gene-prediction decision.

---

## Added 2026-11-18 — the congener tier, measured

Came out of `projects/wobbegong`, a from-scratch classifier built as a learning exercise.
Two measurements that bear on `release/02_reference`, neither of which was run before.

### 1. The stated rationale for db1 is wrong for divergent genera. FIX THE SENTENCE.

`release/02_reference/README.md` currently says:

> a reference absent from the database is invisible and unrecoverable because the reads go
> silently to the nearest taxon present

The second clause does not survive measurement. On full CDS sets at k=31, congeneric
*Puccinia* share **0.10 to 0.50%** of their k-mer space (rust against Ascomycete is 87 k-mers
in 14 million, i.e. chance). Holding one species out of an index containing only its
congeners, and simulating 20,000 error-free 150 bp reads:

| held out | >=1 k-mer matched | >=5 k-mers matched |
|---|--:|--:|
| *P. striiformis* | 3.28% | 1.97% |
| *P. triticina* | 6.94% | 4.25% |
| *P. graminis* | 5.97% | 3.69% |

**About 96% of reads from an absent congener match nothing at all.** They are lost, not
misassigned. The first clause ("invisible and unrecoverable") stands; the mechanism in the
second does not. Error-free reads make these upper bounds.

Worse, the interception that does occur is pathological: median 5 of 120 possible k-mers
among reads that matched. A few reads on a few conserved k-mers is the *Phakopsora* signature
(248,504 reads on 1,239 minimizers). Congeners are not rescuing absent species here, they are
manufacturing the artefact the complexity criterion exists to reject.

**Caveat, and it is the important one.** This is *Puccinia*, a divergent genus (~16 ANI points
between species). *Alternaria* sect. Infectoriae is the counterexample already in our own
results: 24 species sharing most of their k-mers, collapsing to the genus node, 2,985 of
wheat's 5,506 above-criterion detections. **Generalising from rusts to all genera would repeat
the error that produced the Alternaria denominator bias.** Do not act on this beyond fixing
the sentence until the per-genus measurement below is run.

Rerun: `projects/wobbegong/experiments/congener_holdout.py`.

### 2. The per-genus blur index — cheap, not yet run

`kraken2-inspect` places a minimizer at the genus node exactly when two or more species in
that genus share it. So for every genus in db_v4:

    blur = genus-node minimizers / (genus-node minimizers + sum of species-node minimizers)

Near zero is *Puccinia*: congeners separate cleanly, db1 buys recall only. Near one is
*Alternaria*: congeners blur, species calls are unreliable, db1 costs resolution. One pass
over `outputs/db_inspect_<name>.tsv`, no new compute, and it turns "is db1 required" into a
per-genus answer defensible in a supplementary table.

### 3. What db1 is actually made of

`research/db1_composition.py`. The scope rule reads USDA at genus, which lumps two different
justifications into one tier. Of 1,744 db1 assemblies across 322 genera:

| | assemblies | species | genera | share |
|---|--:|--:|--:|--:|
| A covers an unresolvable USDA record | 1,032 | 750 | 167 | 59.2% |
| B speculative congener | 556 | 390 | 93 | 31.9% |
| ? genus absent from the USDA screen | 156 | 94 | 62 | 8.9% |

A means USDA named something in that genus that NCBI cannot place, so the real organism may be
any species in it and the congeners are covering a real but unresolvable observation. B means
every USDA record in that genus resolved cleanly and is already represented in db0, so the
additions are congeners never observed on a plant.

**The majority is defensible.** The 8.9% unattributable bucket needs explaining before review:
those genera are not in the screen at all, so they arrived by exception or by a spelling the
screen does not share.

**Case B is dominated by sequencing effort, not biology.** *Coemansia* contributes 154
assemblies on 3 USDA records, nearly 9% of the whole tier, for a Kickxellomycotina genus from
dung, soil and insect guts. *Wallemia* 52 on 1 record. *Trichophyton* 37 on 1 and *Microsporum*
7, both animal dermatophytes. *Hanseniaspora* 17 and *Clavispora* 9, yeasts. This is the same
shape as the FUNGuild problem in `crypt-pathogen-label-not-inherited`: a genus-level rule
inheriting whatever NCBI happens to have sequenced heavily.

Decide before submission whether case B is reported as a deliberate inclusion or trimmed. It
is already built, so this is a defensibility question, not a compute one.

---

## Updated 2026-11-18 (evening) — blur index on the FULL release reference

The `db1p` build finished (7h21m, kraken2-inspect committed
`release/02_reference/data/outputs/db_inspect_db1p.tsv`, 10.09B root minimizers, 99.26% Fungi).
Reran `db1_blur.py` and `db1_blur_impact.py` on it. The silver-subset numbers held to within
hundredths of a percent, confirming blur is a property of the taxa, not the db version:

  Puccinia 0.45%, Alternaria 32.66%, Fusarium 27.70%, Aspergillus 18.27%, Rhynchosporium 62.05%.

**Impact (figures/db1_blur_impact_db1p.png): Alternaria remains the SOLE genus both common on
wheat and blurred** (2,378 calls, 32.7%). Fusarium 19 calls, Aureobasidium 1, Akanthomyces 1.
Caveat: the call counts are still from the silver classification; re-classifying against db1p
is a separate run. Blur values are final, call counts provisional.

**New finding worth acting on: Coemansia is a double-hit.** It is the single largest SPECULATIVE
congener block (154 assemblies on 3 USDA records, db1_composition.py) AND the most-blurred genus
in the whole reference with real sequence (49.67%, 146 species, 137M minimizers). Weakly
justified and resolution-costing. Strongest single case in the reference for trimming a db1
genus. Morchella (41%), Pseudogymnoascus (36%), Claviceps (32%) are lesser versions.

Next: (a) the Alternaria -> genus rollup in reporting, blur column in the supplement; (b) decide
Coemansia/Morchella/etc trim; (c) rerun db1_blur_impact once the cohort is classified against db1p.

---

## Added 2026-11-19 — minimizer-crossover supplementary figure

Lives in `release/02_reference/figures/crossover/` (moved out of research/, follows the
matrix/tree figure convention: plot_crossover.py + caption.py + crossover_matrix.tsv, with
compute_matrix.py/.slurm as the Setonix provenance). `caption.py` generates
crossover_caption.txt from the data.

Pairwise minimizer sharing among the 254 wheat-row taxa (90 genera, 1,289 assemblies),
l=31/k=35 minimizers from db1p CDS, exact set intersection. Result: crossover is block-diagonal
by genus (mean within-genus 0.015, cross-genus 0.0001). 12 of 25 wheat pathogens are hard to
resolve (crossover >=20% with a congener): Tilletia bunt pair, Fusarium head-blight group,
Alternaria/Pyricularia/Ustilago complexes. The 3 wheat rusts (Puccinia) resolve cleanly and are
NOT flagged, the key separability result. Largest off-diagonal: Alternaria alternata x
tenuissima, 8.83M minimizers.
