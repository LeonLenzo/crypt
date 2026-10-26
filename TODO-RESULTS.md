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
