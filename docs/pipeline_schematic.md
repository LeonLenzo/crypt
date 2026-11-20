# Pipeline schematic

What each step consumes, what it produces, and what depends on it. Regenerated **2026-11-17**
over both trees, after a memory and methods audit found this file still describing the
2026-10-20 single-tree layout and a `02_classify/04_align` path that never existed.

**This document deliberately does not repeat each script's options or output columns.** Those
live in the script's own docstring, where they cannot drift out of step with the code. An
earlier version duplicated them and went stale across two restructures. What belongs here is
the chain: the order, the hand-offs, and where the data physically sits.

**Nor does it carry status.** Whether a module is built, running or waiting is in the module's
own README and in the `crypt-session-state-<date>` memory note. A status line in this file is a
status line that will be wrong.

Run every script from the repository root (`phd/01-review`).

---

## Two trees

```
release/     the reviewer deliverable. Public sources only, rebuilt forward from its own
             inputs, nothing copied in. This is the tree that gets a Zenodo snapshot.
research/    the working tree: the history, the reasoning, and the design docs written
             before a thing is built. Also where anything touching a private database lives.
```

A `release/` module is the boiled-down version of the question its `research/` counterpart
answered, and supersedes it when it lands. They are SIBLINGS; nothing is copied from one to
the other. See `crypt-release-tree` in memory for the four tiers inside a release module and
why `inputs/` differs from `raw/`.

The numbering does not line up across the trees, and that is not a mistake:

| release | research | why they differ |
|---|---|---|
| `01_cohort` | `01_literature` | renamed: the cohort comes from databases as well as papers |
| `02_reference` | `02_classify/01_refs` | reference scope split out as its own module |
| `03_classify` | `02_classify/02_build`..`05_detect` | build, classify and detect together |
| `04_align` | `03_align` | research has no separate reference module, so it is one lower |

---

## The release chain

```
01_cohort            the cohort: who is in it, and what is known about each run
  inputs/dois.txt          503 DOIs. AUTHORED, the one irreplaceable input
  fetch_title_abstract ->  raw_papers.tsv
  host_screen          ->  curation/host_screen_papers.tsv
  fetch_fulltext       ->  raw_fulltext.tsv + raw/fulltext/*.xml
  bioproject_screen    ->  curation/bioproject_screen_papers.tsv
  curation/declared_pathogens.tsv  what the STUDY declared, per BioProject, with its quote
  fetch_samples / screen_samples / audit
  make_supp_tables     ->  outputs/supp_data_s1_studies.tsv
                           outputs/supp_data_s2_cohort.tsv   <-- the cohort, one row per run
  figures/                 cohort_map, cohort_network
        │
        │  wheat. Related cereal hosts in cohort BioProjects are the host contrast;
        │  lab libraries are the controls. Counts in the module README.
        ▼
02_reference         which organisms the database contains, and which assemblies supply CDS
  inputs/scope_rule.tsv        the criterion, with a pin on every dataset it is applied to
  inputs/scope_exceptions.tsv  admitted outside the master, and what was rejected
  inputs/wheat_pathogens_curated.tsv  the 8 wheat diseases PHI-base does not curate
  fetch_species        ->  raw_usda_fungus_host.csv, raw_taxdmp.zip
  fetch_pathogens      ->  raw_phibase.csv
  screen_species       ->  curation/screen_species.tsv      scope + on_wheat
  screen_merges        ->  curation/screen_merges.tsv       two names, one fungus
  screen_pathogens     ->  curation/screen_pathogens.tsv    the pathogen label
  fetch_assemblies     ->  raw_assembly_summary.txt
  screen_assemblies    ->  curation/screen_assemblies.tsv
  make_supp_tables     ->  outputs/supp_data_s3_reference.tsv  <-- ONE ROW PER TAXON
  plan_build           ->  curation/build_plan.tsv             <-- the db tier per assembly
  plan_extra_cds       ->  curation/extra_cds_db2.tsv          <-- the 135 predicted genomes
  fetch_cds --tier     ->  CDS on Setonix
  figures/                 matrix, tree
        │
        │  the db ladder. kraken2-build --add-to-library is incremental, so each
        │  tier reuses the last and the subtraction between them is clean:
        │
        │     db0  baseline     route in {usda_species, exception}, annotated
        │     db1  +congeners   db0 + route usda_genus, annotated
        │     db2  +predicted   db1 + the 135 gene-predicted genomes (research/predict)
        │
        │  db1-db0 measures competition among already-represented genera;
        │  db2-db1 measures reach into genera nothing could represent before.
        │  PREDICTED CDS IS A COMPETITOR ONLY, never a positive call.
        ▼
03_classify          Kraken2 against the ladder, then the detection calls
                     Supersedes research/02_classify/{02_build,03_cohort,04_classify,05_detect}
                     The THREE FLOORS and the saturation statistic carry across unchanged;
                     see research/02_classify/05_detect/README.md for their derivation
  inputs/thresholds.tsv            the operating point, authored, floor by floor
  inputs/null_biotroph_orders.tsv  which orders the off-host null is built from
  classify.py --arm     ->  data/raw/reports/<arm>/<PRJ>/<Run>.report.tsv   on Setonix
                            data/raw/reports_manifest_<arm>.tsv
  calibrate_floors.py   ->  curation/null_offhost_<arm>.tsv   where the null is exhausted
  detect.py --arm       ->  outputs/detections_<arm>.tsv      the calls
                            data/raw/candidates_<arm>.tsv     every candidate, on scratch
  make_supp_tables.py   ->  outputs/supp_data_s5_detections.tsv  one row per call
                            outputs/supp_data_s6_runs.tsv        one row per cohort run
  compare_arms.py       ->  outputs/arm_resolution.tsv        db2 - db1p, by assigned rank
  figures/resolution/                                         the purpose test

                     TWO ARMS, db1p and db2. db0 is not one: the copy on disk has no host
                     sink, so db1p-db0 would measure the sink. Coverage divides by each
                     arm's own inspect table, so the floor is re-anchored per arm and a
                     coverage figure is never restated across arms by arithmetic.
        │
        │  this is the SCREEN: which organisms are present and worth characterising.
        │  It is explicitly NOT a prevalence. See CLAUDE.md.
        ▼
04_align             NOT BUILT. Competitive assignment, then function
                     Design and reasoning: research/03_align/README.md
                     bowtie2 -k against a combined per-host tagged CDS index.
                     Species-unique reads only; NO EM.
```

---

## The research tree

`research/` is not a dead archive. Three parts of it are load-bearing right now.

| | |
|---|---|
| `01_literature/` | the curation audit chain. Every curated value traces to a quoted sentence, and `seeds/curation.tsv` + `seeds/provenance.tsv` are AUTHORED and irreplaceable |
| `02_classify/` | the db_v4 results the release tree's scope decisions were measured against, and the derivations (three floors, saturation, the rarefaction, the off-host null) |
| `predict/` | the gene-prediction arm, moved out of `01_refs` on 2026-11-09. vipsania over 135 genomes, DONE and validated. Supplies db2 |
| `03_align/` | the align DESIGN doc. No code |
| `superseded/` | four retired modules, on disk and untracked. See `SUPERSEDED.md` |

`research/02_classify`'s chain, for reading its outputs:

| step | reads | writes |
|---|---|---|
| `01_refs/routes,survey,select,compare` | seeds, USDA, NCBI datasets | `silver/ref_candidates.tsv` |
| `02_build/fetch_cds,build_db` | `ref_candidates.tsv` | `silver/db/db_v4/`, `db_inspect.tsv` |
| `03_cohort/resolve,fetch_depth,fetch_reads` | the `01_literature` tables | `bronze/reads/<BioProject>/` |
| `04_classify/classify.py` | reads + db_v4 | `bronze/reports/<BioProject>/` |
| `05_detect/detect.py` | reports + `db_inspect.tsv` + curation | `silver/_classified.tsv`, `gold/` |

- **`db_inspect.tsv` is the hinge.** Every coverage figure divides by a number from it, which
  is why **coverage is not comparable across database versions**: db_v3 against db_v4, and
  equally db_v4 against db0/db1/db2.
- **`select_refs.py` queries nothing.** It filters the survey's own record, so the selection
  cannot disagree with the survey that justified it.
- **`fetch_mode` is not a denominator.** It records intent. What was classified is counted from
  the reports.

---

## Shared root

| | |
|---|---|
| `research/_paths.py` | artefacts that cross module boundaries. Holds the `SUPERSEDED` prefix |
| `research/_util.py` | domain-free helpers: `atomic_write`, `_Tee`, `make_log_dir`, `link_latest` |
| `research/_exclusions.py` | unusable BioProjects, with evidence |
| `tools/dataflow.py` | which script reads and writes which file |
| `tools/docxedit/` | the manuscript is a .docx in OneDrive open in Word; edit it through here |

Each research module has its own `_layout.py` (paths) plus shared modules beside it.

**Named `_layout`, not `_paths`.** The root owns `_paths.py`, and a second module of that name
shadows it on `sys.path`. That collision broke `review.py` and destroyed a table once.

Nothing imports by module name: every cross-module import inserts the repository root on
`sys.path`. A leading digit in a directory name would otherwise break `import`.

Release modules do not use `_layout.py`. Each script resolves its own module root from
`__file__` and reads `data/` beneath it, because a release module is meant to be readable and
runnable on its own.

---

## Where the data actually lives

| | |
|---|---|
| repo | code, small tracked TSVs, figures (PDF), log history. **No research data** |
| `release/*/data/{inputs,curation,outputs}` | TRACKED. These are the deliverable |
| `release/*/data/raw/` | fetched bytes, md5-pinned where it matters. Large bodies gitignored, manifests tracked |
| Setonix `/scratch` | `bronze/`: CDS, reads, reports, sketches, genomes. **21-day purge on atime** |
| Acacia | the CDS tree, project bucket. ~1 TB, so reads stay on Lustre |
| `_refdata` | the shared reference store, outside the repo. Search it before downloading |
| `superseded/` | the four retired modules, on disk and untracked |

- **The purge clock is atime**, and `/scratch` is Lustre mounted without `noatime`, so reading
  4 KB refreshes a file as well as reading all of it. The keep-alive covers the cohort tree
  only.
- **Predicted CDS lives on Setonix scratch only** and is not mirrored to Acacia.
- **Never compute on a login node.** `sbatch` anything that does work.
