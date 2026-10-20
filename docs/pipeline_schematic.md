# Pipeline schematic

What each step consumes, what it produces, and what depends on it. Regenerated **2026-10-20**
against the two live modules, after the four earlier modules were retired to `superseded/`.

**This document deliberately does not repeat each script's options or output columns.** Those
live in the script's own docstring, where they cannot drift out of step with the code. An
earlier version of this file duplicated them and went stale across two restructures. What
belongs here is the chain: the order, the hand-offs, and where the data physically sits.

Run every script from the repository root (`phd/01-review`). A numbered directory is a pipeline
step; an unnumbered one (`seeds/`, `slurm/`, `report/`, `utilities/`, `data/`, `logs/`) is
support.

---

## The chain

```
01_literature         the papers are the sampling frame
  retrieve/           Undermind + NCBI -> candidate studies
  studies/            the papers themselves, one directory each (manual drop zone)
  curate/             every value traced to a quoted sentence
  seeds/              AUTHORED: curation.tsv, provenance.tsv  <-- irreplaceable
  data/{bronze,silver,gold}
  release/            cohort_samples.tsv, cohort_bioprojects.tsv, evidence/*.md
        │
        │  7,688 runs / 6,222 BioSamples / 41 BioProjects
        ▼
02_classify
  seeds/              AUTHORED: who belongs in the database, and hosts
  01_refs             routes -> survey -> select -> compare        4,271 assemblies
        │             predict/  gene prediction for unannotated species (benchmarked)
        ▼
  02_build            fetch_cds -> build_db                        db_v4, 43 GB
        │                                                          db_inspect.tsv
        │                                                          ^ every coverage denominator
        ▼
  03_cohort           resolve -> fetch_depth -> fetch_reads        reads/<BioProject>/
        │             migrate_reads (one-shot, done)  freshen (keep-alive)
        ▼
  04_classify         kraken2 --confidence 0.15 --minimum-hit-groups 3
        │                     --report-minimizer-data
        │                                              reports/<BioProject>/<Run>.txt
        ▼
  05_detect           absolute -> coverage -> saturation -> f.sp.
        │                                              detections.tsv, coinfection.tsv
        ▼
  04_align            SPECIFIED, NOT BUILT. Competitive read assignment, to resolve
                      congener clusters and low-coverage/high-saturation calls
```

`01_literature` feeds `02_classify` through **five** files, all named in
`02_classify/_layout.py` so a layout change breaks in one place:

| | |
|---|---|
| `release/cohort_samples.tsv` | the cohort, one row per BioSample |
| `release/cohort_bioprojects.tsv` | per-project summary |
| `data/silver/runs.tsv` | run-level detail, including `declared_pathogens` and `BioSample` |
| `seeds/curation.tsv` | what each study declared, as curated names |
| `seeds/provenance.tsv` | host, tissue, location, year, each traced to a quote |

**`cohort_samples.tsv` does NOT carry `declared_pathogens` or `BioSample`.** Both are in
`runs.tsv`. Joining on the wrong one makes every detection look undeclared, with no error.

---

## Step by step

### 01_literature — build the cohort from the papers

The frame. It exists because the original frame was a detector: screen SRA for pathogen
signal, then read the papers. That cannot support a prevalence claim, because the gate selects
for the thing being measured and its pass rate is host-dependent (wheat 3.42%, soybean 0.55%).

| script | reads | writes |
|---|---|---|
| `retrieve/*.py` | Undermind reports, NCBI, ENA portal | `data/bronze/`, candidate lists |
| `curate/scaffold.py` | a study directory | an `adapter.yaml` stub |
| `curate/apply_curation.py` | `seeds/curation.tsv` | `data/silver/` |
| `curate/import_provenance.py` | `superseded/02_literature` provenance mining | `seeds/provenance.tsv` |
| `report/*.py` | the silver tables | `data/gold/`, `report/figures/` |
| `release/build_release.py` | everything above | `release/` |

- **`seeds/` is authored and irreplaceable.** It sits outside `data/` deliberately: `data/`
  means "delete it and it rebuilds", and these are the only files for which that is false.
  Losing `seeds/` costs re-reading every paper.
- **`release/` is the reviewer deliverable** and is generated, never hand-edited: flat TSV plus
  one Markdown evidence file per project.
- Read `curate/CURATING.md` before curating anything.

### 02_classify — the database, the reads, and the calls

Each step has its own README with the decisions that are not visible in the code. In outline:

| step | reads | writes |
|---|---|---|
| `01_refs/routes.py` | `seeds/`, holobase, USDA Fungus-Host | membership, by route |
| `01_refs/survey.py` | NCBI datasets | `silver/assemblies.tsv`, `assembly_survey.tsv` |
| `01_refs/select_refs.py` | `assemblies.tsv` | `silver/ref_candidates.tsv` |
| `01_refs/compare_refs.py` | `ref_candidates.tsv`, db_v3's | the diff |
| `02_build/fetch_cds.py` | `ref_candidates.tsv` | `bronze/cds/` (Setonix, ~200 GB) |
| `02_build/build_db.py` | `bronze/cds/` | `silver/db/db_v4/`, `db_inspect.tsv` |
| `03_cohort/resolve.py` | the five `01_literature` files | `silver/run_list.tsv`, `fetch_plan.tsv` |
| `03_cohort/fetch_depth.py` | ENA portal, one request per BioProject | `silver/run_depth.tsv` |
| `03_cohort/fetch_reads.py` | `run_depth.tsv` | `bronze/reads/<BioProject>/` |
| `04_classify/classify.py` | reads + db_v4 | `bronze/reports/<BioProject>/` |
| `05_detect/detect.py` | reports + `db_inspect.tsv` + curation | `silver/_classified.tsv`, `gold/` |

- **`db_inspect.tsv` is the hinge.** Every coverage figure divides by a number from it, which is
  why **coverage is not comparable across database versions**.
- **`select_refs.py` queries nothing.** It filters `assemblies.tsv`, which is the survey's own
  record, so the selection cannot disagree with the survey that justified it.
- **`resolve.py` cross-checks the two `01_literature` tables against each other** and fails if
  they disagree about which runs exist, rather than letting one win.
- **`fetch_mode` in `fetch_plan.tsv` is not a denominator.** It records intent. What was
  actually classified is counted from the reports.

---

## Figures, and which step owns them

All under `02_classify/report/figures/`, each a `prep_*.py` that writes a TSV plus an `.R` that
draws it, plus a **generated caption** carrying a FULL and a SHORT version.

| figure | what it shows |
|---|---|
| `db_contents` | the db_v4 ring: 1,547 species by route and class |
| `predict_counts` | prediction against annotation, three panels |
| `predict_gain` | what prediction would add, by route |
| `rarefaction` | depth against detection; read length costs nothing beyond Gkmer |
| `saturation_cloud` | the two-dimensional presence surface |
| `detect` | the detection calls, per species |
| `detect_runs` | the same, per run |

- **Generate captions from data.** A caption hand-edited to match a figure is a caption that
  will describe a figure that no longer exists; a silent `str.replace` renders the code stale.
- The PDF is tracked and the PNG is not: the PDF is what goes into Illustrator and regenerates
  the PNG at any resolution.
- Ask whether labels should be stripped for assembly before drawing. Several of these are
  printed and annotated separately.

---

## Shared root

| | |
|---|---|
| `_paths.py` | artefacts that cross module boundaries, named once. Holds the `SUPERSEDED` prefix |
| `_util.py` | domain-free helpers: `atomic_write`, `_Tee`, `make_log_dir`, `link_latest` |
| `_exclusions.py` | unusable BioProjects, with evidence |
| `manifest.py` | a CLI over one `_util` helper: what exists remotely, without the data |

Each live module also has its own `_layout.py` (paths) and shared modules beside it:
`02_classify/_reads.py`, `_reports.py`, `_assembly.py`, `_taxonomy.py`, `_seeds.py`.

**Named `_layout`, not `_paths`.** The repo root owns `_paths.py`, and a second module of that
name shadows it on `sys.path`. That collision broke `review.py` and destroyed a table once.

Nothing imports by module name: every cross-module import inserts the repository root on
`sys.path`. A leading digit in a directory name would otherwise break `import`.

---

## Where the data actually lives

| | |
|---|---|
| repo | code, small tracked TSVs, figures (PDF), log history. **No research data** |
| `02_classify/data/silver/*.tsv` | TRACKED from 2026-10-19. Tables, not layers: the twelve are 6.6 MB while `silver/db/` alone is 139 GB |
| Setonix `/scratch` | `bronze/`: CDS, reads, reports, sketches, genomes. **21-day purge on atime** |
| Setonix cohort reads | `02_classify/data/bronze/reads/<BioProject>/`, 4.1 TiB over 28 projects |
| Acacia | the CDS tree, project bucket. ~1 TB total, so reads stay on Lustre |
| `_refdata` | the shared reference store, outside the repo. Search it before downloading |
| `superseded/` | the four retired modules, on disk and untracked. See `SUPERSEDED.md` |

- **The purge clock is atime**, and `/scratch` is Lustre mounted without `noatime`, so reading
  4 KB refreshes a file as well as reading all of it. `02_classify/slurm/freshen_reads.slurm`
  is the keep-alive, and it covers the cohort tree only.
- The flat pre-migration tree holds ~9.1 TB of non-cohort STAT-gated runs. Not freshened,
  meant to expire.
- **Never compute on a login node.** `sbatch` anything that does work.
