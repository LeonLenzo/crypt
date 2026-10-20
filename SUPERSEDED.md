# superseded/ — the retired modules

Four modules were retired on 2026-10-20 and moved to `superseded/`. They are **on disk and out
of the repository**: 709 files untracked, 74 MB of tracked content removed, 60 MB of which was
figures. The git history keeps every version of all of it.

This file is tracked and `superseded/` is not, so this is the only place in the repo that
records what each module did and what replaced it. Read it before reaching into that directory.

## Why they were retired rather than deleted

Each one is an evidence chain. `01_literature`'s curated values trace to quoted sentences that
`02_literature` mined, and db_v4's reference set was chosen partly by comparison with db_v3's.
Deleting the modules would break those traces. Retiring them ends a different problem: six
module directories in the working tree where two are live, four of which describe frames the
project has moved off.

## The four

### `01_stat` — NCBI STAT k-mer screening

Screened 608,368 SRA runs on NCBI's own k-mer profiles and retained 10,995. It was the original
entry point: find runs with pathogen signal, then read the papers.

**Replaced by `01_literature`**, which inverts that. The cohort is built from the literature
first and carries no pathogen-signal condition, because a signal-gated denominator cannot
support a prevalence claim. STAT's pass rate is also host-dependent (wheat 3.42%, soybean
0.55%) and its blind spot is measured, not suspected: 0% eukaryotic signal on 15 runs where
Kraken2 and kallisto both find 65 to 68% *Puccinia striiformis*.

Still read by live code:

- `03_filter/data/runs.tsv`, the STAT-retained set, by `01_literature/retrieve/triage.py`,
  `report/stat_queue.py` and `report/truncation_check.py`. It is used to ask which cohort runs
  STAT would have found, which is a question about the old frame, not a gate on the new one.

`02_fetch/data/stat_cache.jsonl` is 2.3 GB and was never tracked. It stays on disk.

### `02_literature` — NCBI search and full-text LLM classification

Resolved BioProjects to papers, gated on retrievable full text, and classified 6,467 BioSamples
by host, setting and declared pathogen. `03_classify/data/samples.tsv` was the sample universe
the whole pipeline ran on.

**Replaced by `01_literature`.** The full-text gate was the module's defining constraint and
also its ceiling: only projects with retrievable text could be classified at all.
`01_literature` curates from the paper outward instead, with every value traced to a quoted
sentence, and reached host and tissue 100%, location 99.8%, year 97.2% on a 7,688-run cohort.

Still read by live code:

- `03_classify/data/cohort_provenance.tsv` and `02_text/data/supp_provenance.tsv`, by
  `01_literature/curate/import_provenance.py` and `retrieve/geocode_localities.py`. This is the
  provenance-mining result that raised location coverage before curation plateaued, and it
  enters the curation chain with `locus=02_literature provenance mining`.
- `02_text/data/text_cache.jsonl` and `02_text/data/manual_pdfs/`, via
  `01_literature/_layout.py`. The PDFs are publishers' copyrighted files and were never
  tracked.

### `03_kraken` — db_v3, classification and the detection filter

Built db_v3 (2,115 assemblies, 1,043 taxa), classified 3,220 of 3,225 runs, called hosts from
reads, and produced the detection filter and co-occurrence networks.

**Replaced by `02_classify`.** Not a rename and not a port: every input changed, and two
dependencies on the superseded frames were live gates rather than stale paths.
`kraken_select.py` resolved BioSample to Run through `01_stat`'s STAT-retained table, so only
1,227 of 7,688 cohort runs could be selected at all; and `05_filter/detect.py` skipped any
taxon outside PHI-base's 632 fungal taxa, which under db_v4 would mean every resident.

Still read by live code:

- `01_search/data/ref_candidates_v3.tsv`, by `02_classify/01_refs/compare_refs.py`, as the
  db_v3 baseline a database change is measured against.
- `05_filter/data/_classified.tsv`, by `02_classify/report/figures/prep_rarefaction.py`, as the
  db_v3 arm of the rarefaction comparison.

**Two traps in that table.** Its `score` and `cls` columns are the superseded statistic: a
fitted accumulation curve, with `cls == "real"` meaning score >= 0.5 and no coverage floor,
which is the setting that gave wheat 80% co-infected. `05_detect` writes its own
`_classified.tsv` with `call` and `reason` against saturation instead. And **coverage is not
comparable between the two**, because the denominator is each build's own inspect table.

`03_select/freshen.py` was the last piece of this module still executing. It hardcoded the flat
pre-migration reads directory, which `migrate_reads.py` had drained, so it kept reporting a
healthy tree while freshening nothing the cohort used. Replaced by
`02_classify/03_cohort/freshen.py`, which takes its root from `_layout.READS`.

### `04_kallisto` — competitive read assignment

Six strata, one kallisto index per host, intended to confirm Kraken2 calls by competitive EM
and to calibrate the detection thresholds.

**Superseded by the planned `02_classify/04_align`**, which answers the same question against
db_v4's reference set. Its measurements are worth carrying across:
`02_build/data/neighbours.tsv` holds the containment between reference species, which is how
cross-talk was quantified without going through taxonomy.

## Rules

**Do not re-run anything in `superseded/` without fixing its paths first.** Those scripts use
repo-root-relative strings such as `Path("03_kraken/05_filter/data/three.tsv")`, which now
resolve one directory too high. They were left as they were rather than rewritten, because
rewriting a retired script invites running it.

**Live code reaches these outputs through `_paths.py` only.** Every retired constant there is
re-rooted on one `SUPERSEDED` prefix, so a second move is a one-line edit. Do not add a new
hardcoded `superseded/...` string anywhere.

**One historical string was deliberately left stale.**
`01_literature/curate/import_provenance.py:194` records `02_literature/02_text/supp_provenance.py`
inside an audit `locus=` message. That is a statement about where a curated value came from at
the time it was curated, so it is correct as written and must not be updated.

## Reads on scratch

The flat tree `superseded/03_kraken/03_select/data/reads` is a symlink target on Setonix holding
about 9.1 TB: the 2,041 non-cohort STAT-gated runs. They were held open for the db_v3-versus-db_v4
comparison, which is closed. **They are no longer freshened and are meant to expire** under the
21-day purge. Every one is a public SRA accession, re-pullable from `run_list.tsv`.

The live cohort tree is `02_classify/data/bronze/reads`, one directory per BioProject, 4.1 TiB
over 28 projects, kept alive by `02_classify/slurm/freshen_reads.slurm`.
