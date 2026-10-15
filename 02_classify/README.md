# 02_classify — DESIGNED, NOT BUILT

This README is the design, written before the code. Nothing here exists yet except
`seeds/`, which is already built and moves across unchanged.

It replaces `03_kraken`. Not a rename and not a port: **every input to that module has
changed**, and two of its dependencies on the superseded frames turned out to be live gates
rather than stale paths.

## Why a rewrite rather than a rewire

I argued against this twice, from line count: 03_kraken is ~9,850 lines and rewriting looks
like waste. That argument conflated **method** with **wiring**, and the split is what matters:

| | lines | |
|---|---:|---|
| `01_search` | 2,167 | wiring. PHI-base seeded; superseded by `seeds/` |
| `03_select` | 697 | wiring. **This is where the STAT gate lives** |
| `04_assign` | 992 | wiring. Reads samples.tsv, STAT runs.tsv, phibase_db.json |
| `05_filter` | 840 | **method.** Carried across, not rewritten |
| `utilities` | 3,394 | **method/QC.** Carried across; no STAT dependency at all |
| `02_build` | 365 | essentially the kraken2-build call; survives either way |

The rotten half is the wiring, and every one of its inputs is now different: the cohort
(7,688 literature-first runs, not 3,225 STAT-gated), the reference database (db_v4 from
`seeds/`), the declared set (curated, not LLM), the pathogen definition (holobase, not
PHI-base) and the scope (fungi only, six hosts).

### The two gates, because they are the reason

Neither was a stale path. Both silently constrained results.

**STAT gated what gets classified.** `03_select/kraken_select.py` picks target BioSamples from
`02_literature`'s LLM cohort, then resolves BioSample → Run through
`01_stat/03_filter/data/runs.tsv` — the STAT-*retained* set, 10,995 of 608,368 runs screened.
**Only 1,227 of the 7,688 literature cohort runs (16%) appear in it; 6,461 (84%) do not.** A
selected BioSample whose runs are not STAT-retained gets no runs, with no warning.
`04_assign/kraken_report_analysis.load_run_to_biosample()` has the same dependency on the
output side. This is why the literature-first cohort "has never been pushed through kraken": with
this wiring it structurally cannot be. And the STAT gate selects for pathogen presence at a
host-dependent rate (wheat 3.42%, soybean 0.55%), which is the exact bias `01_literature` exists
to escape.

**PHI-base gated what gets scored.** `05_filter/detect.py:41`:

    if str(tid) in hosts or str(tid) not in pathogens or ...: continue

`pathogens` is PHI-base's 632 fungal taxa. db_v3 holds 877 fungal species and db_v4 will hold
~1,100, **including every resident taxon the 2026-10-13 scope work added**. None of
*Cladosporium*, *Aureobasidium*, *Papiliotrema* or the phyllosphere yeasts is a PHI-base
pathogen, so `detect.py` would silently skip all of them. That is not a stale path; it is a
definition of "pathogen" narrower than the database and narrower than the question.

## What this module does

Three jobs, in order. Its output is one sentence per sample: *these organisms are present, at
this confidence, and this many of them were never reported by the study that generated it.*

## Layout

```
02_classify/
  README.md     this file
  seeds/        MOVED from 03_kraken/seeds, unchanged. The membership authority: what is in
                and out of the CDS database, authored vs generated, and THE SCOPE RULE
                (plant-associated, ANY host; never scope by the cohort's own hosts).
  01_refs/      NEW.      select assemblies from seeds/ and download their CDS
  02_build/     moved.    kraken2-build over the CDS, + the inspect table
  03_cohort/    REWRITTEN. resolve the cohort to a run list and fetch reads. NO STAT.
  04_classify/  REWRITTEN. run kraken2; call the host from reads
  05_detect/    moved.    the detection filter: real signal vs classifier artefact
  utilities/    moved.    busco, saturation, annotate (miniprot), split
  slurm/        adapted.  the Setonix job chain
```

Step names say the job rather than the tool, like the module name. Renamed from
search/build/select/assign/filter on the same principle (`04_kallisto` becomes `03_align`).

## The input contract, stated once

Everything comes from three places and nothing else:

    01_literature/release/cohort_samples.tsv    Run, BioProject, crop, host, setting,
                                                sampling_selection, tissue, location, date
    01_literature/data/silver/runs.tsv          Run -> BioSample for all 33,958 runs. The
                                                ungated replacement for the STAT run table:
                                                it covers 7,688/7,688 cohort runs over 6,222
                                                BioSamples
    01_literature/.../declared_pathogens        what the STUDY declared. COMPLETE as of
                                                2026-10-15: 41/41 projects, 7,688/7,688 runs
    02_classify/seeds/                          what the database contains
    _refdata/holobase.db                        pathogen / host / guild LABELS
    Setonix                                     CDS tree + reads, keyed by ACCESSION

**No data moves on Setonix.** Both the CDS tree (`cds/pathogen/GCA_000002495.2/…`) and the
reads are keyed by accession, so the new module points at the same directories. That removes the
main objection to a rewrite. The old databases stay (and are mirrored on Acacia).

Every path goes through `_paths.py`. No 03_kraken script used it — all seven hardcoded their
paths, which is exactly the debt `_paths.py` was written for, and exactly how a gate hid in a
string literal.

## The four design changes, each with its reason

**1. No STAT, anywhere.** Run ↔ BioSample comes from `01_literature`. The cohort comes from
`01_literature`. There is no run table that can silently drop 84% of the cohort.

**2. Pathogen status is a LABEL applied after detection, never a filter applied before it.**
From holobase: `taxon_guild` (1,835 Plant Pathogen taxa from FUNGuild), `host_association`
(`association_type='pathogenic'`, which is the PHI-base content itself) and the NCBI tree's
`kingdom` for the host test. One reference backbone instead of two. The consequence is that
residents get scored and reported like anything else, and `pathogen | resident | host` becomes a
column rather than a gate. This is the single most important code change in `05_detect`, and it
is small.

**3. `declared` comes from the curated field, not the LLM.** `declared_pathogens` stores NAMES;
taxids resolve in Python (an LLM recalling a taxid produces confident wrong numbers).
`declared_species_set()` already promotes a declared forma specialis to its species node, which
matters because the values are at f. sp. rank and detections are at species rank.

**4. Fungi only, six hosts.** Oomycetes are out (Leon 2026-10-13): zero of db_v3's 29 oomycete
species pass both `05_filter` floors on any cereal host. Hosts are the four cohort cereals plus
*Hordeum* (monocot outgroup) and *Arabidopsis* (dicot outgroup); `04_classify`'s host call RANKS
rather than thresholds, so the outgroups are what keep a mislabelled non-cereal detectable.

## One open decision, deliberately not taken

`04_assign/kraken_report_analysis.py:321` excludes a whole run from the cryptic rate when
`coinf_intent == "coinf_experiment"`, on the grounds that "a deliberate multi-pathogen experiment
is not cryptic by definition". A correct multi-species `declared_pathogens` set makes that
largely redundant: both detections are marked declared, so nothing is cryptic anyway. But not
identically — `continue` also drops the run from the denominator and from `coinf_runs`.

PRJNA418364 is the live case: 45 runs declaring *A. flavus* **and** *F. verticillioides*.

My view is to drop the exclusion and let the declared set do the work: one mechanism instead of
two, and the exclusion shrinks a denominator the prevalence question needs. It changes a reported
number, so it is Leon's call.

## What must be carried across, not rewritten

`05_detect` is `05_filter` with one change (the pathogen filter becomes a label). Its README
records a documented arc of what was tried and REJECTED — a 1D line residual, a GMM cluster, a
rarefaction test, and automated weight fits that either let coverage dominate or over-rejected
obvious reals. That is knowledge, and it is why the weights are a stated provisional operating
point rather than a fitted model. Keep the README with the code.

Also carried: the `unbiased` f. sp. feature (a detection dominated by one sub-taxon is
penalised), the per-host application of the test, and the two-threshold discipline
(`score >= 0.7` is the artefact filter, `coverage >= 1%` is the abundance floor; they answer
different questions and must not be conflated).

`utilities/` moves wholesale: the BUSCO scan with its tracked cache, the FracMinHash saturation
analysis, the miniprot annotation pipeline, and `split/`. None touches STAT.

## Traps that still apply, carried from CLAUDE.md

- **`_classified.tsv` has NaN scores.** 3,382 of 24,467 rows, all `cls=unclassified`. Test the
  floors POSITIVELY (`>=`); `if score < FLOOR: continue` ADMITS them, because NaN compares False
  against everything. It inflated a species count from 32 to 44 before being caught.
- **Report read counts, not percentages.** A percentage moves whenever the dominant organism's
  abundance moves, conflating two things.
- **Coverage is measured against the DATABASE, not the genome**, so reference incompleteness
  largely cancels but the floor is relative.
- **Coverage is not an artefact filter** and we do not claim it is. The cleanest cross-mapping
  null (obligate rusts on hosts they cannot infect) reaches 13% coverage at its 75th percentile.
- **Polyploid hosts lose k-mer space to the LCA.** Dropping *T. turgidum* recovers the 27.7% of
  *Triticum* minimizers currently parked at the genus node.
- **`--workers 8`, not 16**, for NCBI downloads; `fasterq-dump` dead-second-mate is recurring.
- **Always smoke-test before a full launch.** A 6-run test caught `--run-list` taking the whole
  TSV line as the run ID, which had made every run return `no_local_reads`.

## 01_stat and 02_literature can now go to legacy, and 01_stat's gate is why

Both are superseded as CODE. Everything `01_literature` still reads from them is **data**:

    02_literature/02_text/data/text_cache.jsonl        39 MB, 711 full texts keyed by DOI
    02_literature/03_classify/data/samples.tsv         12 MB, needed only to re-run the
                                                       one-time provenance import
    02_literature/02_text/data/supp_provenance.tsv     224 KB
    02_literature/03_classify/data/cohort_provenance.tsv  328 KB
    02_literature/02_text/data/manual_pdfs/            publisher PDFs
    01_stat/02_fetch/data/stat_cache.jsonl + 03_filter/data/runs.tsv   (triage.py, vestigial)

So the legacy move is: **relocate those data artefacts to a neutral home, repoint
`01_literature/_layout.py` and `curate/import_provenance.py`, then move both modules' code to
`legacy/`.** Two 01_literature scripts are STAT-frame tools and should go with them:
`report/stat_queue.py` and `report/truncation_check.py`.

After that the numbering is correct in one pass: `01_literature`, `02_classify`, `03_align`,
with no collisions and one Setonix migration.

## Build order

1. `seeds/` moves across (already built; it is the only finished piece).
2. `01_refs` + `02_build` — db_v4. Sized: ~530 assemblies, ~5.2 B minimizers, ~28 GB index,
   comfortable on a Setonix standard node. Annotation is the constraint, not size.
3. `03_cohort` — resolve 7,688 runs, fetch the ~4,500 reads not already on scratch.
4. `04_classify` — run it, including the **1,120 sorghum runs that have never been classified**.
5. `05_detect` — carried over, pathogen filter becomes a label, recalibrate nothing yet.

Expect the co-infection rate to MOVE, and only downward on the resident question: adding
competitors can only take reads away from pathogen calls. Wheat 47% is an upper bound until this
is run. Expect inflated Pleosporales calls on sorghum specifically, because sorghum's dominant
residents are exactly what db_v3 lacks.
