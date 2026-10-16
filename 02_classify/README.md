# 02_classify — 01_refs DONE, 02_build RUNNING

This README began as the design, written before the code, and the design sections below still
read that way. Status as of 2026-10-16:

| step | state |
|---|---|
| `seeds/` | moved across unchanged from `03_kraken`, already built |
| `01_refs` | **DONE.** 4,286 annotated assemblies / 1,547 species, measured against NCBI, not estimated |
| `02_build` | `fetch_cds.py` + `build_db.py` written; the fetch is running on Setonix |
| `utilities/saturation` | carried across, repointed, demoted to a QC record |
| `03_cohort`, `04_classify`, `05_detect` | not started |

Not started also means not designed past the sections below. Two pieces of `01_refs` are
deliberately unfinished: the **miniprot** pass over the 2,138 species that have assemblies but
no gene model, and **host genomic** references (`seeds/hosts.tsv`; `kraken_hosts.py` has not
moved). Host CDS *is* in the candidate set — 43 assemblies across the six species — which is
what db_v4 actually needs, because these are RNA-seq libraries.

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

`pathogens` is PHI-base's 632 fungal taxa. db_v3 holds 877 fungal species and db_v4 holds
**1,547** (measured, not the ~1,100 first guessed here), **including every resident taxon the
2026-10-13 scope work added**. So the gate would now skip rather more than it did. None of
*Cladosporium*, *Aureobasidium*, *Papiliotrema* or the phyllosphere yeasts is a PHI-base
pathogen, so `detect.py` would silently skip all of them. That is not a stale path; it is a
definition of "pathogen" narrower than the database and narrower than the question.

## What this module does

Three jobs, in order. Its output is one sentence per sample: *these organisms are present, at
this confidence, and this many of them were never reported by the study that generated it.*

## Layout

Laid down 2026-10-15. **The code layout and the data layout are separate decisions, and
03_kraken got the second one wrong.**

```
02_classify/
  README.md     this file
  _layout.py    every path, and which layer may be regenerated. NOT _paths.py: the repo root
                owns that name and a second module of it shadows the root for whichever
                script imports first, which has already cost a lost constant once.
  seeds/        MOVED from 03_kraken/seeds, unchanged. AUTHORED and irreplaceable, so it sits
                OUTSIDE data/ and is tracked by default. The membership authority: what is in
                and out of the CDS database, authored vs generated, and THE SCOPE RULE
                (plant-associated, ANY host; never scope by the cohort's own hosts).
  data/         everything DERIVED, one medallion layer per directory. Ignored except gold/.
    bronze/     immutable external captures: the CDS tree, the reads, raw kraken reports,
                NCBI assembly summaries, BUSCO lineages. ~200 GB, on Setonix scratch.
    silver/     generated tables: ref_candidates, db_inspect, run_list, host_calls,
                _classified. Rebuildable from bronze + seeds; hand-editing one is a mistake.
    gold/       small, generated, TRACKED: detections, coinfection, host_summary. The numbers
                a reader wants to see change over the life of the project.
  01_refs/      NEW.      select assemblies from seeds/ and download their CDS
  02_build/     moved.    kraken2-build over the CDS, + the inspect table
  03_cohort/    REWRITTEN. resolve the cohort to a run list and fetch reads. NO STAT.
  04_classify/  REWRITTEN. run kraken2; call the host from reads
  05_detect/    moved.    the detection filter: real signal vs classifier artefact
  utilities/    moved.    busco, saturation, annotate (miniprot), split
  slurm/        adapted.  the Setonix job chain
  logs/         one place, not six
```

**The step directories hold `.py` and nothing else.** That is the change from 03_kraken, where
every numbered step carried its own `data/`, `logs/` and `figures/`: twenty-odd directories for
eight scripts, with the tables split across five sibling `data/` dirs so nothing said which was
authoritative or which could be rebuilt. That is the soil the seven hardcoded paths grew in, and
it is how a live gate came to hide in a string literal.

**The numbering stays, because here the order is real.** Refs before build, build before
classify, classify before detect, and the SLURM chain enforces exactly that with
`--dependency=afterok`. When a job dies overnight the first question is which step. 01_literature
dropped its numbering for the opposite reason: nothing there runs in a fixed order, the same
entry point is re-run dozens of times, and its growth item is one resolver per awkward study, so
grouping by role localised the growth. Numbering it would have been a lie about the order. What
IS carried across from that reorganisation is the data discipline, which is the half that
applies: *if losing a file costs a refetch or a rerun it is bronze or silver; if it costs
re-deciding what belongs in the database it is seeds.*

Step names say the job rather than the tool, like the module name. Renamed from
search/build/select/assign/filter on the same principle (`04_kallisto` becomes `03_align`).

## The input contract, stated once

Everything comes from three places and nothing else:

    01_literature/release/cohort_samples.tsv    Run, BioProject, crop, host, setting,
                                                sampling_selection, tissue, location, date
    01_literature/data/silver/runs.tsv          Run -> BioSample for all 33,958 runs, AND
                                                `declared_pathogens` + its _source in the same
                                                table (columns 22-23): what the STUDY declared,
                                                COMPLETE 2026-10-15, 41/41 projects. The ungated
                                                replacement for the STAT run table: it covers
                                                7,688/7,688 cohort runs over 6,222 BioSamples
    02_classify/seeds/                          what the database contains
    _refdata/holobase.db                        pathogen / host / guild LABELS
    Setonix                                     CDS tree + reads, keyed by ACCESSION

Two files, not three: the declared set is a pair of columns in `runs.tsv` rather than a table of
its own, so the run resolution and the declaration arrive together and cannot disagree about
which runs exist. `_layout.py` names both, and `LIT_CURATION` / `LIT_PROVENANCE` alongside them,
because a declared value is an evidence claim traceable to a quote and not a lookup.

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
2. `01_refs` + `02_build` — db_v4. **DONE and RUNNING respectively.** The sizing originally
   written here (~530 assemblies, ~5.2 B minimizers, ~28 GB index) was an estimate from
   holobase availability counts and was wrong by roughly 8x on assemblies, which is the third
   time an estimate from availability has missed: see CLAUDE.md, "size the next build from a
   `--select-only` run, not from holobase counts". Measured instead:

   | | db_v3 | db_v4 |
   |---|---:|---:|
   | annotated assemblies | 2,063 candidates / 2,081 built | **4,286** |
   | fungal species | 877 | **1,547** |
   | CDS sequence | ~25.5 GB | ~90 GB |
   | distinct minimizers | 4,800,752,701 | **not estimated** |
   | `hash.k2d` | 27.4 GB | **not estimated** |

   The two blanks are deliberate. Every projection this project has made from availability
   counts has been high, so the minimizer count comes from `kraken2-inspect` after the build
   and from nowhere else. The linear upper bound is ~17 B minimizers and ~97 GB of hash, which
   still fits one 250 GB `work` node, and `highmem` offers 1 TB if it does not. What is NOT
   acceptable is `--max-db-size`: it makes Kraken2 subsample the minimizer set to fit, silently
   discarding k-mer space and invalidating every coverage denominator in `05_detect`.

   Of the 4,286, **1,753 were already downloaded for db_v3** and sit in Acacia's project
   bucket, so the NCBI gap is 2,533. Annotation remains the constraint rather than size: 2,138
   further species have assemblies with no gene model and are miniprot candidates.
3. `03_cohort` — resolve 7,688 runs, fetch the ~4,500 reads not already on scratch.
4. `04_classify` — run it, including the **1,120 sorghum runs that have never been classified**.
5. `05_detect` — carried over, pathogen filter becomes a label, recalibrate nothing yet.

Expect the co-infection rate to MOVE, and only downward on the resident question: adding
competitors can only take reads away from pathogen calls. Wheat 47% is an upper bound until this
is run. Expect inflated Pleosporales calls on sorghum specifically, because sorghum's dominant
residents are exactly what db_v3 lacks.
