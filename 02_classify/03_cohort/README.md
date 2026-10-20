# 03_cohort — which runs, where they live, and getting them

Five scripts. **State: BUILT.** The cohort itself is not complete: 1,184 of 7,688 runs are on
disk and classified.

| | |
|---|---|
| `resolve.py` | the plan: `run_list.tsv` (what we want) and `fetch_plan.tsv` (per project) |
| `migrate_reads.py` | one-shot, already run: flat tree -> `reads/<BioProject>/` |
| `fetch_depth.py` | ENA base counts and file URLs, one request per BioProject |
| `fetch_reads.py` | the fetcher |
| `freshen.py` | keeps the cohort tree out of the 21-day /scratch purge |

## The layout: one directory per BioProject

```
data/bronze/reads/<BioProject>/<Run>_1.fastq.gz      paired
                              /<Run>_2.fastq.gz
                              /<Run>.fastq.gz        single-end, OR a paired run whose second
                                                     mate died in fasterq-dump
data/bronze/reads_tmp/<Run>/                         in flight; NOT presence
```

Leon's call, and the reason is navigability at 7,688 runs: a flat directory of ~15,000 fastq
cannot be reasoned about, and a per-project tree makes "is this project down yet" a `ls`.

**Presence is a non-empty file in its final location.** Never a directory, never a name in a
table. That rule is in `_reads.py` because 03_kraken answered the same question three ways and
the three disagreed, and because the 2026-09 CDS tree had 1,997 accession DIRECTORIES against
1,162 real FASTAs, so a directory-existence check counted all 1,997 as held and the build ran
on a third less sequence than it reported.

## The migration, and why it was instant

4.05 TiB / 2,362 files / 25 projects, in **7 seconds**. A rename within one Lustre filesystem
moves no data. `migrate_reads.py` verifies size, inode AND **atime** per rename and aborts on
the first mismatch, because atime is the purge clock and a move that reset it would have
started a fresh 21 days silently in the wrong direction.

39 orphaned `_<RUN>_sra` prefetch directories were cleared at the same time. An `iterdir()`
loop yields those as runs.

## Fetching from ENA, not SRA

ENA serves `fastq.gz` directly, which skips `prefetch` + `fasterq-dump` + `pigz` entirely.
Validated over 400 runs: 100% of pilot runs have fastq.gz available, byte counts are uniformly
0.937 of ours (range 0.907 to 0.948, a compression-level difference), and **zero mate-count
disagreements**.

- **`DEFAULT_WORKERS = 24`, not 8.** The "8 not 16" rule is an *NCBI* rule about connection
  contention. EBI is latency-bound: each stream sits near 0.7 MB/s regardless of how many are
  open, so concurrency is the only lever.
- **`pick_files()` must choose.** ENA lists a third orphan-reads file beside `_1`/`_2` for some
  runs. `_1`/`_2` win.
- **Timeouts scale with size:** `budget = max(900, nbytes/(2 MB) + 300)`. A flat timeout broke
  both the genome download and classification previously.
- Files land in `reads_tmp/<Run>/` and are renamed into place only after a byte-count check and
  `gzip -t`.

## Pilots, and what they are for

Most of the volume is a few large physiology projects. `resolve.py` flags any project over
`PILOT_THRESHOLD = 150` eligible runs or `PILOT_GIB = 300` estimated volume and fetches
`PILOT_RUNS` from it first. If a project yields no fungal reads it will not be included, so
fetching all of it first would be buying the answer at full price.

- **`--limit` is a smoke test, not a sample.** It takes the first N in accession order, which
  is submission order. Change the sample size with `--pilot N`, which redraws against
  `SEED = 20261019`.
- `ELIGIBLE_STRATEGIES = {"RNA-Seq"}` drops 55 AMPLICON and 66 ncRNA-Seq runs.
- **`fetch_mode` is not a denominator.** It records what the plan intended. What was actually
  classified is counted from the reports by `05_detect`.

## Three projects flagged, deliberately not fixed here

PRJNA383416 (740 runs / 1 BioSample), PRJNA721299 (36 / 1) and PRJNA1181034 (12 / 1) collapse
to one BioSample. Any per-sample statistic over them is wrong, and the flag is carried rather
than silently deduplicated, because which is right depends on the question.

## The depth confound, measured

Gkmer, the composite a detection actually depends on, is `bases x (L-34)/L`: 0.33 at 51 nt, 0.66
at 101, 0.77 at 150. Across the cohort depth spans **300x and aligns with crop**: 60% of rice
is 51 nt / 0.057 Gb. **Rank on Gkmer, never on library layout**, and remember the coverage
denominator itself varies 24,000x across taxa.

Rarefaction settled the related question: **read length costs nothing beyond Gkmer.** That is a
negative result and it is useful, because it makes Gkmer a sufficient single statistic.

## freshen.py, and the bug it exists because of

`/scratch` is Lustre mounted without `noatime`, so the purge clock runs on atime and reading
4 KB refreshes a file as well as reading all of it: a full pass over 2,442 fastq is ~10 MB of
I/O, not 4 TB.

The previous keep-alive was `03_kraken/03_select/freshen.py`, which hardcoded the flat
pre-migration directory. After `migrate_reads.py` drained that tree the script still ran, still
reported, and still exited 0, against a tree the cohort had left. **A keep-alive that silently
covers the wrong files is worse than none, because the report looks healthy.** This version
takes its root from `_layout.READS`.

The flat tree still holds ~9.1 TB of non-cohort STAT-gated runs. They are **not** freshened and
are meant to expire; every one is a public SRA accession.
