# 04_classify — Kraken2 against db_v4

One script. **State: BUILT and RUN**, 1,184 runs classified.

```
kraken2 --db data/silver/db/db_v4 --confidence 0.15 --minimum-hit-groups 3 \
        --report-minimizer-data
```

## Decisions

- **`--report-minimizer-data` is not optional.** The `distinct` minimizer column is the
  numerator of both coverage and saturation. Without it `05_detect` has nothing to divide, and
  the failure is a column-count change rather than an error.
- **`--confidence 0.15` with `--minimum-hit-groups 3`.** Carried from db_v3 unchanged, so a
  rate difference between builds is attributable to the database rather than to the caller.
- **Reports are written per BioProject**, matching the reads tree:
  `data/bronze/reports/<BioProject>/<Run>.txt`.
- **Write to `.partial`, rename on success.** `done()` tests
  `report.exists() and size >= MIN_REPORT_BYTES` (200), because a kraken2 process killed by a
  timeout leaves a short, parseable, wrong file. A flat timeout silently truncated 224 files
  in 2026-09 and the build consumed them without complaint.
- **Report counts, not percentages.** A count supports a claim about genome-level coverage of a
  secondary organism; a percentage moves whenever the dominant organism's abundance moves,
  conflating two things.

## Parsing a report, which is where the traps are

`_reports.py` owns this, ported by copy from `03_kraken/04_assign/kraken_report_analysis.py`
rather than imported, so the retired module can be retired.

- **Six-column and eight-column reports are told apart by COUNT, not by content.**
  `--report-minimizer-data` adds two columns in the middle. Blind indexing reads the minimizer
  count as the rank code, which does not match any known rank, which drops **every species row**
  while exiting 0.
- **Host removal is unnecessary.** Tested directly: PST read counts were ratio 1.00x raw versus
  HISAT2-split across all five ground-truth runs. Removal changes the percentage, not what
  Kraken2 finds.
- **The host call is a ranking, not a threshold**, because every plant library has a host.
  Ranking by reads agrees with the LLM 91.2% of the time. It is also **one-directionally biased
  in rust surveys**, so the barley in any host chart is mostly wheat. Keep coverage as a flag:
  a `low_coverage` host call locates a database gap, not a weak library.
