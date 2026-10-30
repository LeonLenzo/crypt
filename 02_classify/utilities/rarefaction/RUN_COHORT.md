# Cohort-wide rarefaction — how to run it

Written 2026-10-30, ready to submit once the outstanding fetch lands. Nothing here has been
run yet. The design and its reasoning are in `select_cohort.py`'s docstring; this file is the
operating procedure and the preconditions.

## Why it is waiting

757 analysed runs have no reads on disk and so are not classified: 295 maize (PRJEB67964),
290 sorghum (PRJEB83049 and PRJNA527782) and 5 others. 585 of them clear 2 Gkmer at a median
of 4.13, so they land squarely inside the selection below and sorghum is the crop that most
needs them. Running before they arrive measures a curve whose central question — does the
wheat result generalise — is still unanswerable.

Jobs queued 2026-10-30, each chained `afterany` on the one before:

    51742589  fetchctl  the three mock-arm controls (running)
    51848379  fetch     5 projects, 2.32 TB
    51848380  fetch     continuation, in case A hits its 24 h wall
    51848381  classify  sweeps whatever landed

## Preconditions, all of which must be true before submitting

1. `squeue -u llenzo` is empty, or at least 51848381 has completed.
2. The reports are synced and `05_detect` has been re-run:

       rsync -a setonix:/scratch/pawsey1168/llenzo/crypt/02_classify/data/bronze/reports/ \
             02_classify/data/bronze/reports/
       python3 02_classify/05_detect/detect.py

3. The supplementary tables are rebuilt, because the selector reads S2 and S3 for the
   analysed flag and the depth:

       python3 02_classify/report/figures/depth_funnel/prep_depth_funnel.py
       python3 01_literature/report/tables/build_s2.py
       python3 01_literature/report/tables/build_s3.py

4. The audit is clean — every analysed run has reads and a report:

       python3 - <<'EOF'
       # the audit from 2026-10-30; it should report 0 in both directions
       EOF

## Submitting

The selector writes a run list; `rarefy.py` turns it into a worklist and the SLURM array
runs one task each. `rarefy.py` is UNCHANGED — its measurement, its two arms, its seed
replicates at the shallow rungs and its refusal to rarefy upward are all as they were for the
five-library version, which is what makes the two curves comparable.

    cd /scratch/pawsey1168/llenzo/crypt

    # cohort arm
    RUNS=$(python3 02_classify/utilities/rarefaction/select_cohort.py)
    python3 02_classify/utilities/rarefaction/rarefy.py --select --runs "$RUNS"
    # read the printed task count before submitting the array

    # control arm, same rung ladder
    CRUNS=$(python3 02_classify/utilities/rarefaction/select_cohort.py --controls)

Submit the array over the worklist, then `rarefy.py --collect`.

## What to check in the result

- **Do the three crops give the same knee?** This is the question the whole exercise exists
  for. If they do, the floor generalises and the caveat currently in the Methods can be
  deleted. If they do not, the floor has to be per-crop and the cohort definition changes.
- **Where does the control arm sit against the cohort arm?** The cohort arm measures recovery
  of calls we would make at full depth, which is self-referential. The control arm measures
  recovery of a known inoculated pathogen. If the control curve is shallower, the cohort arm
  is flattering itself.
- **Does 0.5 Gkmer survive?** It was set at the knee of the five-library curve on 2026-10-30,
  replacing 0.463, which had no derivation recorded anywhere. Expect it to move.

## What this cannot answer

The cohort arm's reference is each run's own full-depth calls, so it can say how much depth is
needed to recover what we would otherwise find and nothing about what full depth itself
misses. Only the control arm speaks to that, and only for three pathogens on wheat.

The floor is on TOTAL depth, not fungal depth. Effective fungal depth is bimodal with a 150x
gap and cannot become an inclusion criterion, because unlike total depth it is an outcome of
the classification rather than a property of the library.
