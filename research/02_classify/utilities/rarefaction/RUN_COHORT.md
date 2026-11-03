# Cohort-wide rarefaction — how to run it

Written 2026-10-30, ready to submit once the outstanding fetch lands. Nothing here has been
run yet. The design and its reasoning are in `select_cohort.py`'s docstring; this file is the
operating procedure and the preconditions.

## Why it is waiting

757 analysed runs have no report and so contribute nothing: 295 maize (PRJEB67964), 290
sorghum (PRJEB83049 and PRJNA527782) and 5 others. 585 of them clear 2 Gkmer at a median of
4.13, so they land squarely inside the selection below and sorghum is the crop that most needs
them. Running before they arrive measures a curve whose central question -- does the wheat
result generalise -- is still unanswerable.

**State 2026-11-02.** The reads are now almost all down: the weekend fetches landed 740 runs,
taking the cohort from 4,121 to 4,861 runs with reads, and 3,524 of the 3,541 analysed runs
have them. 17 do not and are not worth waiting for: 15 are transient ENA failures in
PRJEB67964 (`curl` 56 on an SSL short read, and three 403s) and 2 are MARPLE runs in
PRJEB31334. A fetch retry is queued for them but the rarefaction does not depend on it.

None of those 740 are classified yet. The weekend's classify job died in 9 seconds because it
was submitted without `--export=ALL,PROJECT=...`, which `classify.slurm` requires and which
has no default. **Always pass PROJECT.** Two jobs are queued today, fetch then classify, the
fetch first so a single classify pass covers all 757:

    51974436  fetch     the 17 stragglers, 3 h
    51974437  classify  all five projects, 5 h, afterany on the fetch

The walltimes are cut from the scripts' 24 h because a 24 h request crosses tomorrow's 07:00
maintenance reservation and sits `PENDING (ReqNodeNotAvail)` indefinitely rather than starting.
24 h was never needed: 8 workers classify concurrently at an effective 6 s per run, measured
over 1,258 runs in job 51651279, so 740 runs is ~1.2 h plus a 200 s hash stage.

## Preconditions, all of which must be true before submitting

1. `squeue -u llenzo` is empty, or at least the classify job has completed.
2. The reports are synced and `05_detect` has been re-run. **Filter the rsync to reports.**
   An unfiltered pull of the parent directory once dragged down 18 GB of control reads:

       rsync -a --include='*/' --include='*.report*' --exclude='*' \
             setonix:/scratch/pawsey1168/llenzo/crypt/02_classify/data/bronze/reports/ \
             02_classify/data/bronze/reports/
       python3 02_classify/05_detect/detect.py

3. The supplementary tables are rebuilt, because the selector reads S2 and S3 for the
   analysed flag and the depth:

       python3 02_classify/report/figures/depth_funnel/prep_depth_funnel.py
       python3 01_literature/report/tables/build_s2.py
       python3 01_literature/report/tables/build_s3.py

4. The audit is clean -- every analysed run has a report:

       python3 02_classify/utilities/rarefaction/audit_analysed.py

   It exits non-zero while any analysed run is unclassified, and prints the outstanding
   BioProjects as a `PROJECT="..."` string ready to paste into the resubmit. Pass
   `--reads-list` to separate "needs a fetch" from "needs a classify". The reports it counts
   as outside the cohort are expected and are kept deliberately: a run excluded on design or
   depth was often classified before the criterion existed, and two criteria have already been
   revisited.

## Submitting

The selector writes a run list; `rarefy.py` turns it into a worklist and the SLURM array
runs one task each. `rarefy.py` is UNCHANGED — its measurement, its two arms, its seed
replicates at the shallow rungs and its refusal to rarefy upward are all as they were for the
five-library version, which is what makes the two curves comparable.

    cd /scratch/pawsey1168/llenzo/crypt

`select_cohort.py` now REFUSES to emit a list while any selected run is unclassified, and
that guard is not theoretical: on 2026-11-02, 585 of the 2,908 selected runs -- a fifth of the
denominator -- were the outstanding maize and sorghum ones. Each would have contributed zero
recovered detections at every rung, which is indistinguishable from a run that genuinely lost
them, so the curve would have been pulled down by 20% with nothing in the output to show it.
Pass `--no-require-report` only to inspect a selection, never to measure one.

    # cohort arm
    RUNS=$(python3 02_classify/utilities/rarefaction/select_cohort.py)
    python3 02_classify/utilities/rarefaction/rarefy.py --select --runs "$RUNS"
    # read the printed task count before submitting the array

    # control arm, same rung ladder
    CRUNS=$(python3 02_classify/utilities/rarefaction/select_cohort.py --controls)

Submit the array over the worklist, then `rarefy.py --collect`.

## The selection splits on whether there is anything to recover

Found 2026-11-02, after the 757 runs were classified and before anything was submitted. Of the
2,908 runs that clear the 2 Gkmer ceiling, **1,593 (54.8%) have no fungal detection at full
depth**, and they are not spread evenly: wheat 1,010 of 1,012 carry one, against maize 215 of
974 and sorghum 90 of 922.

This matters because the cohort arm's reference is each run's own full-depth calls. A run with
none contributes 0/0 at every rung. Selecting on depth alone therefore produced a denominator
balanced across crops and a CURVE still 77% determined by wheat, which is the precise
assumption the cohort-wide run exists to retire. The five-library problem would have survived
the fix meant to remove it.

So the selector now emits two kinds of row, and the measurement reads them separately:

- **recovery rows**, 1,315 runs with at least one full-depth call (wheat 1,010, maize 215,
  sorghum 90). Sensitivity per rung. 90 sorghum and 215 maize runs is not parity with wheat,
  but it is enough to ask whether the knee differs by crop, which 0 was not.
- **spurious-gain rows**, a 202-run sample of the zero-detection runs. A run with no full-depth
  call that GAINS one at a shallow rung is a false discovery, so these give a false-discovery
  rate per rung to set against the recovery rate. Keeping them is a better use than discarding
  them, and a few hundred pin the rate down as well as 1,593 would at half the compute.
  `--all-zero` keeps all of them; `--null-sample 0` measures recovery only.

The sample is per-crop and deterministic (sorted, then strided) so it cannot itself be
crop-skewed and needs no seed recorded.

Cost after the split: 30,441 tasks, 10,216 Gkmer, roughly 4 h on 64 cores, against 59,556
tasks and 8.3 h before it.

## What to check in the result

- **Do the three crops give the same knee?** This is the question the whole exercise exists
  for. If they do, the floor generalises and the caveat currently in the Methods can be
  deleted. If they do not, the floor has to be per-crop and the cohort definition changes.
- **Where does the control arm sit against the cohort arm?** The cohort arm measures recovery
  of calls we would make at full depth, which is self-referential. The control arm measures
  recovery of a known inoculated pathogen. If the control curve is shallower, the cohort arm
  is flattering itself.
- **Where does the false-discovery rate cross the recovery rate?** The floor should sit where
  recovery is still high and spurious gain is already low. This is the first version of the
  curve that can show both at once.
- **Does 0.5 Gkmer survive?** It was set at the knee of the five-library curve on 2026-10-30,
  replacing 0.463, which had no derivation recorded anywhere. Expect it to move.

## What this cannot answer

The cohort arm's reference is each run's own full-depth calls, so it can say how much depth is
needed to recover what we would otherwise find and nothing about what full depth itself
misses. Only the control arm speaks to that, and only for three pathogens on wheat.

The floor is on TOTAL depth, not fungal depth. Effective fungal depth is bimodal with a 150x
gap and cannot become an inclusion criterion, because unlike total depth it is an outcome of
the classification rather than a property of the library.
