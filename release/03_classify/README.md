# 03_classify

Classifies the cohort's reads against the database ladder and decides which assignments are
real organisms. The output is one row per call: this organism is present in this run, at this
coverage and this shape, and this is whether the study that generated the run reported it.

**This module is the SCREEN.** It says which organisms are present and worth characterising.
It is not a prevalence: the cohort cannot support one, so no rate out of here is quoted as a
rate. The finding is carried by `04_align`, which confirms presence read by read and asks
whether the organism is expressing anything.

## The two arms

Both databases hold the plant host sink, so a subtraction between them measures the tier and
not the sink.

| arm | contents | hash | species | clade minimizers |
|---|---|---|---|---|
| `db1p` | USDA scope at genus, annotated, plus the 3 hosts | 57.5 GB | 2,308 | 8,727,920,674 |
| `db2` | `db1p` plus 135 gene-predicted genomes | 59.8 GB | 2,443 | 9,010,023,071 |

`db2 - db1p` is the purpose test: the 135 predictions are the only difference. They are NOT
mostly genus-orphans, and the measurement has to respect that. The v2 worklist is three tiers,
118 + 3 + 14, and 112 of the 135 land in genera `db1p` already represents:

    A   118   in the reference scope already, no gene model, so present and uncallable.
              Giving them CDS both makes them reachable and puts a competitor in their genus.
    B     3   a genus the database did not hold at all. The only tier "do reads move DOWN to
              species" is about.
    C2   14   a competitor for a taxon that had none, so Kraken2's LCA can decline. Its
              intended effect is reads moving UP to the genus node.

So the question db2 - db1p answers is whether the predictions take reads off species that were
holding them for want of a competitor. Pooled across the tiers the aggregate read movement is
upward, which is tiers A and C2 doing their job, not tier B failing.

**`db0` is not an arm.** The copy on disk was built with no plant CDS at all, so `db1p - db0`
would measure the host sink rather than the congener tier. The congener question was answered
separately and does not need the arm: on held-out congeneric *Puccinia* at k=31, about 96% of
reads from an absent species match no reference k-mer at all.

**The 135 predicted species are competitors only and can never be a positive call.** They are
in `db2` to absorb reads that would otherwise be forced onto a congener, and a detection
resting on predicted CDS would be an artefact of our own build. `detect.py` marks them
`competitor_only` from `02_reference/data/curation/build_plan.tsv` and refuses to call them,
in both arms, so the two tables stay comparable.

## Coverage is not comparable between the arms

Every coverage figure divides by the taxon's clade minimizers in the build that produced the
report, and those move when the build moves. *Triticum aestivum* holds 52,100,304 clade
minimizers in `db1p` and 52,096,895 in `db2`: adding the predictions lifts a few of its
minimizers to a higher LCA. The effect is small for the host and is not small for a species
whose new neighbour is close.

So each arm is scored against its own `db_inspect_<arm>.tsv`, each arm gets its own null, and
a coverage figure from one arm is never restated against the other by arithmetic.

## The call

Four tests, all of which must hold, tested in this order and always positively. A missing
value fails the test rather than passing it: `None >= floor` is written out, never
`not (value < floor)`, because a missing value compares False against everything and the
negated form admits it.

| test | what it asks | |
|---|---|---|
| `distinct >= distinct_floor` | is there enough signal to divide at all | absolute |
| `coverage >= coverage_floor` | is the DNA there | relative, per arm |
| `saturation >= saturation_floor` | are the reads spread or piled on a few fragments | shape |
| `fsp_frac <= 0.80` | is one forma specialis carrying the call | inert, see below |

    redundancy = reads * kmers_per_spot / M        M = the taxon's clade minimizers in the arm
    coverage   = distinct / M
    saturation = coverage / (1 - exp(-redundancy))

`kmers_per_spot` comes from the library and not from a constant. `base_count / read_count` in
`supp_data_s2_cohort.tsv` is bases per spot, so a paired library's per-mate length is half of
it, and at k=35 a 51 nt mate yields 17 k-mers against 116 for a 150 nt mate. That is a factor
of 7 on the expectation.

Below `min_redundancy` the saturation statistic is silent and the call rests on the two floors.
Saturation is not scale-free: for small x, `1 - exp(-x)` is approximately x, so coverage
approximates redundancy mechanically and everything scores near 1. It must be read against
redundancy, never as a scalar.

**The absolute floor is not optional, and the arms make the case plainly.** Species ceilings
span 2,921 to 70,294,921 clade minimizers in `db1p`, so a relative floor is 2,308 different
floors spanning four orders of magnitude. 115 species in `db1p` and 121 in `db2` have ceilings
under 100,000, where 1% is fewer than 1,000 k-mers and cross-mapping supplies that for free.
`db2`'s smallest ceiling is 897, where 1% is 9 k-mers.

**The forma specialis test is inert in both arms and the run says so at startup.** Neither
`db1p` nor `db2` contains a single sub-species node, so `fsp_frac` is 0 on every candidate and
the test never fires. It is retained rather than deleted because what it catches, a detection
whose reads pile onto one sub-taxon, is a property of the build: admit one strain-level
assembly and it matters again. A test that silently does nothing is how knowledge rots.

## Where the coverage floor comes from

It is anchored to the ceiling of an off-host null, and the null is rebuilt per arm by
`calibrate_floors.py` because the denominator changed.

The null is obligate rusts scored in runs whose host genus they carry no USDA record for: real
organism DNA that cannot be infecting anything in that sample, so whatever coverage it reaches
is cross-mapping plus deposited material. The floor goes where the null is exhausted, not at a
percentile ratio. Triticale is excluded from it: *x Triticosecale* matches no USDA genus token,
so genuine wheat rusts on it score as off-host and contaminate the tail.

Under `db_v4` that construction put the floor at 2%, where zero of 10,314 null rows cleared it.
**That number does not transfer to these arms and this module does not assume it.**
`calibrate_floors.py` measures; `data/inputs/thresholds.tsv` records what was then decided. The
split is deliberate and is the module's own decide/act discipline: a floor that a script sets
from the data it is about to filter is a floor tuned to a result.

What the floor does and does not exclude, carried over because it is knowledge and not code:

- **Cross-mapping is excluded, and the null is the evidence.**
- **Drift and deposition are not excluded and are not excludable.** Deposited material is real
  DNA off a real genome and spreads k-mers normally, so no abundance or shape statistic
  separates it from infection. That is answered by the host axis and by obligate biotrophy, and
  is disclosed rather than thresholded away.
- **Saturation is not the artefact filter either.** On the `db_v4` null the artefacts scored
  *higher* than genuine calls, because at low N the expectation collapses toward N*L.

## The pipeline

    inputs/thresholds.tsv       the operating point, authored, with the evidence for each floor

    classify.py          -> raw/reports/<arm>/<BioProject>/<Run>.report.tsv   on Setonix
                            raw/reports_manifest_<arm>.tsv   run -> bytes, species rows, TRACKED
    calibrate_floors.py  -> curation/null_offhost_<arm>.tsv  the null, and where it is exhausted
    detect.py            -> outputs/detections_<arm>.tsv     the calls
                            raw/candidates_<arm>.tsv         every candidate with its metrics
    make_supp_tables.py  -> outputs/supp_data_s5_detections.tsv  one row per call, labelled
                            outputs/supp_data_s6_runs.tsv        one row per run
    compare_arms.py      -> outputs/arm_resolution.tsv       db2 - db1p, by assigned rank
    figures/resolution/                                      the purpose test

Labels are joined after detection, never applied before it: `on_wheat`, `is_plant_pathogen` and
`is_wheat_pathogen` are columns on `supp_data_s3_reference.tsv` and are read off a call. The
study's own declaration comes from `supp_data_s1_studies.tsv`. Filtering on who is already
known on wheat would make the finding impossible to make.

`outputs/` is the deliverable and the only thing the figures read, the same contract as
`01_cohort` and `02_reference`. The full candidate table stays in `raw/` on scratch: it is one
row per run and taxon, it is rewritten on every threshold change, and the repo holds code and
small tables.

## Running it

Classification is the only expensive step and it runs on Setonix.

    sbatch release/03_classify/scripts/slurm/classify.slurm --arm db1p --project PRJEB39201
    python release/03_classify/scripts/classify.py --arm db1p --list

**One job per project group per arm, with the hash staged into `/dev/shm`.** The expensive
operation is loading the database, not classifying against it, so a per-run array task would
pay one 58 GB Lustre load per run. A long-lived job stages once and every worker mmaps that
single RAM-resident copy.

**`--memory-mapping` is correct only against `/dev/shm`.** Against Lustre it is catastrophic
and measured: 13 minutes, zero reports, kraken2 pinned at 1.6% CPU, because every hash lookup
becomes a filesystem RPC. `classify.py` does not default it on; the SLURM script decides.

**The memory request is the CPU bill**, so it is sized to the hash and no further: about 100 GB
for either arm, against `db_v4`'s 80 GB, because the hash is 57.5 or 59.8 GB rather than 41 GB.
A request sized for `db_v4` does not fit.

**Every step after classify.py also runs on Setonix**, because the reports are 1,560 files per
arm on scratch and are never copied down. So `detections_<arm>.tsv`, `null_offhost_<arm>.tsv`
and the supplementary tables are written there and reach the local clone through
`tools/repo_sync.sh publish`, by the machine that produced them. That is the same
one-author-per-file rule the manifests already run on, and it is why a partial table should be
deleted rather than left sitting dirty: two authors for one tracked file is what the rule
exists to prevent.

**No depth filter, and no selection of any kind.** Every run on disk is classified. A detection
is valid at any depth; only absence is depth-dependent, and that caveat attaches to a
denominator rather than to a call. A pathogen recovered from a 0.019 Gkmer library had to be
abundant to appear at all, which makes it a stronger observation and not a weaker one.

**Resume is a non-empty report in its final location.** A kraken2 process killed by a timeout
leaves a short, parseable, wrong file, and existence alone would mark it done forever. Reports
are written to a per-job temporary name and renamed on success, and the per-job temporary
directory is a correctness requirement rather than tidiness: two overlapping jobs sharing one
partial path produce a full-size report of interleaved output that no size check can catch.

## Two things no threshold will settle

**Congener clusters are a flag, never a filter.** One run calls fourteen *Alternaria* species
at coverage 5.6 to 7.3% and saturation 0.070 to 0.080; fourteen genuinely co-infecting
congeners would differ wildly in abundance, and the uniformity is the diagnostic. But every one
of the 96 runs a collapse rule would flip is *Puccinia*, 82 of them *striiformis* with
*triticina*, which differ 46.9% against 6.9% coverage and are separable by exact k-mers across
a 16 point ANI gap. Collapsing those discards a real co-infection. It is a property of the
database's congener density, which is deliberate, and `04_align` is what resolves it.

**Presence is not activity.** Spore drift, surface inoculum and mixed samples all deposit real
organism DNA, and coverage cannot separate that from infection. `04_align` confirms presence
read by read before anything is characterised.

## State

**Being built, 2026-11-20.** Nothing classified yet. The inputs are in place: both arms are
built and backed up, their inspect tables are tracked in `02_reference/data/outputs/`, and the
cohort's 1,560 runs are complete on scratch under `01_cohort/raw/reads/`.
