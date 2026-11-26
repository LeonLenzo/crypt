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
| | **the host is exempt from both abundance floors below** | |
| `coverage >= coverage_floor` | is the DNA there | relative, per arm |
| `saturation >= saturation_floor` | is it abundant relative to what this depth allowed | depth-normalised |
| `fsp_frac <= 0.80` | is one forma specialis carrying the call | inert, see below |

    redundancy = reads * kmers_per_spot / M        M = the taxon's clade minimizers in the arm
    coverage   = distinct / M
    saturation = coverage / (1 - exp(-redundancy))

`kmers_per_spot` comes from the library and not from a constant. `base_count / read_count` in
`supp_data_s2_cohort.tsv` is bases per spot, so a paired library's per-mate length is half of
it, and at k=35 a 51 nt mate yields 17 k-mers against 116 for a 150 nt mate. That is a factor
of 7 on the expectation.

**Saturation is a depth-normalised ABUNDANCE floor and NOT a shape test, measured 2026-11-24.**
It was specified as a shape test and cannot be one. Over db2's 430,820 candidates, every row
saturation alone rejects has redundancy >= 0.534 and every row the coverage floor alone
rejects has redundancy <= 0.480: the two partition at r = 0.51, where `0.05*(1-e^-r)` crosses
0.02, with no overlap. Above r = 3 the median `saturation / coverage` is 1.007 and above 10 it
is 1.000, because `1-exp(-r)` has gone to 1. So it only ever acts where the shape information
has vanished, and where it would carry shape information the coverage floor has already bound.

The two are therefore ONE floor, `coverage >= max(coverage_floor, saturation_floor*(1-e^-r))`,
and read that way it is worth having: it asks for 5% of what the depth ALLOWED, so a shallow
library gets a proportionally lower bar and a deep one a flat 5%. **The genuine shape
question, where in a transcript the reads fell, needs positions this module does not have and
belongs to `04_align`.** The two do not overlap, which is why this stays here.

Below `min_redundancy` saturation is silent and the call rests on the other two floors.

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

    figures/samples/      samples.png        the per-run raster, one column per run
    figures/runs/         runs_1..4.png      four of those columns opened up
    figures/multiplicity/ multiplicity.png   co-infection multiplicity, per BioProject
    figures/declared/     declared.png       what each study declared
    figures/resolution/   resolution_reads.png  where the reads went, in reads
                          resolution.png        the same, run by run, SUPPLEMENTARY

Each figure directory is a `prep_*.py` that DECIDES, a plotting script that DRAWS, and a
`caption.py` that writes the caption from the tables between them. `figures/multiplicity/`
and `figures/resolution/` run locally off `outputs/`; `figures/runs/` runs on Setonix,
because it reads the candidate table.

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

## A co-occurrence network is REJECTED, and the measurement is why

It was proposed as the module's second figure on 2026-11-25, measured before being drawn,
and it is not being built. Over the 1,347 drawn field-wheat runs, the marginal odds ratio
for the commonest wheat-pathogen pairs against the Mantel-Haenszel odds ratio stratified by
BioProject, where `k` is the number of projects contributing an informative stratum:

    pair                                  marginal OR   within-project OR    k
    P. striiformis + P. triticina                0.28                4.43    4
    P. graminis    + Z. tritici                  2.61                0.26    2
    M. nivale      + P. striiformis              1.08                0.15    4
    P. triticina   + Z. tritici                  1.80                1.10    8
    B. graminis    + P. triticina                2.23                1.46    7
    P. striiformis + Z. tritici                  0.84                1.68    7
    B. graminis    + P. striiformis              0.41                0.53    6
    B. graminis    + Z. tritici                  1.24                1.27    9

**Three of the eight reverse direction outright** and the two that keep their sign sit near
1. Every top pair is about half concentrated in a single BioProject and spread over 2 to 10
projects, so a test over runs reports a p-value for a handful of independent units.

This is worse in kind than the prevalence error the module already refuses, not just worse
in degree: a network LOOKS like a claim that two organisms associate, when what it draws is
where the surveys sampled, and a reviewer who conditions on project gets the opposite
answer from the figure. If a pairwise view is ever wanted it has to be within-project or
stratified, and 13 nodes is too few for a network to say anything a matrix would not say
more honestly.

**`figures/multiplicity/` is what replaced it.** Co-infection multiplicity per project
carries the same subject with no Simpson exposure, because each project is one observation
of its own rate: 537 of the 1,347 runs (39.9%) carry two or more wheat pathogens, and across
the 19 projects with at least 10 runs the rate runs from 0.0% to 84.6%. The spread between
studies is the finding, and it is not a depth artefact; the panel's second column shows the
projects with no co-infection sitting among the deepest in the cohort.

The old co-occurrence design is still on disk in the research tree and reads as live work.
It is not.

## State

**Both arms are classified, 1,560 / 1,560, and the module has run end to end** (2026-11-20).

| | db1p | db2 |
|---|---|---|
| candidates † | 409,459 | 430,830 |
| calls, the host aside | 6,629 | 6,380 |
| host calls | 1,488 | 1,488 |
| rejected on the absolute floor | 392,538 | 378,982 |
| rejected on coverage | 8,648 | 7,545 |
| rejected on saturation | 261 | 245 |
| competitor_only | n/a | 36,295 rows, 35.2 M reads |
| runs with a fungal call | 1,406 | 1,403 |
| off-host null | 8,234 rows, max 1.964% | 8,233 rows, max 1.956% |

Call counts are read from the current `outputs/detections_<arm>.tsv`, which `detect.py`
rewrote on 2026-11-25 when `host_coverage_floor` went live (host calls 1,560 -> 1,488).
**† The candidate and rejection rows predate that re-run**: they come from `detect.py`'s own
log on Setonix and the candidate table never leaves scratch, so they cannot be refreshed from
the repo. Re-read them from the next run's log rather than trusting them to the digit.

The coverage floor is **anchored on each arm's own null** and both land on 2%, independently
and on different denominators, which is where db_v4's null also sat (max 1.970%). No null row
clearing the distinct floor reaches 2% in either arm.

**db2 removes 249 calls net.** Residents go 4,402 -> 4,209, wheat pathogens 2,185 -> 2,159
and calls labelled a pathogen not recorded on wheat 42 -> 12; the host is 1,488 in both arms.
The clearest case is
*Stemphylium lycopersici*, a tomato pathogen called in 21 runs under `db1p` and in 1 under
`db2`: its reads were Pleosporaceae reads with nowhere better to go, and with *Alternaria
triticina* in the database they go there instead, taking it from about 0.07 coverage to 0.013.

**That case is also the module's sharpest open question.** *A. triticina* is a wheat pathogen
that was in the reference scope and could never be called for want of a gene model, which is
exactly what tier A exists to fix. But predicted CDS is a competitor only, so `db2` can now
remove the false *Stemphylium* call without being able to make the true *Alternaria* one. The
rule is right for a novelty claim resting on our own gene prediction; whether it should hold
for a tier-A taxon whose assembly is public and whose prediction passed BUSCO is a decision
nobody has taken yet. `04_align` is where such a call would be confirmed from reads.

**Declarations are complete, 31 of 31 BioProjects** (2026-11-25), so no detection reads
`unknown` any more. The nine projects that declare no pathogen are all virus surveys, 114
field-wheat runs, and `study_target` is a curated column on `supp_data_s2_cohort.tsv`: a study
that declared nothing because it was looking for a virus is not a study whose curation is
missing.
