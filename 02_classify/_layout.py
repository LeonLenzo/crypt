#!/usr/bin/env python3
"""Every path in 02_classify, and which layer may be regenerated.

The module keeps its numbered steps because the order here is REAL: refs before build, build
before classify, classify before detect, and the SLURM chain enforces exactly that with
`--dependency=afterok`. 01_literature dropped numbering because nothing there ran in a fixed
order; that is not true of a pipeline.

What 01_literature's 2026-10-06 reorganisation got right, and what IS carried across, is the
separation of the DATA layout from the CODE layout. 03_kraken gave every numbered step its own
`data/`, `logs/` and `figures/`: twenty-odd directories for eight scripts, with the tables split
across five sibling `data/` dirs so nothing said which one was authoritative or which could be
rebuilt. That is the soil the seven hardcoded paths grew in, and it is how a live gate came to
hide in a string literal (`01_stat/03_filter/data/runs.tsv`, which silently dropped 84% of the
cohort). So here the step directories hold `.py` and nothing else, and all data is one tree:

    seeds/            AUTHORED and irreplaceable: the membership authority for the database.
            Which fungi and hosts are in, which are out and why, authored vs generated, and
            THE SCOPE RULE (plant-associated on ANY host; never scope by the cohort's own
            hosts). TRACKED by default. Half of it is a reproducible holobase extraction that
            `build_seeds.py` rewrites; the other half must never be overwritten. seeds/README.md
            says which is which.
    data/bronze/      immutable captures of somebody else's data: NCBI assembly summaries, the
            downloaded CDS tree, the BUSCO lineage caches, the raw per-run Kraken2 reports.
            Deleting one costs a refetch or a recompute and nothing else. Ignored.
    data/silver/      generated tables: the selected assemblies, the inspect table, the resolved
            run list, the host calls, `_classified.tsv`. Rebuildable from BRONZE + SEEDS, so
            ignored. Hand-editing one is always a mistake.
    data/gold/        small human-readable outputs: the detections, the co-infection summaries,
            the per-host rates. Generated, but TRACKED, because these are the numbers a reader
            wants to see change over the life of the project.

The rule in one line: **if losing a file costs a refetch or a rerun it is BRONZE or SILVER; if
it costs re-deciding what belongs in the database it is SEEDS.**

Most of the heavy data lives on Setonix scratch rather than in the checkout, so the repo-side
tree is mostly small tables. That is an argument for getting the split right, not for skipping
it: the same constants resolve on both machines because the module is at the same relative path.

Named `_layout` and NOT `_paths`, for the reason 01_literature found the hard way: the repo root
owns `_paths.py` for cross-module artefacts, and a second module of that name shadows it for
whichever script imports first. That broke review.py once, which resolved `_paths` to the root
module and lost its SEEDS constant. Scripts here need both, so the two names must stay distinct.

Scripts must take paths from here rather than deriving them from `Path(__file__).parent`. The
~50 hardcoded `HERE / "data"` expressions in 01_literature were the actual debt: the layer split
existed in the directory listing but not in the code, which made moving any file unsafe and let
a script write to the wrong layer unnoticed.
"""

import sys
from pathlib import Path

MODULE = Path(__file__).resolve().parent      # 02_classify/
ROOT = MODULE.parent                          # phd/01-review/

sys.path.insert(0, str(ROOT))
from _paths import REFDATA  # noqa: E402  the shared reference store, defined once at the root

# Two top-level ideas, as in 01_literature: seeds/ is AUTHORED, data/ is DERIVED. seeds/ stays
# OUT of data/ deliberately, because data/ means "delete it and it rebuilds" and that is false
# of nothing else in the module. Merging it back in would put the irreplaceable files inside the
# one directory every .gitignore and backup rule in this repo treats as disposable.
SEEDS = MODULE / "seeds"
DATA = MODULE / "data"
BRONZE = DATA / "bronze"
SILVER = DATA / "silver"
GOLD = DATA / "gold"
LOGS = MODULE / "logs"
# Figures are CODE-adjacent, not data: the repo policy is code, small figures and manifests,
# so report/ sits outside data/ and is tracked by default.
REPORT = MODULE / "report"
FIGURES = REPORT / "figures"

# ── the input contract: three places and nothing else ───────────────────────────────────────
# Stated in README.md and enforced here. Everything this module consumes from outside itself is
# named in this block, so a change in 01_literature's layout breaks in one file rather than in
# seven string literals.

# The cohort, final: 7,688 runs / 6,222 BioSamples / 41 BioProjects. Run, BioProject, crop,
# host, setting, sampling_selection, tissue, location, collection_date, lat, lon + the _rule
# column behind each curated value.
COHORT_SAMPLES = ROOT / "01_literature/release/cohort_samples.tsv"
COHORT_PROJECTS = ROOT / "01_literature/release/cohort_bioprojects.tsv"

# Run -> BioSample for ALL 33,958 runs, and `declared_pathogens` + `declared_pathogens_source`
# in the same table (columns 22-23). This is the UNGATED replacement for the STAT run table:
# it covers 7,688/7,688 cohort runs, where the STAT-retained set covered 1,227 (16%). Nothing in
# this module may resolve a run through any other table; that is the whole point of the rewrite.
LIT_RUNS = ROOT / "01_literature/data/silver/runs.tsv"

# What each STUDY declared, as curated names traced to a quote in 01_literature's provenance.
# COMPLETE as of 2026-10-15: 41/41 projects. Read from LIT_RUNS above; the seed table is the
# authority behind it and is listed here because `declared` is an evidence claim, not a lookup.
LIT_CURATION = ROOT / "01_literature/seeds/curation.tsv"
LIT_PROVENANCE = ROOT / "01_literature/seeds/provenance.tsv"

# ── seeds: the membership authority ─────────────────────────────────────────────────────────
# Authored, never written by a script:
SOURCES = SEEDS / "sources.tsv"
HOSTS = SEEDS / "hosts.tsv"
RESIDENTS = SEEDS / "residents.tsv"
EXCLUSIONS = SEEDS / "exclusions.tsv"
# Generated by build_seeds.py, each from a different kind of truth:
FUNGI_USDA = SEEDS / "fungi_usda.tsv"            # holobase: WHICH fungi are plant-associated
FUNGI_ANNOTATED = SEEDS / "fungi_annotated.tsv"  # NCBI: which of them have a gene model, so CDS

# ── bronze: external captures ───────────────────────────────────────────────────────────────
# The CDS tree is keyed by ACCESSION on both machines, which is why no data moves on Setonix
# for this rewrite: the new module points at the same directories as the old one.
CDS = BRONZE / "cds"
CDS_PATHOGEN = CDS / "pathogen"
CDS_HOST = CDS / "host"
READS = BRONZE / "reads"
# Raw Kraken2 output, one report per run. Large and regeneratable.
REPORTS = BRONZE / "reports"
NCBI_SUMMARY = BRONZE / "ncbi"          # assembly summaries, as fetched
BUSCO_CACHE = BRONZE / "busco"          # lineage downloads + the scan cache
# FracMinHash sketches, one .npy per accession. A cache: expensive to compute (one pass over
# every CDS FASTA), worthless to keep if the sketch parameters change, so BRONZE and ignored.
SKETCHES = BRONZE / "sketches"

# ── silver: generated tables ────────────────────────────────────────────────────────────────
REF_CANDIDATES = SILVER / "ref_candidates.tsv"   # 01_refs: the assemblies chosen, and why
# 02_build: one row per WANTED accession with where its CDS actually came from (disk, acacia
# or ncbi) and how many bytes arrived. This is the table that makes the build auditable: the
# 2026-09 tree had 1,997 accession DIRECTORIES on scratch but only 1,162 non-empty FASTAs, and
# a directory-existence check counted all 1,997 as held. Presence is a non-empty file, always.
CDS_MANIFEST = SILVER / "cds_manifest.tsv"
HOST_CANDIDATES = SILVER / "host_candidates.tsv"
DB = SILVER / "db"                               # the Kraken2 index itself (db_v4)
# taxid -> minimizers: the denominator every coverage figure is measured against. The root
# `_paths.DB_INSPECT` is db_v3's and still exists; a script that needs both must alias one,
# because the names are identical and only the module prefix separates them.
DB_INSPECT = SILVER / "db_inspect.tsv"
RUN_LIST = SILVER / "run_list.tsv"               # 03_cohort: the 7,688 runs, resolved
READS_MANIFEST = SILVER / "reads_manifest.tsv"   # what is actually on scratch
HOST_CALLS = SILVER / "host_calls.tsv"           # 04_classify: the host called FROM READS
CLASSIFIED = SILVER / "_classified.tsv"          # 05_detect: one row per candidate detection

# ── gold: small, generated, tracked ─────────────────────────────────────────────────────────
DETECTIONS = GOLD / "detections.tsv"             # the presence calls this pipeline turns on
COINFECTION = GOLD / "coinfection.tsv"           # per-sample: declared, undeclared, cryptic
HOST_SUMMARY = GOLD / "host_summary.tsv"         # per-host rates, the headline numbers
# The saturation QC record: how much k-mer space each species' assemblies actually add, and
# what depth costs in congener misassignment. TRACKED because it is the EVIDENCE for the
# no-per-species-cap decision, which a reviewer may reasonably challenge, and because a
# decision defended by a number nobody can see is a decision defended by assertion.
SATURATION_CURVE = GOLD / "saturation_curve.tsv"
SATURATION_SUMMARY = GOLD / "saturation_summary.tsv"
SATURATION_BIAS = GOLD / "saturation_depth_bias.tsv"
SKETCH_INDEX = GOLD / "saturation_sketch_index.tsv"

# ── outside the repo ────────────────────────────────────────────────────────────────────────
# holobase supplies the LABELS (pathogen | resident | host) that 05_detect applies AFTER
# detection. It replaces PHI-base, which was applied BEFORE it as a filter over 632 fungal taxa
# and would silently skip every resident db_v4 adds. REFDATA comes from the root _paths.py so
# there is one definition of where the shared store lives.
#   taxon.ncbi_taxid is the real key; taxon.taxon_id is an internal surrogate that joins to 0%
#   without erroring. taxon.kingdom spells plants `Viridiplantae`, not `Plants`.
#   assembly.is_annotated is 1 on all 36,709 rows and cannot be used to ask about annotation.
HOLOBASE = REFDATA / "holobase.db"

# Acacia (Pawsey S3). The CDS downloaded for db_v3 are backed up here, and pulling from Acacia
# is free and internal where re-pulling from NCBI is neither, so 02_build always reconciles
# against this prefix BEFORE it calls `datasets`.
#   * **Two buckets exist and only one is right.** `pawsey1168-llenzo-kraken-db` sits on the
#     personal 100 GiB allocation, was filled on 2026-09-01 and holds 1,019 accessions;
#     `pawsey1168-kraken-db` sits on the project allocation and holds 1,994, a superset of
#     both the personal bucket and live scratch. Use the PROJECT bucket.
#   * **The profile name differs by machine**, which is why it is an argument and not a
#     constant: Setonix authenticates as `pawsey1168`, and the local `acacia` profile cannot
#     read the project bucket at all. It fails as `argument of type 'NoneType' is not a
#     container` — aws-cli 2.13 crashing while formatting a Ceph denial whose <Message> is
#     empty, not a network or SSL fault.
#   * **The --endpoint-url is mandatory.** Omitting it gives a misleading InvalidAccessKeyId.
#   * On Setonix the CLI is a Singularity wrapper that overrides --home, so it reads
#     $MYSOFTWARE/setonix/<ver>/.shpc_home/.aws, NOT ~/.aws. Credentials in the real home are
#     invisible to it and the symptom is "Unable to locate credentials".
ACACIA_ENDPOINT = "https://projects.pawsey.org.au"
ACACIA_BUCKET = "pawsey1168-kraken-db"
ACACIA_CDS_PREFIX = "db/search/cds/pathogen"
