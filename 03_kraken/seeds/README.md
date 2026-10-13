# seeds/ — the single source of truth for what is in the CDS database

Before this directory existed, membership was implicit: it was whatever `kraken_search.py`
happened to select from PHI-base, and `01_search/data/` held eleven tables with no statement of
which one was authoritative. A taxon could be absent for three unrelated reasons — nobody
decided to include it, the selector dropped it on a quality gate, or no assembly exists — and
nothing distinguished them. That is how the resident-community gap went unnoticed through two
database builds.

This directory makes membership explicit. Nothing may be dropped from the build without a row
in `exclusions.tsv` saying why.

## The scope rule

> **Plant-associated, ANY host. Never scope by the cohort's own hosts.**

USDA Fungus-Host Distributions is host-**association** data. It is used here to ask *is this
fungus ever recorded on any plant*, which is a scope gate. It must **never** be used to ask *is
it recorded on wheat, maize, rice or sorghum*, which is a plausibility filter.

CLAUDE.md already rejects host-plausibility filtering of calls, because the observed
fungus-host edges record only what has already been seen and the study exists to find
unreported associations. Moving that same filter earlier, into reference selection, is **worse
than filtering calls**: a filter on calls can be relaxed afterwards, but a missing reference is
invisible and unrecoverable without a rebuild, because the reads go silently to the nearest
taxon that *is* present and the detection looks genuine on every axis. See
`../01_search/README.md`.

It is a one-line change away from being wrong. `build_seeds.py` filters hosts on
`kingdom = 'Viridiplantae'` and nothing narrower.

## Authored versus generated

| file | status | holds |
|---|---|---|
| `sources.tsv` | **authored** | each seed source, its role, and why it is here |
| `hosts.tsv` | **authored** | the 6 plant species, and why each one is in |
| `exclusions.tsv` | **authored** | explicitly OUT, with a reason per row |
| `residents.tsv` | **authored** | the abundance-driven genus list, with evidence per row |
| `fungi_usda.tsv` | *generated* | `build_seeds.py`, from holobase: WHICH fungi are plant-associated |
| `fungi_annotated.tsv` | *generated* | `build_seeds.py --annotate`, from NCBI: which have a GENE MODEL |

Authored files are irreplaceable: losing one costs re-reading papers and re-making decisions.
A script may never overwrite them. `build_seeds.py` writes only `fungi_usda.tsv`. This split is
the same rule `01_literature/_layout.py` states, for the same reason.

## What is in scope

**FUNGI ONLY** as of 2026-10-13, Leon's call. Oomycetes are out of scope: db_v3's 29 oomycete
species (185 M minimizers, 3.8% of the database) come out, and zero of them pass both
`05_filter` floors on any cereal host, so nothing measurable is lost. It also avoids inheriting
a holobase gap — holobase holds no oomycete assemblies at all, so USDA's 738 plant-associated
oomycete taxa would each have needed separate NCBI resolution. CLAUDE.md's Scope section still
says "fungi and oomycetes" and must be updated.

**Fungi, primary — USDA, 1,486 taxa** with a scaffold-plus assembly, of which 432 are already
in db_v3, so **1,054 are new**.

**But only 179 of those 1,054 can be added today**, holding **369 unique annotated assemblies**
(`fungi_annotated.tsv`, from `--annotate`). **849 — 80% — have a scaffold-plus assembly but no
gene model at all**, and 26 have nothing. Since the database is CDS-only, annotation is not a
detail: it is the gate, and it closes on four fifths of what USDA otherwise offers. Miniprot is
therefore not a side quest but the main lever on breadth.

PHI-base was db_v3's basis, via `01_stat`, and is dropped as the scope source: `01_stat` is out
of the pipeline, and `pathogenic` resolves to only **340 fungi** against USDA's 20,204 on
plants. `recorded_on` is the right relation precisely because it reaches the saprotrophs,
endophytes and epiphytes that `pathogenic` structurally cannot. PHI-base remains *reference
data* — `phibase_db.json` still labels which taxa are declared pathogens, which `05_filter`
uses — just not the scope.

**Fungi, secondary — `residents.tsv`, 30 genera.** USDA is a herbarium records database, so a
leaf yeast is barely recorded as "fungus on host X": *Papiliotrema*, *Dioszegia*,
*Cystofilobasidium* and *Udeniomyces* appear at 1-2 records each, and *Candida* and *Mrakia* not
at all. But the abundance literature says **seven of the top 20 wheat genera are yeasts**, and
*Papiliotrema* is rice's single most abundant genus at 14.5%. USDA alone would reach
*Cladosporium* and *Nigrospora* and miss that whole layer. Two sources, because neither is
sufficient.

**Hosts — 6 species.** The four cohort crops plus *Hordeum vulgare* (monocot outgroup) and
*Arabidopsis thaliana* (dicot outgroup). db_v3 held 38 plant species from an era when the
cohort spanned 36 hosts; 32 are dropped, which is 413 M minimizers or 8.6% of db_v3. The
outgroups are kept because `04_assign/host` *ranks* reads rather than thresholding, so with
only cohort hosts present every library is called as one of them and a mislabelled non-cereal
library becomes undetectable. That ranking is how *Fragaria*'s absence was found.

## Sizing: availability is not the build

Two things in holobase cannot answer the size question, and both have already produced
confidently wrong numbers:

- **`assembly.is_annotated` is 1 on all 36,709 rows.** It cannot filter for annotation. The
  database is CDS-only (`kraken_build.py` skips any assembly without a gene model), so
  annotation is the binding constraint and holobase cannot see it.
- **`taxon.taxon_id` is an internal surrogate, not an NCBI taxid** — the schema says so, and
  the real key is `ncbi_taxid`. Joining db_v3's taxids against `taxon_id` returns a confident
  **0% overlap** rather than an error.

So the only way to size a build is to measure it. Two measurements, both of which cut the
estimate hard:

- `--annotate` against NCBI: of 1,054 new taxa, **179 are addable**, not 1,054.
- `--select-only` on the resident set: availability said 508 assemblies, the selector took
  **161**, because the depth floor acts as a cap for any genus STAT never detected.

Every estimate made from availability counts alone has been too high by roughly 3x.

### The build is close to size-neutral

Once annotation is measured rather than assumed, db_v4 barely grows:

| | assemblies | fungal species | minimizers |
|---|---:|---:|---:|
| db_v3 | 2,063 | 877 | 4.80 B |
| − 32 hosts | | | −413 M |
| − 29 oomycete species | | −29 | −185 M |
| + USDA addable | +369 | +179 | |
| + residents | +161 | +80 | |
| **db_v4** | **~2,590** | **~1,107** | **~5.2 B** |

About +26% assemblies and +26% fungal species, for roughly +8% minimizers after the cuts. The
earlier 6.9-8.9 B projection assumed all 1,054 new taxa were addable; 80% are not. A ~5.2 B
index is ~28 GB, comfortably inside a Setonix standard node's 230 GB, and the CDS download is
~530 assemblies rather than ~3,000.

**So size was never the constraint, and nor is download. Annotation is.** The 849 miniprot taxa
are where the remaining breadth lives.

## Regenerating, and the order to work in

```
python 03_kraken/seeds/build_seeds.py --check     # report, write nothing
python 03_kraken/seeds/build_seeds.py             # -> fungi_usda.tsv
```

Then size the build before downloading anything:

```
python 03_kraken/01_search/kraken_search.py --skip-seeds \
    --extra-clades 03_kraken/seeds/residents.tsv \
    --sweep-rank genus --select-only --version residents --no-symlink
```

`--skip-seeds` makes a standalone run: no PHI-base seed query, no seed-derived sweep targets,
and `--retain` forced off. `--no-symlink` is **required** for any side run, because
`ref_candidates.tsv` is a symlink that `utilities/busco.py` and the build read through, so
repointing it at a partial table silently makes that partial table "the database".

Column 1 of `residents.tsv` may be a **taxid** rather than a name, and must be for a taxonomic
homonym: the `datasets` CLI refuses an ambiguous name, and `datasets_query()` used to return
`[]` on the non-zero exit, so the clade selected nothing with no warning. *Cryptococcus* is both
a fungus (5206) and a scale insect (79213); *Rhizophagus* is both a fungus (1129544) and a
beetle (295919). Between them that cost 42 assemblies and 17 species. `datasets_query()` now
warns.

## Open decision

`select_rank_sweep` requires `scaffold_plus(a) and has_annotation(a)`. For six genera those
never co-occur, because the only annotated assembly is **contig-level**: *Boeremia* (3),
*Dioszegia* (2), *Papiliotrema* (1), *Clohesyomyces* (1), *Leucosporidium* (1), *Phaeosphaeria*
(1). That excludes *Papiliotrema*, rice's most abundant genus, and *Dioszegia*, a wheat core
yeast. Relaxing the gate for residents recovers 9 assemblies across 6 species. Arguments both
ways are in `../01_search/README.md`. `SCAFFOLD_PLUS` in `build_seeds.py` would have to move in
step. Not taken: it gates the whole database and db_v3 was built under it.
