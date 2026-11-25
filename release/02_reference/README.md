# 02_reference

Decides which organisms the classifier's database contains, and selects the assemblies that
can supply CDS for them. This file states the scope rule and the evidence for it.

The cohort is **1,560 runs over 31 BioProjects** (`release/01_cohort`), of which 1,450 are
*Triticum aestivum* and 1,379 of those are field. Wheat is the analysed host; maize, sorghum
and rice are out of the chapter as of 2026-11-17.

The other 110 runs are retained deliberately and are two different things. **84 are field runs
on a related host from a cohort BioProject** (*Hordeum vulgare* 57, *Triticosecale* 10,
*T. turgidum* 8, *Triticum* sp. 8, *Secale cereale* 5, *H. jubatum* 4, *T. dicoccon* 1,
perennial grass 1): same deposit, same protocol, a different host, which is what makes them
usable as a host contrast. **86 are lab**, 71 of them *T. aestivum*, and these are the axenic
and inoculated controls that anchor the detection floors. Neither set is noise.

## The rule

One scope decision and two labels, no routes.

    SCOPE   every fungus USDA ever recorded on a plant, read at GENUS, with an annotated
            assembly (data/inputs/scope_exceptions.tsv adds three genera USDA misses).
    LABEL 1 recorded on wheat?        USDA, host = Triticum              -> on_wheat
    LABEL 2 a known pathogen?         PHI-base, a NAMED DISEASE on any plant
                                                                 -> is_plant_pathogen
            a known wheat pathogen?   PHI-base, a named disease on Triticum,
                                      or the curated supplement  -> is_wheat_pathogen

The whole database in one sentence: every fungus found on a plant, with a note of when it was
recorded on wheat and when it was shown pathogenic on wheat. A detection then reads off those
columns: a species with none of them is a resident, one with `is_wheat_pathogen` is an expected
pathogen, and one that is `is_pathogen` but not `on_wheat` is a candidate unreported interaction.

**Scope is ANY plant, never wheat.** The study exists to find unreported host associations, and
a reference absent from the database is invisible and unrecoverable: its reads are simply not
classified. Contrary to the common assumption, they do not transfer to the nearest relative
present; measured on held-out congeneric *Puccinia* at k=31, ~96% of reads from an absent
species match no reference k-mer at all, and the few that do carry only a handful of conserved
k-mers (the pile-up signature, not a detection). So absence costs sensitivity, not false
attribution. Either way `on_wheat` and the PHI-base labels are columns, never filters:
filtering the database to what is already known on wheat would make the finding impossible to
make.

**Genus, not species.** Read at species rank the master loses 903 of 5,973 fungal detections
on field wheat (15.1%), led by *Dioszegia hungarica* and by unnamed strain designations USDA
cannot hold for want of a binomial (*Alternaria* sp. BMP 0032, *Cladosporium* sp. NU13). None
is pathogen-labelled, so the genus reading protects the resident layer, not the pathogen layer.

**Why PHI-base and not FUNGuild/FungalTraits for the pathogen label.** FUNGuild and FungalTraits
assign "plant pathogen" at GENUS by name substring, so they inherit every homonym and label
every species in a genus from one member: 37 of 39 pathogen labels db_v4 carried were genus
priors landing on chytrids, pollen parasites and the industrial yeast *Komagataella phaffii*.
PHI-base is curated per species from pathogenicity experiments and ships the NCBI taxid of both
pathogen and host, so `is_pathogen` and `is_wheat_pathogen` are a join, not a guess. The cost is
recall (PHI-base is smaller), but the label it gives is one you can defend species by species.

## The pipeline

    inputs/scope_rule.tsv        the criterion, and a pin on every dataset it is applied to
    inputs/scope_exceptions.tsv  admitted outside the master, and what was rejected

    fetch_species.py     -> raw/raw_usda_fungus_host.csv   827,645 records, md5-pinned
                            raw/raw_taxdmp.zip             NCBI Taxonomy, dated archive
    fetch_pathogens.py   -> raw/raw_phibase.csv            PHI-base v4-19, the pathogen label
    screen_species.py    -> curation/screen_species.tsv     99,059 names, scope + on_wheat
    screen_merges.py     -> curation/screen_merges.tsv       22 provisional NCBI names
    screen_pathogens.py  -> curation/screen_pathogens.tsv   150 PHI-base fungal pathogens
    fetch_assemblies.py  -> raw/raw_assembly_summary.txt    26,826 GenBank fungal assemblies
    screen_assemblies.py -> curation/screen_assemblies.tsv  21,413 in-scope assemblies
    fetch_hosts.py       -> raw/raw_ensembl_plants_species.txt  267 genomes, release-63
    screen_hosts.py      -> curation/screen_hosts.tsv        11 cohort hosts, 3 kept
    fetch_host_cds.py    -> raw/host_cds/<species>/*.cds.all.fa.gz   21 CDS sets, 1.06 GB
                            raw/host_cds_manifest.tsv        resolved file, bytes, md5
    plan_build.py        -> curation/build_plan.tsv          5,323 selected assemblies,
                                                             tiered: db0, db1, host
    fetch_cds.py         -> raw/cds/<accession>/cds_from_genomic.fna   5,272 CDS, 98 GB
                            raw/cds_fetch_manifest.tsv       accession -> status, bytes
    build_db.py          -> data/db/<name>/                  the Kraken2 index
                            outputs/db_inspect_<name>.tsv    minimizers per taxon
    make_supp_tables.py  -> outputs/supp_data_s3_reference.tsv  5,846 rows: 5,843 fungal
                            taxa + the 3 hosts, ONE ROW EACH

    raw/raw_sources_manifest.tsv   provenance for all four sources, TRACKED

Each fetch pairs with its screen: scope (species), pathogenicity (pathogens), buildability
(assemblies), the sink (hosts). A fetch gets pinned bytes; a screen decides. `screen_merges.py` sits between
them and reconciles the two naming systems; see "Two names, one fungus".

The module ends in TWO artefacts, and the names say which is which. `build_db.py` builds the
index; `make_supp_tables.py` writes the published table, and is named for what it does rather
than for the thing it describes (it was `build_reference.py`, which reads as the database
builder and was taken for one). `plan_build.py` decides and `build_db.py` acts, the same split
as a screen and a fetch, and they stay apart for a reason: `build_plan.tsv` must be IDENTICAL
across the arms of the ladder for a tier subtraction to isolate one effect, and a plan
regenerated inside a six-hour build can differ between arms with nothing saying so. That file
therefore lives in `curation/`, with the other decisions, and not in `outputs/`.

**`outputs/` is the deliverable and the only thing the figures read.** `make_supp_tables.py`
is the last step: it joins the four curation tables and the taxdump into one published row
per taxon, `supp_data_s3_reference.tsv`, and both panels are drawn from that and nothing
else. This follows `01_cohort`, whose figures read `supp_data_s2_cohort.tsv` rather than
reaching back into `curation/`. The curation tables remain the audit trail; a reviewer
holding only S3 can regenerate both panels.

The table is self-sufficient on purpose. It carries the six parent ranks per taxon, so panel
B's topology is rebuilt from those columns and nothing under `figures/` opens the taxdump.
Its columns are the two figure axes (`usda_row`, `phibase_col`), the five-way `category`
that is panel B's tip colour, the counts (`n_assemblies`, `n_annotated`, `annotation_ratio`),
the labels with their provenance (`wheat_pathogen_source`, `plant_pathogen_evidence`), the
`route` into scope, `merged_from`, `n_predicted` for the gene models predicted here rather
than downloaded, and `in_tree` marking the rows panel B draws.

**The pathogen label is a SEPARATE table on purpose.** It never enters the Kraken2 database,
which is pure taxonomy; it is an annotation joined to a detection afterwards. Keeping it out of
`screen_species.tsv` also gives one clean row per pathogen species rather than the label
inheriting across every USDA name variant.

Four public sources and nothing else (USDA, NCBI Taxonomy, NCBI assemblies, PHI-base). Nothing
is resolved against a private database: holobase was used during design and was **wrong twice**,
once in each direction, so the shipped pipeline does not read it. See "What holobase got wrong".

## The funnel

    99,059  distinct fungus names in USDA
    21,071  resolve to a fungal NCBI taxid          (28,858 names; 7,787 are synonyms
                                                     of a name already in the list)
     4,135  genera, which is what the master reads
     5,843  taxa inside those genera with any assembly
     2,314  of them carry a gene model

69.5% of USDA names have no NCBI node at all. That is a property of NCBI, not a matching
failure: USDA is a herbarium records database reaching back to 1960s literature and most of
those fungi have never been sequenced. A name with no node cannot have an assembly either,
because assemblies are keyed to taxids, so the 69.5% could never have contributed CDS.

**Synonymy is load-bearing and comes only from NCBI.** USDA ships nine columns and none is an
accepted-name field, so it records whatever name each publication used across six decades.
9,285 of the 28,858 resolutions arrive through a name NCBI files as a synonym, and **2,902
taxa are reachable by no other route**: drop synonymy and you lose 14% of the master outright,
not merely 9,285 spellings. `authority` and `type material` are excluded as name classes,
carrying the author string and a strain designation respectively.

**The unit is a TAXON, not a described species, and the figures say `taxa` for that reason.**
1,373 of the 3,465 genus-route taxa are NCBI `sp.` nodes: unidentified isolates with a strain
code and no binomial (*Alternaria* sp. BMP 0032). They are legitimate reference entries,
because an unidentified *Alternaria* gives Kraken2 as real a competitor as a named one, but
they are not species in any taxonomic sense, and 5,843 printed as "species" invites a
diversity reading the number cannot carry. The observation side does not have the problem:
it is 99.4% named binomials, because USDA only indexes formal names.

## Two names, one fungus

USDA and NCBI can hold the same fungus under two names, and then the reference counts it
twice and builds it once. *Cercospora flagellaris* is the case that exposed it: USDA has the
binomial with 33 host records and NO genome, while every genome is filed under the
provisional node *Cercospora* cf. *flagellaris*, which USDA never names. The name arrived as
an observation with nothing to build and the genome on the genus route, and the two
never met.

`screen_merges.py` folds the provisional node onto the binomial. The criterion is a token
match, never a substring: the NCBI name must parse as `<Genus> cf. <epithet>` with an
optional trailing strain code, and `<Genus> <epithet>` must be a `usda_name` in the master.
It writes `curation/screen_merges.tsv`, one row per affected node with its action, and that
table is the whole record: delete a row and the next run un-merges that taxon.

    merge     16 taxa,  61 assemblies   the binomial has its own NCBI node; move the
                                        assemblies onto it. Four of the targets had no
                                        genome at all before the merge.
    reroute    1 taxon,   2 assemblies  USDA names it, NCBI has no node for the binomial,
                                        so the provisional node is promoted in place
                                        (*Cercospora* cf. *sigesbeckiae*)
    keep       5 taxa,    8 assemblies  `aff.`, NOT merged

**`aff.` is deliberately not merged.** `cf.` ("compare with") says the isolate is probably
that species; `aff.` ("has affinity to") says it is probably a DIFFERENT, undescribed species
related to it. Merging an `aff.` node would assert the identity its own name denies, and an
undescribed relative is exactly what the genus route is for. The five are written out with
`action = keep` so the decision is visible rather than silent.

The merge also repaired a resolution fault worth recording. NCBI lists `Haematonectria sp.`
and `Neocosmospora sp.` as *equivalent names* of the single node *Fusarium* cf. *solani*, and
USDA carries both of those placeholders, so 24 assemblies were being counted as an
observation of a named species on the strength of two genus-level placeholders pinned to one
arbitrary isolate. They now sit under *Fusarium solani*, which is what USDA's 1,225 records
actually describe.

## Why one list and not five

The research tree reached its database through five routes: `target` (FUNGuild + PHI-base),
`congener`, `host`, `plant_associate` (USDA) and `host_associate` (literature residents).
Measured against db_v4's calls on the **1,135 classified field-wheat runs**:

| | detections |
|---|---:|
| species is in USDA | 5,053 of 5,973 fungal (84.6%) |
| species is not, but its GENUS is | 903 (15.1%) |
| **outside USDA entirely** | **17 (0.3%)** |

The 17 are *Malassezia* and nothing else: *restricta* 15, *arunalokei* 1, *globosa* 1.

The other four routes were then checked for what they uniquely contribute. **FungalTraits**
is fully redundant: 1,069 genus-keyed plant-pathogen genera, no species, and FUNGuild already
holds 1,068 of them, the exception being *Heterospora*. **PHI-base** adds no genera and 2
buildable species, with USDA a strict superset at genus rank. **FUNGuild** adds 13 species
USDA has no plant record for, of which only *Zygosaccharomyces bailii* and *Tilletiaria
anomala* have an annotated assembly; both were built into db_v4, both received reads across
hundreds of runs, and both produced zero detections, every one of their 845 rows reading
`below absolute floor`. **Residents** contributes the exceptions file. **Congener** was never
a scope list: it expands the genus of a known pathogen so the LCA has somewhere to park a
read, and a genus-level master does that for every USDA genus rather than only pathogen
genera.

These four comparisons were made in the research tree against holobase and have not been
re-derived from the public sources. The conclusions they support are each corroborated by the
detection evidence above, which has been.

### The genus tier is not a weak tier, and is no longer a row

**Re-derived on db2 and the 1,379 field-wheat runs, 2026-11-23**, independently of the db_v4
measurement above and agreeing with it to a tenth of a percent:

| | detections | |
|---|---:|---:|
| `route = usda_species` | 5,166 | 84.6% |
| `route = usda_genus` | 918 | 15.0% |
| `route = exception` | 21 | 0.3% |

Same leader as before, *Dioszegia hungarica* at 353 of 1,379 runs, and the exceptions are
still *Malassezia* and nothing else (*restricta* 19, *globosa* 1, *arunalokei* 1). Two
databases and two cohorts giving the same split is the strongest evidence here that reading
the scope at genus is doing real work rather than padding the reference.

**`usda_row` used to collapse `route = usda_genus` to the word `congener`, and the word was
wrong.** In this project's vocabulary "congener" means a relative admitted so Kraken2's LCA
has somewhere to park a read, with no expectation that the organism is present. Applying it
to a tier carrying 15.0% of all field-wheat detections invited the reader, and invited one
of this repo's own figures, to treat those calls as artefacts and drop them. **One did: a
draft of `03_classify/figures/samples` removed *Dioszegia hungarica* on the strength of the
label, which is the single detection the genus reading exists to keep.**

**Settled 2026-11-24, Leon's call: the ROUTE IS NOT A ROW AT ALL.** `usda_row` has two
values, `wheat_resident` and `plant_resident`, and they say only whether USDA names that
species on wheat. Three attempts to put the route on the axis were tried and dropped in one
sitting, and the reason each failed is the same: how a taxon got into the reference is a
fact about the BUILD, and drawing it as a row makes the figure assert something about the
ORGANISM that the route cannot support.

| attempt | why it was dropped |
|---|---|
| `congener` row | the word means something else here, and it read as "artefact" |
| `genus_only` row | honest name, but still a row for a build fact |
| `wheat_genus` + `plant_genus` rows | split the 3,465 by whether USDA records the genus on wheat (1,130 / 2,335). It put 1,130 taxa in the wheat half of the figure that are not wheat-associated organisms |

The route survives where it belongs, in the `route` column of S3, which distinguishes
`usda_species`, `usda_genus` and `exception` per taxon. Any analysis that needs the tier can
read it there; the figure no longer draws it. `figures/matrix/caption.py` still derives and
states the genus-route counts in words, so nothing is lost to the reader.

**What the tier actually mixes, and what nothing here can separate.** USDA Fungus-Host
records fungi OBSERVED ON PLANTS, not pathogens, so an absent species-level record is weak
evidence about anything. Two different cases arrive at `usda_genus` and no column
distinguishes them:

- a genuine plant associate USDA has not recorded under that name, which is what
  *D. hungarica* and *Sporobolomyces ruberrimus* are, both common wheat phyllosphere yeasts;
- a species of a plant-associated genus that is not itself plant-associated.

**Counted, that mix splits 2,865 / 600**, and `figures/matrix/caption.py` derives both at
write time rather than quoting them:

| | taxa |
|---|---:|
| genus is represented in the reference by a USDA-recorded SIBLING species | 2,865 |
| genus qualified on USDA records that could not be built into a row | 600 |
| genus USDA never names at any rank | 0 |

Only the first group is a congener in the strict sense. The 600 are genera whose USDA-named
species exist but have no annotated assembly, so the record that qualified the genus never
became a row of its own; they are lichens, mushrooms and mycorrhizas (*Cladonia* with 7
resolved USDA names, *Glomus* 34, *Agaricus* 15). The third row being zero is a property of
the scope rule, and the caption raises rather than absorbs a non-zero count, so a future
scope leak fails loudly.

**Match USDA names through `screen_species.tsv`, never by splitting the raw string.** A raw
token match reports *Dioszegia* as absent from USDA, which is how this was first miscounted;
it is there as *Bullera crocea* -> *Dioszegia crocea*. Synonymy is load-bearing here exactly
as it is for the master, and only 4 USDA names in the whole master resolve to a genus node,
so "USDA named only the genus" is not a case that occurs at any scale.

No available source resolves it: FungalTraits is genus-keyed (1,069 genera, no species) and
so cannot speak to a species its genus already covers, and that is why it was found redundant
above. The honest position is that **the label describes USDA's coverage of a NAME and is
never evidence about the organism**; whether the organism is there is what the detection and
the off-host null answer, independently.

Panel B is unaffected by the rename: `figures/tree/` keys on `category`, which never carried
the word.

## The pathogen label, measured

PHI-base v4-19 holds **150 fungal pathogen taxa; 118 pathogenic on a plant, 26 on wheat**
(18 from PHI-base, 8 from the curated supplement). The wheat set is what a reviewer would
expect: the *Fusarium* head-blight complex (*graminearum*, *culmorum*, *avenaceum*,
*asiaticum*, *pseudograminearum*), all three rusts (*Puccinia graminis*, *striiformis*,
*triticina*), *Zymoseptoria tritici*, *Bipolaris sorokiniana*, *Gaeumannomyces graminis*,
*Pyrenophora tritici-repentis*, *Pyricularia oryzae* (wheat blast), *Parastagonospora
nodorum* and *Blumeria graminis*, plus the eight canonical wheat diseases PHI-base does not
curate (`inputs/wheat_pathogens_curated.tsv`: *Tilletia caries*, *T. controversa*, *Ustilago
tritici*, *Oculimacula yallundae*, *Microdochium nivale*, *Claviceps purpurea*, *Rhizoctonia
solani*, *Alternaria alternata*). `wheat_pathogen_source` tags every call's provenance.

**A plant host is not enough; the record must name a disease.** PHI-base scores an
interaction against whatever host the experiment used, and the plant is not always the host
being harmed. In a tripartite insect-fungus-plant assay, or an endophyte colonisation study,
the plant is present but is not the victim, and those rows carry `Disease_name = "no data
found"`. Requiring at least one NAMED disease separates them. Five fungi in this reference
were labelled plant pathogens on disease-less rows alone, and every one is a plant
ASSOCIATE rather than a pathogen:

| | what it is | why PHI-base has a plant row |
|---|---|---|
| *Epichloe festucae* | the textbook grass endophyte, a mutualist | colonisation of *Lolium*/*Festuca* |
| *Aspergillus ochraceus* | storage fungus, ochratoxin A | grain and pear assays |
| *Leptographium clavigerum* | beetle-vectored blue stain on pine | vector association |
| *Metarhizium acridum* | locust biocontrol agent | tripartite assay on rice |
| *Metarhizium robertsii* | insect pathogen and root endophyte | root colonisation of barley |

That set is coherent in a way that argues the condition is sound rather than lucky. It is
also the same error this file already indicts FUNGuild for, arriving by a different route:
the join is exact, but what it joins on is the experiment's host, not the disease.

The condition is applied to the wheat label too, and is a **no-op** there: all 18 PHI-base
wheat taxa already carry at least one named disease, so the wheat count does not move.
Nothing is dropped from the table either way. A taxon whose only plant rows lack a disease
keeps its row with `plant_pathogen_evidence = host_only`, recording why it was not labelled.

**The label is applied at SPECIES rank**, because PHI-base keys pathogens by species taxid.
`screen_pathogens.tsv` is one row per pathogen species, so *Fusarium oxysporum* appears once as a
wheat pathogen rather than once per forma specialis. A detection resolves only to species and is
joined to this table on the species taxid, which is the matching granularity.

**No known wheat pathogen falls outside the scope.** Zero of the 26 wheat pathogens are
absent from the USDA master. So the broad scope loses nothing the narrow label cares about,
which is the test the whole design has to pass.

**No taxon admitted at genus rank carries a pathogen label at all.** After the merge and the
disease condition, all 3,465 genus-route taxa sit in the non-pathogen column, none of them a
known plant or wheat pathogen. That is a structural property, not a coincidence. Since the
2026-11-24 row change it is no longer visible as an empty row, so the caption states it in
words instead. A fungus PHI-base
scores with a named disease is a studied fungus, and USDA records studied fungi, so the
pathogen annotation can never reach a taxon admitted only because its genus qualified. No
detection can be called a pathogen on the strength of genus membership.

## What the rule drops

Of db_v4's 1,547 built taxa the rule retains 1,477 and drops 70.

    usda_species     964   kept
    usda_genus       494   kept
    exception         19   kept
    outside           70   DROPPED

Those 70 absorbed 8,815,071 reads, **0.11%** of the 7.9 billion fungal reads classified on
field wheat. The reads do not disappear when the taxa do; they are reassigned to the nearest
surviving taxon. That is the one effect this cannot predict, and only a rebuild settles it.

No pathogen-labelled detection falls outside USDA. Zero of the 17 carry one.

A caveat worth keeping: this measures whether USDA covers what the UNION database found, not
whether a USDA-only database would have found the same. The union is the generous case, so it
is the strongest evidence available rather than a proof.

## What holobase got wrong

holobase is a private tool and is not in this pipeline, which is why. Both errors were found
by checking the public sources directly.

- ***Cryptococcus*** was reported absent from USDA and carried as an exception. It has **39
  USDA records, 18 of them on *Triticum aestivum***. Almost certainly the 5206-fungus /
  79213-scale-insect homonym its own seed file warns about. It is not an exception.
- **The `assembly` table cannot answer whether an assembly exists**, in either direction. It
  recorded zero for *Tilletiaria anomala*, which NCBI holds two annotated assemblies of, and
  three accessions for *Hansenula subpelliculosa*, which NCBI returns nothing for. Its
  `is_annotated` is 1 on all 36,709 rows, so it cannot answer annotation either.

Annotation now comes from `assembly_summary.txt`'s own `annotation_provider`,
`annotation_date` and `protein_coding_gene_count` columns. GenBank only: RefSeq holds
re-annotated copies of assemblies already in that file, and taking both is how db_v4 came to
build 432 genomes twice under paired GCA/GCF accessions, every pair differing in CDS size.
That defect is retired by construction rather than by a dedup rule.

## The host sink: why there is wheat in a fungal database

**db0 and db1 were built with no plant CDS at all** (`Fungi 100%` in `db_inspect.tsv`), and
that was a regression rather than a decision. The scope rule above says "every fungus USDA
ever recorded on a plant", a wheat genome satisfies no clause of it, and when the module was
written around that sentence the host sink silently fell out. db_v3 had held 38 plant species,
568 M minimizers, 11.8% of the database.

Two independent things break without it.

**The sink.** Without host CDS a plant read is absorbed by the nearest fungal CDS and reported
as a fungal hit. These are plant RNA-seq libraries, so the great majority of reads are wheat,
and that is a false-positive mechanism aimed straight at the claim this chapter makes.

**The host call is a result, not a diagnostic.** `04_assign` ranks reads across hosts to decide
which plant a library came from, and agrees with the LLM-curated host 91.2% of the time. With
no host CDS there is no host call to make.

Neither is contradicted by the measured finding that host-read REMOVAL is unnecessary (PST read
counts ratio 1.00x raw against HISAT2-split across five ground-truth runs). That tested
deleting host reads BEFORE classification, which changes the percentage and not what Kraken2
finds. Putting host CDS IN the database is the opposite operation: it gives a plant read a
correct home instead of a fungal one.

### Three hosts, and the exclusions matter more than the inclusions

`screen_hosts.py` is a SCREEN OVER THE COHORT, not an authored list. The research tree's
`seeds/hosts.tsv` was hand-written for a four-cereal cohort and still names *Oryza*, *Zea* and
*Sorghum*, none of which has a run left.

    analysed   Triticum aestivum     1,450 runs   19 Ensembl cultivars
    contrast   Hordeum vulgare          57 runs   reference only (MorexV3)
    outgroup   Arabidopsis thaliana      declared  reference only (TAIR10)

**ENSEMBL PLANTS, NOT NCBI, and the difference is decisive.** NCBI GenBank holds the wheat
pangenome as ASSEMBLIES and not as gene models: of its twelve *T. aestivum* entries carrying
an `annotation_provider`, **eight have `protein_coding_gene_count = 0`** -- the NIAB cultivars
(Alchemy, Xi19, Robigus, Soissons, Brompton, Rialto, Claire) are assembled and not annotated
-- leaving only Chinese Spring, annotated twice. So "all available wheat CDS" from NCBI is two
annotations of one cultivar plus three fragments, and `annotation_provider` is not a test for
"has genes": selecting on it would queue eight 15 Gb downloads that each yield an empty CDS
file. Ensembl Plants publishes **19 *T. aestivum* cultivars with CDS**, which is the breadth a
sink actually wants.

Ensembl is also the better provenance story. It ships immutable numbered releases, so plants
release-63 is a stable artefact a reviewer re-fetches byte-identically and `expected_md5` is
pinned; the fungal `assembly_summary.txt` is rebuilt daily and cannot be pinned at all. Note
the two numbering systems: Ensembl Genomes numbers plants releases separately from Ensembl
core, so `summary.txt` calls release-63 "Release 116" because it loads into Ensembl v116
databases, and `/pub/release-116/` does not exist.

**BREADTH SCALES WITH COHORT REPRESENTATION.** The analysed host takes every cultivar, because
that is where the reads are and a cultivar-specific transcript missing from the sink leaks onto
a fungus. A contrast or outgroup host takes its reference only: barley is 57 runs and
Arabidopsis is pure QC, so the 77 barley cultivars Ensembl carries would be weight with no
question behind it.

Taking all 19 wheat cultivars is safe in a way one-species-per-genus is not, and for the same
reason: **they all share taxid 4565**, so they consolidate on ONE node and no genus-level LCA
bleed can occur, and Kraken2 stores each distinct minimizer once so the redundancy adds no
weight to the taxon.

**ONE SPECIES PER GENUS, and it is the load-bearing clause.** It drops hosts the cohort really
contains: *T. turgidum* (8 runs), *Triticum* sp. (8), *T. dicoccon* (1) and *H. jubatum* (4).
Adding them would be actively harmful. db_v3 held both hexaploid *T. aestivum* and tetraploid
*T. turgidum*, and 27.7% of the *Triticum* clade's 51.8 M minimizers migrated to the GENUS
node, so a correct wheat call scored low coverage; db_v4 fixed that by consolidating.
Re-adding a congeneric host would re-create a defect already paid for. Those runs call as
*T. aestivum*, which is taxonomically honest and costs nothing. Note the floor is a property of
the GENUS (*Hordeum* is 61 runs across *vulgare* and *jubatum*) while a database entry is a
SPECIES, so the floor qualifies a genus and one representative is picked from it.

### Two traps that live in this screen

**Ensembl files barley under 112509 (*H. vulgare* subsp. *vulgare*), not the cohort's 4513.**
Hosts cannot be matched on taxid equality, so `screen_hosts.py` walks `nodes.dmp` and accepts
any Ensembl genome at or below the cohort's taxon, then tags the CDS with the COHORT's species
taxid so a host call resolves at the rank the cohort records and the chapter reports. Same
reasoning as using 1047171 for *Z. tritici* rather than the IPO323 strain node 336722.

**The CDS filename cannot be constructed from the species table.** Its `assembly` column and
the FASTA's assembly token disagree -- Landmark is `PGSBv2.0_Landmark` in the table and
`PGSBv2.1` in the file, Weebil is `EI_WeebilV1` against `WeebilV1` -- so building URLs 404s on
2 of 19 cultivars, silently, because one missing cultivar of nineteen looks like nothing.
`fetch_host_cds.py` lists each `cds/` directory instead and records the resolved filename in
`host_cds_manifest.tsv`. It takes `.cds.all.fa.gz`, the gene-build CDS, never
`.cds.abinitio.fa.gz`, which is unreviewed prediction Ensembl excludes from its gene set.

*Triticosecale* (10 runs) and *Secale cereale* (5) fall below a 1% floor on the cohort, and
**the sink survives their exclusion** because Triticeae CDS are near-identical: a rye or
triticale read lands on wheat or on a Triticeae ancestor, and either way not on a fungus,
which is the whole job of the sink.

**The outgroup is declared, not derived, and that is deliberate.** Host assignment ranks rather
than thresholds, so with only cohort hosts present every library is called as one of them by
construction and a mislabelled or contaminated library cannot be seen. A distant dicot gives
that signal somewhere to land for ~10 M minimizers. It is how *Fragaria*'s absence was found:
30 strawberry runs ranking onto *Populus*, *Arachis* and *Glycine*.

Every cohort host gets a row in `screen_hosts.tsv` with its `decision` and `reason`, kept or
dropped, following `screen_merges.tsv`. A host crossing the run floor changes its row rather
than appearing or vanishing without trace.

### What this does to S3 and to the figures

**S3 is now the DATABASE CONTENTS, not the fungal scope.** Host rows carry `route=host`,
`kingdom=Viridiplantae` and empty fungal columns, because no clause of the scope rule and no
PHI-base label means anything for a wheat genome. They are published here rather than only in
`build_plan.tsv` because that is not a `supp_data_*` table, so hosts would otherwise be absent
from the supplement entirely, describing a database that contains wheat as though it did not.

**Anything counting the funnel or drawing a panel filters `kingdom == "Fungi"`.** Both figure
scripts and both captions do; `matrix_numbers.tsv` is byte-identical before and after the
hosts were added, which is the check that the filter works. Every host row also carries
`in_tree=no`, so panel B is unaffected either way.

**The host tier is in EVERY database and is never subtracted.** `build_plan.tsv` tags it
`host` rather than giving it a rung on the db0/db1/db2 ladder: a comparison where one arm has
the sink and the other does not measures the sink, not the tier.

## The gene-prediction gap, handed to research/predict

3,529 in-scope taxa have no gene model. They are not equal:

    genus already has annotated assemblies   2,597   predicting here mostly adds competitors
    genus has NO annotation at all             932   the gap that costs calls

**The 932 was the wrong target and is superseded.** The unit is the GENUS, since one
prediction gives every species in it a k-mer neighbour, and 932 species is only 437 genera --
led by *Leccinum*, *Inocybe* and lichenised taxa that do not occur in a wheat leaf. The live
worklist is **135 genomes** (`research/predict/worklist.py`): 118 silent reportables in the
wheat row or a pathogen column with no gene model, 3 on-wheat dark genera, 14 congener
competitors for callable wheat-relevant taxa that have none. All 135 are predicted and
validated; this section is rewritten when the prediction arm lands here as supplementary
table S4.

Predicted CDS is a **competitor only** and must never support a positive call about its own
species, because homology projection biases toward genes conserved with the guide. Every
predicted species therefore dilutes species-level resolution for the species that can be
reported, and buys nothing reportable itself. That argues for the smallest set that closes
real misattribution gaps. The worked example is *Puccinia hordei*, absent from db_v3 and read
as *P. triticina*.

All 932 are in genera with no annotation, so none has an annotated congener to guide miniprot.
That makes this set ab initio's case.

## Decisions not yet taken

- **The host side is not screened.** The criterion says "recorded on at least one plant host",
  and testing it needs the 64,961 host names resolved too. What is implemented is "recorded in
  USDA at all". The dataset is a plant-pathology compilation so the two nearly coincide, but
  they are not the same statement: *Clohesyomyces*' single record gives the host as `"wood"`.
- **`assembly_summary.txt` is unpinned.** It is rebuilt daily and NCBI publishes no archive,
  so there is nothing stable to pin to. The md5 and the server's `Last-Modified` go into the
  manifest instead, and setting `expected_md5` in `scope_rule.tsv` turns enforcement on.
- **The declared-pathogen set has not been ported.** `target_pathogens` is empty for 53 of the
  63 studies in `release/01_cohort`, so the test that no declared pathogen falls outside the
  master cannot yet be run against the live cohort.

## The figure

Two panels, each owning a directory under `figures/`: its scripts, its data, its caption and
its output live together, following `02_classify/report/figures/<name>/`. Every number is
read from the curation tables at draw time; none is transcribed.

    figures/matrix/   plot_matrix.py  -> matrix.png + matrix_numbers.tsv
                      caption.py      -> matrix_caption.txt
    figures/tree/     prep_tree.py    -> tree_{tree.nwk,tips,clades,params}
                      tree.R          -> tree.png
                      tree.R          -> tree_arcs.tsv (the shape per class)
                      tree_legend.R   -> tree_legend.png
                      caption.py      -> tree_caption.txt

Both read `data/outputs/supp_data_s3_reference.tsv` and nothing else. The two caption
scripts additionally read `curation/screen_merges.tsv` and `screen_pathogens.tsv`, for the
two facts S3 cannot carry: the `aff.` nodes deliberately not merged, and a labelled taxon
with no in-scope assembly. Neither is a row in the reference.

Run, after the screens:

    python scripts/make_supp_tables.py
    python figures/matrix/plot_matrix.py && python figures/matrix/caption.py
    python figures/tree/prep_tree.py && Rscript figures/tree/tree.R \
      && Rscript figures/tree/tree_legend.R && python figures/tree/caption.py

**Panel A is a matrix because the two facts are independent.** A taxon carries where USDA
recorded it and what PHI-base scored it as, and those do not nest. Drawn as a Venn, the
awkward combination (a plant pathogen USDA records on wheat but PHI-base never scores on
wheat) needs an overlap lens with no name, and every box ends up labelled with the inclusive
set name while holding the exclusive count. Crossed axes draw as a matrix, where that
combination is simply a cell, and the counts are disjoint and sum to the cohort.

**One cell is empty and is not drawn.** The gap is the result: no wheat pathogen in scope is
off-wheat in USDA. The captions say
which gap means what, because an absent box on its own reads as a layout choice.

**Five colours, one per occupied cell, shared by both panels.** Hue is the PHI-base column,
saturation the USDA row (`scripts/_palette.py` `MATRIX`). Panel B's five tip categories ARE
those five cells, which is why `prep_tree.py` emits five and not four: a plant pathogen USDA
records on wheat is a different cell from one it does not. The hexes are repeated in
`tree.R` and `tree_legend.R` and the three must be kept in step.

**A PREDICTED GENE MODEL PUTS A TAXON IN THE TREE, as of 2026-11-25.** `in_tree` was
"observation with an NCBI annotation", which excluded all 135 taxa `research/predict` ran
on: every one of them had no annotation, which is exactly what put it on the worklist, so
the panel could not show the arm built to close that gap. The rule is now "observation with
a gene model, downloaded OR predicted" and the tree is 1,300 tips rather than 1,172. 128 of
the 135 are observations and became tips; the other 7 arrived on the genus route and stay
off, on the same rule as every other genus-route taxon.

Two thin rings at the tip base carry the two facts and must not be conflated. The INNER ring
is NCBI's annotation ratio, the purple ramp, and still measures the public gap. The OUTER
ring is a black tick where the model is ours. They cannot both be bright for one taxon: a
predicted taxon has no NCBI annotation at all, so its inner ring is at the dark end and its
whole bar is predicted. `caption.py` asserts that rather than trusting it. A first attempt
put the ticks outside the bars, where most tips never reach and the marks read as dust.

**Sized to be placed, not cropped.** All three files render at their final print size at
300 dpi and tile an **A4 portrait** page (8.27 x 11.69 in): the fan on top at full width,
the matrix as a landscape band beneath it with the key beside it.

    width   0.30 + tree 7.50 + 0.30                                  =  8.10
    band    0.30 + matrix 3.85 + 0.20 + legend 3.67 + 0.30           =  8.32
    height  0.30 + 7.50 + 0.20 + max(matrix 2.17, legend 2.89) + 0.30 = 11.19

Re-measured 2026-11-25. The matrix lost a row and gained a statistic per pill (3.80 x 3.06
-> 3.85 x 2.17) and the key gained a row for the predicted-model ring, so the band is now
as tall as the KEY rather than as the matrix.

Stacking beats side by side here. A square fan placed next to the matrix wastes whatever
height the matrix does not use, and the fan is then limited by the width the matrix leaves
it. Stacked, the fan is limited only by the page, and the band fills the strip left over.
`tree.R` scales every stroke and type size by `OUT_IN / 15`, because the fan was designed at
15 in and shrinking that PNG into the layout would take 12 pt numerals down to about 5 pt.

**The tree carries no legend and no title.** `tree_legend.R` draws the key as its own file,
after `01_cohort`'s `cohort_network_legend`, so the fan can be placed and scaled without
dragging a legend block with it. It is two columns, bar colour and arc symbol, because one
column would be about 5 in tall against a 3.3 in band.

**The clade arcs are marked with a SYMBOL, not a class name.** Set radially, a 15-character
name reads at a different angle in every quadrant, half of it needs flipping, and at this
size it competes with the bars it annotates. A point shape is orientation-free. `tree.R`
writes the assignment to `tree_arcs.tsv` as it draws, and both the key and the caption read
it back, so neither can drift from the panel.

**Captions are generated, never hand-edited.** Each `caption.py` reads the tables and
interpolates every figure, after the drift recorded in `crypt-figure-caption-discipline`.
Paragraphs are wrapped at write time so a number changing width cannot leave a ragged line.

**Superseded.** `plot_partition.py` (the nested Euler), `plot_schematic.py` (the
source-banded flow), `plot_funnel`, `plot_collapse`, `plot_butterfly` and `plot_legend` are
in `research/superseded/02_reference_figs/` with their outputs. PDFs are no longer written.

## Layout

    data/inputs/     what we supply: the criterion, the pins, the exceptions
    data/raw/        what the three sources return. Large files gitignored, manifest tracked
    data/curation/   decisions, by hand or from a stated criterion (incl. screen_hosts.tsv)
    data/outputs/    THE DELIVERABLE: supp_data_s3_reference.tsv, what the figures read
    scripts/         the pipeline steps, plus _palette.py and _caption.py
    figures/matrix/  panel A: its script, its numbers, its caption, its png
    figures/tree/    panel B: same, plus the separate legend
