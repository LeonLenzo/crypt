# metadata — Study Design Enrichment and LLM Classification

## Rationale

A co-infection detection in isolation carries limited biological meaning without study design context. The central confound is that not all positive detections represent ecologically authentic co-infections: a run from a co-inoculation experiment, an abiotic stress trial with a pathogen treatment, or a controlled greenhouse inoculation may show secondary signal by design, not by incidental co-infection. Conversely, a field-collected sample showing secondary signal is far more likely to represent a genuine unreported co-infection.

This module enriches the STAT detections from `runs.tsv` with three layers of context: (1) BioProject and BioSample metadata (including ENA/DDBJ BioSamples, which return empty from NCBI efetch), (2) linked literature and full manuscript text, and (3) LLM-based study design classification. Together these allow stratification of detections by study intent, setting, tissue, and stated pathogen/host — the necessary prerequisite for any epidemiological interpretation.

## Methods

### BioProject and BioSample metadata + literature linkage (`meta_search.py`)

`02_literature/01_search/meta_search.py` resolves BioProject/BioSample metadata and literature identifiers in one atomic pass per BioProject, matching NCBI's actual resolution cascade: BioProject XML `<Publication>` field → PMC full-text search → PubMed title search — all NCBI-XML derived, including bare-DOI-no-PMID cases. Only if all three find nothing does it fall through to a Serper web search, then finally DOI extraction/page scrape/CrossRef/PMC-by-DOI. This replaces the old `ncbi_metadata.py` + `web_metadata.py` split.

BioSample XML is parsed for a harmonised field set (`geo_loc_name`, `tissue`, `collection_date`, `isolation_source`, `dev_stage`, `lat_lon`, `host`). Coverage is uneven — the SRA submission process does not mandate these fields, and a nontrivial fraction of populated fields are NCBI placeholder values (`missing`, `not applicable`, `not collected`) rather than real data — see Results below for filtered coverage numbers. `meta_search.py` also fetches ENA/DDBJ BioSamples (`SAME*`/`SAMD*` accessions, which return empty from NCBI efetch) via the EBI BioSamples API; its `bs_description` field (e.g. "RNA-Seq of a field sample of Wheat Yellow Rust") is a particularly strong signal fed directly into the LLM classification prompt.

### Full-text retrieval (`meta_text.py`)

`02_literature/02_text/meta_text.py` retrieves full manuscript text, in cascade order: PMC OA full-text XML (free, complete, no PDF-parsing artifacts, when `meta_search.py` already resolved a PMCID) → Unpaywall → PDF download → `pdfminer` extraction → manual PDF fallback (`--ingest-manual`, for hand-downloaded PDFs matched against the failed-DOI list). `--apply` writes `full_text` back into `bioprojects.json` in place, which gates entry into classification below — **only BioProjects with retrievable full text are classified at all**, avoiding the old pipeline's failure mode of guessing study design from a title alone.

### LLM study design classification (`meta_classify.py`)

`02_literature/03_classify/meta_classify.py` replaces the old keyword-based classifier (`filter_kw.py`) and single-pass LLM classifier (`llm_classify.py`) entirely — keyword classification is not used at all in the current pipeline; the old approach's ~51% misclassification rate for `host_study` (pathogen/disease language co-occurring with host biology language) made it unreliable as anything but a discarded baseline.

Six LLM calls per BioProject (`gpt-4o-mini`): five independent judgment dimensions — `stress` (biotic/abiotic/none), `study_setting` (field/greenhouse/growth_chamber/detached_leaf_assay/in_vitro/unclear), `tissue` (aerial/non-aerial/unclear), `coinfection_intent` (single_pathogen_focus/not_disease_focused/intentional_multi_pathogen), and `hostpath` (named pathogens + named hosts, each as a list with its own confidence) — each with its own confidence + rationale, plus one fact-extraction call for symptom status, exposure type, geographic location, library prep, host cultivar, and host resistance. `--focus {stress|setting|tissue|coinfection|hostpath|extract}` reruns a single dimension.

A separate per-BioSample host disambiguation pass (`--disambiguate-hosts`) resolves which specific host a BioSample belongs to when a BioProject's `named_hosts` has more than one value (e.g. a rust fungus's alternate-host life cycle spanning two plant genera) — `named_hosts` itself is extracted once per BioProject and would otherwise be duplicated identically across every BioSample in a multi-host study, which is structurally wrong for per-sample analysis. The disambiguation prompt uses per-BioSample metadata (`bs_host`, `bs_description`, `tissue`, `isolation_source`) to pick the best-matching candidate, with its own confidence + rationale, and its own cache keyed by BioSample (not BioProject) — including the exact candidate list each resolution was made against, so a later prompt improvement that reshapes the candidates correctly invalidates stale resolutions rather than silently keeping them.

Author-named pathogens/hosts are resolved to NCBI taxids deterministically in plain Python afterward (`_util.resolve_taxon_name()`), never asked of the LLM — an LLM recalling taxids from memory produces confident-looking wrong numbers.

### Provenance enrichment (`supp_provenance.py`, `cohort_provenance.py`, `geocode_localities.py`, `provenance_tracker.py`)

BioSample XML carries usable geography for 56% of samples and a collection date for 46%, which
is too thin for spatial or temporal analysis. Papers' supplementary tables often carry the
per-sample detail the SRA submission omitted, so a fourth layer recovers it.

`02_text/supp_provenance.py` mines per-sample provenance from supplementary spreadsheets.
Adding a supplement is a CONFIG entry (file, sheet, column names), not code; the header row is
found by content, so banner rows above it need no offset. Journals publish whichever accession
flavour they prefer (ERS/SRR/SAMEA), so each is resolved to a BioSample through ENA and cached.
ENA's comma-batched accession query silently returns nothing, so accessions are queried one at
a time.

`03_classify/cohort_provenance.py` merges three location sources in priority order: supplement >
BioSample `geo_loc_name` > `llm_geographic_location`. It also normalises country names
(`COUNTRY_ALIAS`, `NOT_A_COUNTRY`), which collapsed 61 raw strings to 53 real countries by
unifying UK/United Kingdom, mapping US states that had leaked into the country field, and
dropping "South America" as not a country.

`03_classify/geocode_localities.py` geocodes locality strings via Nominatim at 1 request/second,
cached. It parses the three naming conventions in play (NCBI `Country: region, locality`, most
specific last; LLM `locality, country`, most specific first; supplement bare name) and rejects
hits whose class is not a settlement (road, building, amenity) or whose country disagrees with
the record. `accept-language=en` is mandatory, or Nominatim returns local-language country
names and every country check fails.

`03_classify/provenance_tracker.py` is the worklist, ranking BioProjects by missing provenance
and classing them DONE / ok / TARGET / SETTING? / NO DOI. `SETTING?` means the title indicates
controlled or in-vitro work, so a single location is already correct and chasing its supplement
cannot help.

Three supplements are configured. The largest by far is Adams et al. (2021) *BMC Genomics*, the
rust expression browser, whose metadata table spans 12 studies and covers 898 of our BioSamples.
Together they raised location coverage from 81.4% to 86.0% and collection year from 57.3% to
69.7%; PRJEB31334 went from a single "Ethiopia" label to 75 distinct localities, PRJEB15280 from
none to 69.

Two cautions, both learned the hard way. **Not every supplement describes the sequenced
samples**: PRJNA1217477's PNAS supplement tabulates the origins of the inoculum *isolates*
(Californian vineyards), not the glasshouse plants that were sequenced, so merging it would have
stamped vineyard coordinates onto *Arabidopsis*. Check the table keys against a sample accession
before configuring it. And a collector's surname in a Location field geocodes happily to a
street address.

The gain has plateaued. Almost all of it came from two rust surveillance corpora; the rest of
the unresolved cohort is controlled experiments where a single location is correct rather than
missing. There is little point chasing further supplements for geography.

## Results

The metadata module enriches all 1,285 BioProjects and 9,002 BioSamples from `runs.tsv`. Of those, 732 BioProjects (6,467 BioSamples) pass the full-text gate and are LLM-classified.

**Literature resolution coverage (1,286 BioProjects):**

| Coverage tier | BioProjects | % |
|--------------|-----------|---|
| Has PMID | 707 | 55.0% |
| DOI-only (no PubMed record) | 39 | 3.0% |
| **Any publication identifier** | **746** | **58.0%** |
| Has abstract | 721 | 56.1% |
| Has full manuscript text | 733 | 57.0% |

**BioSample XML field coverage (9,002 BioSamples, excluding NCBI placeholder values like `missing`/`not applicable`):**

| Field | BioSamples | % |
|-------|-----------|---|
| `geo_loc_name` | 5,035 | 55.9% |
| `collection_date` | 4,172 | 46.3% |
| `tissue` | 4,443 | 49.4% |

**LLM study design classification (732 full-text BioProjects):**

| Stress | BioProjects | | Coinfection intent | BioProjects |
|--------|-------------|-|---------------------|-------------|
| Biotic | 601 | | Single-pathogen focus | 530 |
| Abiotic | 94 | | Not disease-focused | 133 |
| None | 37 | | Intentional multi-pathogen | 69 |

The dominance of single-pathogen-focus studies (530/732, 72%) reflects the sampling design: both MAL and HAL query by known PHI-base pathogen species, selecting for experiments with a defined pathogen target. The 69 intentional-multi-pathogen BioProjects are flagged (`llm_coinfection_intent == "intentional_multi_pathogen"`) for exclusion from co-infection rate calculations, as their secondary detections are experimental rather than incidental.

**Setting effect on co-infection rate** (see `02_literature/03_classify/figures/sample_funnel_v3.py`, 6,467 classified BioSamples): field-collected BioSamples show an 11.1% biotic-only cryptic co-infection rate versus 4.0% in greenhouse and 7.8% in other controlled settings (growth chamber, detached-leaf assay, in vitro) — the field rate is roughly 2.8x the greenhouse rate, consistent with the ecological hypothesis that field samples encounter ambient pathogen pressure absent from controlled environments.

## Limitations

**Per-sample provenance is not verified, and this threatens the chapter's headline.** The
setting contrast above rests on `llm_study_setting == "field"` meaning field-collected, and
hand-checking shows it does not reliably mean that. Three fields each conflate two different
things:

- *Setting conflates study design with stated location.* PRJNA1217477 (211 samples, 8% of the
  analysed cohort) is a *Botrytis cinerea* inoculation atlas at UC Davis; PRJNA526829 (60) is
  *in vitro* on artificial surfaces; PRJNA328045 (39) is controlled compatible/incompatible
  inoculation. All are classified `field`. A title-keyword screen flags roughly 353 samples
  (13.4%) as plausibly non-field. The LLM judged these per-BioProject from full text and still
  got them wrong, so re-reading titles is not the fix; this needs a sample-level check.
- *Location conflates the submitting institution with the collection site.* "USA: California,
  Davis" (223 samples) is UC Davis, not a paddock. Geocoding cannot tell the two apart, because
  only the study design distinguishes them, and plotting institutions on a map as if they were
  field surveys is exactly the artefact that would manufacture the result.
- *Date conflates collection with isolation and deposit.* 83 samples are dated before 2010, back
  to 1925, with `isolation_source` values like "Federation" (a wheat cultivar released in 1901)
  and "spores". These are culture-collection and archival isolates, not field observations.

Where a supplement supplies an author-declared Field/Lab column, 10 samples the authors
themselves call Lab are classified `field` by the pipeline. Supplements are currently the only
source that catches this.

The designed fix is a per-sample verification assigning each BioSample three independent
classes rather than trusting one LLM field: `design` (field-collected | inoculated | in-vitro |
archival-isolate | mapping/other), `location_kind` (collection-site | institution |
country-only | none), and `date_kind` (field-collection | isolate-origin | submission-proxy).
The field co-infection rate would then be computed on `design == field-collected` only, and the
maps drawn on `location_kind == collection-site` only. The counts will move; that is the point.
It is designed, not built. Until it is, **do not quote the field co-infection rate, the
country/locality map, or the 2,643-sample denominator as settled.** The supplement mining and
geocoding above are the evidence that verification will draw on, not a substitute for it.

**Where the two location sources disagree.** For the 515 BioSamples with both a BioSample
`geo_loc_name` and an LLM-extracted location, 202 disagree (`geo_agreement`). BioSample is
currently trusted by default, but the direction of the disagreement, study site versus
submitting institution, has not been characterised.

**Geocoding misses.** 199 of 467 locality strings did not geocode. Most are plot codes and road
names that the settlement-class filter correctly rejected, but some are real regions it
over-rejected (for example "Southern New South Wales"). The nulls in `geocode_cache.json` are
worth a manual pass.

**`map_localities.R` subtitle is wrong.** It recomputes the country-only/no-location split
differently from the data and miscounts it. The 388 sites and 1,450 samples in the title are
correct. Fix or drop the subtitle before the figure is used.

**LLM classification errors.** GPT-4o-mini classification is not verified against ground truth for the full corpus. Each of the five judgment dimensions carries its own `llm_*_confidence`/`llm_*_rationale` pair (in `samples.tsv`) and should be consulted when individual BioProject/BioSample assignments are used in analysis, rather than trusting the label alone.

**Publication coverage.** 42.0% of BioProjects (540/1,286) have no publication identifier at all, and classification is further gated on full-text retrieval succeeding (57.0% of BioProjects) — so the analysable population (732 BPs) is a real subset of the full screened corpus, not all of it. Unresolved/no-full-text submissions skew toward data-only repositories, unpublished surveillance datasets, and multi-omics portals that do not link to a primary publication or whose publisher blocks automated + manual PDF retrieval.

**BioSample XML field coverage.** Geographic (56%), temporal (46%), and tissue (49%) metadata are available for only about half of samples (after excluding placeholder values), limiting spatial and temporal analyses of co-infection distribution.

**Host attribution granularity.** `named_hosts`/`named_pathogens` are extracted once per BioProject, not per BioSample — correct for genuinely BioProject-constant judgments (stress, setting, tissue) but only an approximation for the ~147 multi-host BioProjects, where the per-BioSample disambiguation pass (above) is needed for a confident per-sample host assignment; even then, only BioSamples with clear supporting metadata (`bs_host` etc.) get a confident resolution — the rest are correctly left `unresolved` rather than guessed.

## Key output files

| File | Contents |
|------|----------|
| `02_literature/01_search/data/bioprojects.json` | Title, description, submission/pub date, pmid/doi/pmcid, abstract, full_text — 1,286 BioProjects |
| `02_literature/01_search/data/biosamples.json` | BioSample XML attributes — 9,002 samples (incl. ENA/DDBJ via EBI API) |
| `02_literature/02_text/data/failed_dois.tsv` | BioProjects with a DOI but no full text retrieved by any automated strategy |
| `02_literature/03_classify/data/samples.tsv` | **Primary analysis input.** One row per biosample_representative BioSample, full-text BioProjects only — 6,467 rows. See CLAUDE.md's Output schemas section for the full column list. |
| `02_literature/03_classify/data/classify_cache.jsonl` | Per-BioProject LLM classification cache (resumable) |
| `02_literature/03_classify/data/host_disambig_cache.jsonl` | Per-BioSample host disambiguation cache |
| `02_literature/03_classify/figures/sample_funnel_v3.html` | Interactive Sankey: BioSample flow from full-text retrieval through tissue/setting/stress to co-infection outcome |
| `02_literature/02_text/figures/lit_resolution_alluvial.png` | Literature resolution flow through each strategy |
| `02_literature/02_text/data/supp_provenance.tsv` | Per-sample provenance mined from supplementary tables — 1,272 BioSamples |
| `02_literature/03_classify/data/cohort_provenance.tsv` | **Merged per-sample provenance.** Location, country, source, collection year, author-declared sample type — 2,643 rows |
| `02_literature/03_classify/data/provenance_tracker.tsv` | Worklist of BioProjects ranked by missing provenance, classed DONE / ok / TARGET / SETTING? / NO DOI |
| `02_literature/03_classify/data/localities.tsv` | 388 geocoded collection sites, plus the rejected strings |
| `02_literature/03_classify/figures/map_sample_origins.png` | Country-level map: proportional symbols, ranked bars, year histogram |
| `02_literature/03_classify/figures/map_localities.png` | Locality-level point map with a Europe inset (subtitle miscounts, see Limitations) |

All data outputs in this table are gitignored, per the repository's code-not-data policy. The
configured supplements themselves live in `02_text/data/supp_data/` and are likewise excluded:
they are publishers' copyrighted supplementary files.
