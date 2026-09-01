# metadata — Study Design Enrichment and LLM Classification

## Rationale

A co-infection detection in isolation carries limited biological meaning without study design context. The central confound is that not all positive detections represent ecologically authentic co-infections: a run from a co-inoculation experiment, an abiotic stress trial with a pathogen treatment, or a controlled greenhouse inoculation may show secondary signal by design, not by incidental co-infection. Conversely, a field-collected sample showing secondary signal is far more likely to represent a genuine unreported co-infection.

This module enriches the STAT detections from `runs.tsv` with three layers of context: (1) BioProject and BioSample metadata (including ENA/DDBJ BioSamples, which return empty from NCBI efetch), (2) linked literature and full manuscript text, and (3) LLM-based study design classification. Together these allow stratification of detections by study intent, setting, tissue, and stated pathogen/host — the necessary prerequisite for any epidemiological interpretation.

## Methods

### BioProject and BioSample metadata + literature linkage (`meta_search.py`)

`metadata/meta_search.py` resolves BioProject/BioSample metadata and literature identifiers in one atomic pass per BioProject, matching NCBI's actual resolution cascade: BioProject XML `<Publication>` field → PMC full-text search → PubMed title search — all NCBI-XML derived, including bare-DOI-no-PMID cases. Only if all three find nothing does it fall through to a Serper web search, then finally DOI extraction/page scrape/CrossRef/PMC-by-DOI. This replaces the old `ncbi_metadata.py` + `web_metadata.py` split.

BioSample XML is parsed for a harmonised field set (`geo_loc_name`, `tissue`, `collection_date`, `isolation_source`, `dev_stage`, `lat_lon`, `host`). Coverage is uneven — the SRA submission process does not mandate these fields, and a nontrivial fraction of populated fields are NCBI placeholder values (`missing`, `not applicable`, `not collected`) rather than real data — see Results below for filtered coverage numbers. `meta_search.py` also fetches ENA/DDBJ BioSamples (`SAME*`/`SAMD*` accessions, which return empty from NCBI efetch) via the EBI BioSamples API; its `bs_description` field (e.g. "RNA-Seq of a field sample of Wheat Yellow Rust") is a particularly strong signal fed directly into the LLM classification prompt.

### Full-text retrieval (`meta_text.py`)

`metadata/meta_text.py` retrieves full manuscript text, in cascade order: PMC OA full-text XML (free, complete, no PDF-parsing artifacts, when `meta_search.py` already resolved a PMCID) → Unpaywall → PDF download → `pdfminer` extraction → manual PDF fallback (`--ingest-manual`, for hand-downloaded PDFs matched against the failed-DOI list). `--apply` writes `full_text` back into `bioprojects.json` in place, which gates entry into classification below — **only BioProjects with retrievable full text are classified at all**, avoiding the old pipeline's failure mode of guessing study design from a title alone.

### LLM study design classification (`meta_classify.py`)

`metadata/meta_classify.py` replaces the old keyword-based classifier (`filter_kw.py`) and single-pass LLM classifier (`llm_classify.py`) entirely — keyword classification is not used at all in the current pipeline; the old approach's ~51% misclassification rate for `host_study` (pathogen/disease language co-occurring with host biology language) made it unreliable as anything but a discarded baseline.

Six LLM calls per BioProject (`gpt-4o-mini`): five independent judgment dimensions — `stress` (biotic/abiotic/none), `study_setting` (field/greenhouse/growth_chamber/detached_leaf_assay/in_vitro/unclear), `tissue` (aerial/non-aerial/unclear), `coinfection_intent` (single_pathogen_focus/not_disease_focused/intentional_multi_pathogen), and `hostpath` (named pathogens + named hosts, each as a list with its own confidence) — each with its own confidence + rationale, plus one fact-extraction call for symptom status, exposure type, geographic location, library prep, host cultivar, and host resistance. `--focus {stress|setting|tissue|coinfection|hostpath|extract}` reruns a single dimension.

A separate per-BioSample host disambiguation pass (`--disambiguate-hosts`) resolves which specific host a BioSample belongs to when a BioProject's `named_hosts` has more than one value (e.g. a rust fungus's alternate-host life cycle spanning two plant genera) — `named_hosts` itself is extracted once per BioProject and would otherwise be duplicated identically across every BioSample in a multi-host study, which is structurally wrong for per-sample analysis. The disambiguation prompt uses per-BioSample metadata (`bs_host`, `bs_description`, `tissue`, `isolation_source`) to pick the best-matching candidate, with its own confidence + rationale, and its own cache keyed by BioSample (not BioProject) — including the exact candidate list each resolution was made against, so a later prompt improvement that reshapes the candidates correctly invalidates stale resolutions rather than silently keeping them.

Author-named pathogens/hosts are resolved to NCBI taxids deterministically in plain Python afterward (`_util.resolve_taxon_name()`), never asked of the LLM — an LLM recalling taxids from memory produces confident-looking wrong numbers.

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

**Setting effect on co-infection rate** (see `metadata/figures/sample_funnel_v3.py`, 6,467 classified BioSamples): field-collected BioSamples show an 11.1% biotic-only cryptic co-infection rate versus 4.0% in greenhouse and 7.8% in other controlled settings (growth chamber, detached-leaf assay, in vitro) — the field rate is roughly 2.8x the greenhouse rate, consistent with the ecological hypothesis that field samples encounter ambient pathogen pressure absent from controlled environments.

## Limitations

**LLM classification errors.** GPT-4o-mini classification is not verified against ground truth for the full corpus. Each of the five judgment dimensions carries its own `llm_*_confidence`/`llm_*_rationale` pair (in `samples.tsv`) and should be consulted when individual BioProject/BioSample assignments are used in analysis, rather than trusting the label alone.

**Publication coverage.** 42.0% of BioProjects (540/1,286) have no publication identifier at all, and classification is further gated on full-text retrieval succeeding (57.0% of BioProjects) — so the analysable population (732 BPs) is a real subset of the full screened corpus, not all of it. Unresolved/no-full-text submissions skew toward data-only repositories, unpublished surveillance datasets, and multi-omics portals that do not link to a primary publication or whose publisher blocks automated + manual PDF retrieval.

**BioSample XML field coverage.** Geographic (56%), temporal (46%), and tissue (49%) metadata are available for only about half of samples (after excluding placeholder values), limiting spatial and temporal analyses of co-infection distribution.

**Host attribution granularity.** `named_hosts`/`named_pathogens` are extracted once per BioProject, not per BioSample — correct for genuinely BioProject-constant judgments (stress, setting, tissue) but only an approximation for the ~147 multi-host BioProjects, where the per-BioSample disambiguation pass (above) is needed for a confident per-sample host assignment; even then, only BioSamples with clear supporting metadata (`bs_host` etc.) get a confident resolution — the rest are correctly left `unresolved` rather than guessed.

## Key output files

| File | Contents |
|------|----------|
| `metadata/output/meta_search/data/bioprojects.json` | Title, description, submission/pub date, pmid/doi/pmcid, abstract, full_text — 1,286 BioProjects |
| `metadata/output/meta_search/data/biosamples.json` | BioSample XML attributes — 9,002 samples (incl. ENA/DDBJ via EBI API) |
| `metadata/output/meta_text/data/failed_dois.tsv` | BioProjects with a DOI but no full text retrieved by any automated strategy |
| `metadata/output/meta_classify/data/samples.tsv` | **Primary analysis input.** One row per biosample_representative BioSample, full-text BioProjects only — 6,467 rows. See CLAUDE.md's Output schemas section for the full column list. |
| `metadata/output/meta_classify/data/classify_cache.jsonl` | Per-BioProject LLM classification cache (resumable) |
| `metadata/output/meta_classify/data/host_disambig_cache.jsonl` | Per-BioSample host disambiguation cache |
| `metadata/output/figures/sankey/sample_funnel_v3.html` | Interactive Sankey: BioSample flow from full-text retrieval through tissue/setting/stress to co-infection outcome |
| `metadata/output/figures/sankey/lit_resolution_alluvial.png` | Literature resolution flow through each strategy |
