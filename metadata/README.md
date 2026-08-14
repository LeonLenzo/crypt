# metadata — Study Design Enrichment and LLM Classification

## Rationale

A co-infection detection in isolation carries limited biological meaning without study design context. The central confound is that not all positive detections represent ecologically authentic co-infections: a run from a co-inoculation experiment, an abiotic stress trial with a pathogen treatment, or a controlled greenhouse inoculation may show secondary signal by design, not by incidental co-infection. Conversely, a field-collected sample showing secondary signal is far more likely to represent a genuine unreported co-infection.

This module enriches the STAT detections from `runs.tsv` with three layers of context: (1) BioProject and BioSample XML metadata, (2) linked literature and methods sections, and (3) LLM-based study design classification. Together these allow stratification of detections by study intent, setting, and host tissue — the necessary prerequisite for any epidemiological interpretation.

## Methods

### BioProject and BioSample metadata

`metadata/fetch_xml.py` retrieves BioProject XML and BioSample XML records for all accessions in `runs.tsv` via NCBI Entrez POST efetch (batched to avoid URL-length limits). BioSample attributes are parsed for a harmonised field set: `geo_loc_name`, `tissue`, `collection_date`, `age`, `sex`, `disease`, `treatment`, and `description`. Coverage is uneven — the SRA submission process does not mandate these fields — with `geo_loc_name` present in 61% of samples, `tissue` in 56%, and `collection_date` in 42%. The resulting `bioprojects.json` and `biosamples.json` feed all downstream classification steps.

### Literature linkage

Literature resolution uses a DOI-first schema: `primary_doi` is the canonical publication identifier (universal across journals, preprints, and data notes), with `primary_pmid` as a secondary identifier when the DOI maps to a PubMed record. This schema extension from PMID-only was necessary to capture the ~2.5% of BioProjects linked only to preprints or data papers not indexed by PubMed.

`metadata/fetch_lit.py` resolves a primary DOI and/or PMID for each BioProject via a six-strategy short-circuit search, stopping at the first successful hit:

1. BioProject XML `<Publication>` field
2. PMC full-text search by BioProject accession
3. Europe PMC search by BioProject accession
4. ENA XML project record
5. Semantic Scholar title search (requires `S2_API_KEY`)
6. PubMed title search

Where a PMID is resolved, the PMC full-text methods section is retrieved and stored. Three supplementary resolution passes handle BioProjects not resolved by `fetch_lit.py`:

- **`backfill_doi.py`**: batch-fills `primary_doi` for the ~930 BioProjects that had a PMID but no DOI, via NCBI efetch ArticleIdList. Ensures DOI is populated for all PMID-linked records.
- **`serper_resolve.py`**: reads saved Google Serper search results (`metadata/output/serper/serper_results.tsv`), extracts DOIs from academic domain links and snippets, resolves to PMIDs via PubMed → EuropePMC, and enriches DOI-only results via CrossRef. Resolves ~64 additional BioProjects.
- **`serper_scrape.py`**: for BioProjects with Serper hits that did not yield a DOI via link extraction alone, attempts page scraping and title search. Resolves ~39 additional BioProjects.

Where no PMID is available, `primary_doi` alone is sufficient: CrossRef provides title, abstract, and publication date for any DOI, and the abstract is stored for downstream LLM classification. A final manual DOI entry step resolved 31 bot-blocked BioProjects (ScienceDirect and ResearchGate pages) from the Serper unresolved set.

**A critical validation note**: PubMed `esearch` with a `[doi]` term performs fuzzy matching. For ResearchSquare preprint DOIs (pattern `10.21203/rs.3.*`), it returns `count > 1` with unrelated PMIDs. All DOI-to-PMID lookups therefore require `count == 1` before accepting the returned PMID.

### Keyword study design classification

`metadata/filter_kw.py` assigns each BioSample to a treatment category (`single`, `host_study`, `abiotic_stress`, `coinf_experiment`, `surveillance`, `combined_stress`, `unclear`) and setting (`field`, `lab`, `unclear`) using a curated keyword vocabulary applied to BioSample attributes and BioProject title/description. Keyword classification has known limitations: approximately 51% of BioProjects assigned `host_study` by keyword actually belong to other categories, because pathogen and disease language co-occurs with host biology language in study descriptions. Keyword results are retained in `biosample_kw.tsv` but the LLM classification is preferred for downstream analysis.

### LLM study design classification

`metadata/llm_classify.py` submits each BioProject's title, description, and (where available) PMC methods text to GPT-4o-mini for classification into the same treatment and setting vocabulary. The LLM prompt explicitly distinguishes between studies that inoculate plants with a known pathogen (controlled) and studies that collected plants from field environments (field). Classification is cached per BioProject in `llm_classify/data/classify_cache.jsonl` to enable resumable runs and avoid re-billing already-classified entries.

## Results

The metadata module enriches all 1,285 BioProjects and 9,002 BioSamples from `runs.tsv`.

**Literature resolution coverage (9,005 BioSamples; 1,287 BioProjects in runs.tsv):**

| Coverage tier | BioSamples | % |
|--------------|-----------|---|
| Has PMID | 6,279 | 69.8% |
| DOI-only (no PubMed record) | 229 | 2.5% |
| **Any publication identifier** | **6,508** | **72.3%** |
| Has abstract | 6,421 | 71.3% |
| Has methods text (PMC full-text) | 5,393 | 59.9% |
| Unresolved | 2,496 | 27.7% |

The DOI-first refactor expanded total publication coverage from ~65% (PMID-only strategy) to 72.3% of BioSamples. DOI-only BioProjects now have CrossRef abstracts, enabling LLM classification of ~229 additional BioSamples previously inaccessible to the PMC-only methods extraction.

**LLM treatment classification (1,285 BioProjects):**

| Treatment | BioProjects |
|-----------|-------------|
| Single pathogen | 895 |
| Host biology study | 219 |
| Abiotic stress | 94 |
| Co-infection experiment | 21 |
| Surveillance | 6 |
| Unclear | 50 |

The dominance of single-pathogen studies (70%) reflects the sampling design: both MAL and HAL query by known PHI-base pathogen species, selecting for experiments with a defined pathogen target. The 21 co-infection experiments are flagged for manual exclusion from co-infection rate calculations, as their secondary detections are experimental rather than incidental.

**Setting effect on co-infection rate:**

LLM study setting classification reveals a clear ecological signal. Field-classified BioProjects show a co-infection rate of 22.7% (277/1,219 biosample-representative runs) versus 11.4% (531/4,650) in controlled laboratory settings. This approximately twofold difference is consistent with the ecological hypothesis: field samples encounter ambient pathogen pressure absent from controlled environments, and co-infecting organisms present at the time of sampling are captured in the sequencing library alongside the target pathogen.

## Limitations

**LLM classification errors.** GPT-4o-mini classification is not verified against ground truth for the full corpus. Manual review of a subset of BioProjects showed broadly accurate treatment assignments, but edge cases exist — particularly where the BioProject description is sparse and no PMC methods text was available. The `llm_rationale` column in `bioproject_llm.tsv` retains the model's stated reasoning for each classification and should be consulted when individual BioProject assignments are used in analysis.

**Publication coverage.** Despite six primary strategies plus three supplementary passes, 27.7% of BioSamples (2,496 BioSamples; 517 BioProjects) remain unresolved. Unresolved submissions skew toward data-only repositories, unpublished surveillance datasets, and multi-omics portals (OmicsDI, GOLD, SEQout) that do not link to a primary publication. Studies without any publication identifier likely bias the `field` count downward, as field monitoring datasets are disproportionately represented in the unpublished tier.

**BioSample XML field coverage.** Geographic (`geo_loc_name`: 61%), tissue (56%), and temporal (`collection_date`: 42%) metadata are available for only a subset of samples, limiting spatial and temporal analyses of co-infection distribution. The `tissue` field in particular would be valuable for distinguishing whether co-infections are tissue-specific (e.g., root vs. leaf pathogens) but cannot be used systematically given incomplete coverage.

## Key output files

| File | Contents |
|------|----------|
| `metadata/output/fetch_xml/data/bioprojects.json` | BioProject XML records for 1,285 BioProjects |
| `metadata/output/fetch_xml/data/biosamples.json` | BioSample XML attributes for 9,002 samples |
| `metadata/output/fetch_lit/data/literature.json` | Resolved DOIs, PMIDs, abstracts, and PMC methods text (tracked) |
| `metadata/output/fetch_lit/data/lit_cache.json` | Full resolution cache with provenance fields (gitignored; 1.9 GB) |
| `metadata/output/serper/serper_results.tsv` | Raw Google Serper results for 889 unresolved BioProjects |
| `metadata/output/serper/unresolved_hits.tsv` | 198 BPs still unresolved after all strategies; Serper top-3 links per BP |
| `metadata/output/filter_kw/data/biosample_kw.tsv` | 9,002-row biosample table; keyword treatment/setting + metadata |
| `metadata/output/llm_classify/data/bioproject_llm.tsv` | 1,285-row BioProject table; `llm_treatment`, `llm_study_setting`, `llm_rationale` |
| `metadata/output/figures/sankey/lit_resolution_sankey.html` | Interactive Sankey: BioSample flow through each resolution strategy |
