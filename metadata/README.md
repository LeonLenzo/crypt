# metadata — Study Design Enrichment and LLM Classification

## Rationale

A co-infection detection in isolation carries limited biological meaning without study design context. The central confound is that not all positive detections represent ecologically authentic co-infections: a run from a co-inoculation experiment, an abiotic stress trial with a pathogen treatment, or a controlled greenhouse inoculation may show secondary signal by design, not by incidental co-infection. Conversely, a field-collected sample showing secondary signal is far more likely to represent a genuine unreported co-infection.

This module enriches the STAT detections from `runs.tsv` with three layers of context: (1) BioProject and BioSample XML metadata, (2) linked literature and methods sections, and (3) LLM-based study design classification. Together these allow stratification of detections by study intent, setting, and host tissue — the necessary prerequisite for any epidemiological interpretation.

## Methods

### BioProject and BioSample metadata

`metadata/fetch_xml.py` retrieves BioProject XML and BioSample XML records for all accessions in `runs.tsv` via NCBI Entrez POST efetch (batched to avoid URL-length limits). BioSample attributes are parsed for a harmonised field set: `geo_loc_name`, `tissue`, `collection_date`, `age`, `sex`, `disease`, `treatment`, and `description`. Coverage is uneven — the SRA submission process does not mandate these fields — with `geo_loc_name` present in 61% of samples, `tissue` in 56%, and `collection_date` in 42%. The resulting `bioprojects.json` and `biosamples.json` feed all downstream classification steps.

### Literature linkage

`metadata/fetch_lit.py` attempts to resolve a primary PMID for each BioProject via a six-strategy short-circuit search, stopping at the first successful hit:

1. BioProject XML `<Publication>` field
2. PMC full-text search by BioProject accession
3. Europe PMC search by BioProject accession
4. ENA XML project record
5. Semantic Scholar title search (requires `S2_API_KEY`)
6. PubMed title search

Where a PMID is resolved, the PMC full-text methods section is retrieved and stored. This methods text is later used in LLM classification to distinguish study designs that are ambiguous from the BioSample metadata alone (e.g., distinguishing a field transcriptomics study from a controlled inoculation experiment when both use the same pathogen and host).

### Keyword study design classification

`metadata/filter_kw.py` assigns each BioSample to a treatment category (`single`, `host_study`, `abiotic_stress`, `coinf_experiment`, `surveillance`, `combined_stress`, `unclear`) and setting (`field`, `lab`, `unclear`) using a curated keyword vocabulary applied to BioSample attributes and BioProject title/description. Keyword classification has known limitations: approximately 51% of BioProjects assigned `host_study` by keyword actually belong to other categories, because pathogen and disease language co-occurs with host biology language in study descriptions. Keyword results are retained in `biosample_kw.tsv` but the LLM classification is preferred for downstream analysis.

### LLM study design classification

`metadata/llm_classify.py` submits each BioProject's title, description, and (where available) PMC methods text to GPT-4o-mini for classification into the same treatment and setting vocabulary. The LLM prompt explicitly distinguishes between studies that inoculate plants with a known pathogen (controlled) and studies that collected plants from field environments (field). Classification is cached per BioProject in `llm_classify/data/classify_cache.jsonl` to enable resumable runs and avoid re-billing already-classified entries.

## Results

The metadata module enriches all 1,285 BioProjects and 9,002 BioSamples from `runs.tsv`.

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

**PMID resolution coverage.** Despite six search strategies, a substantial proportion of BioProjects cannot be linked to a publication — either because the study is unpublished, the submission did not include a publication link, or the title-based searches returned false positives and were rejected. Studies without PMIDs skew toward surveillance datasets and unpublished data repositories, where field conditions are actually common. Missing literature therefore likely biases the `field` count downward.

**BioSample XML field coverage.** Geographic (`geo_loc_name`: 61%), tissue (56%), and temporal (`collection_date`: 42%) metadata are available for only a subset of samples, limiting spatial and temporal analyses of co-infection distribution. The `tissue` field in particular would be valuable for distinguishing whether co-infections are tissue-specific (e.g., root vs. leaf pathogens) but cannot be used systematically given incomplete coverage.

## Key output files

| File | Contents |
|------|----------|
| `metadata/output/fetch_xml/data/bioprojects.json` | BioProject XML records for 1,285 BioProjects |
| `metadata/output/fetch_xml/data/biosamples.json` | BioSample XML attributes for 9,002 samples |
| `metadata/output/fetch_lit/data/literature.json` | Resolved PMIDs and PMC methods text |
| `metadata/output/filter_kw/data/biosample_kw.tsv` | 9,002-row biosample table; keyword treatment/setting + metadata |
| `metadata/output/llm_classify/data/bioproject_llm.tsv` | 1,285-row BioProject table; `llm_treatment`, `llm_study_setting`, `llm_rationale` |
