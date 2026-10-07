# Curating one BioProject

Written 2026-10-06 because curating a single project was costing ~25 tool calls, most of
them rediscovering the plumbing below. Target is 8-12. Read this file instead of re-reading
`apply_curation.py`.

## Order of operations

1. `python 01a_Literature/project_brief.py PRJNAxxxxxx` - everything known locally, one call.
   It already lists candidate papers from `accession_papers.tsv`. **An empty `primary_paper`
   does NOT mean no paper is known.** Check the brief before planning any search.
2. If no candidate: Europe PMC accession search
   `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%22PRJNAxxxxxx%22&format=json&pageSize=25`
3. Fetch the ONE paper that carries the field methods, as `.../{PMCID}/fullTextXML`, and file
   it at `studies/doi_<doi with / -> _>/`. A second paper on the same deposit needs only its
   data-availability line quoted, not a full copy.
4. Read only: the growth/field setup section, the library construction section, and a
   keyword COUNT scan (not a context dump) for `fungicide|pesticide|disease|pathogen|infect|
   symptom|inocul`. Open context only where a count is non-zero.
5. Write the rows (below). Per-sample values need a resolver; constant ones do not.
6. `apply_curation.py --dry-run`, grep for your own rule ids plus `warn|refus`, then apply.

## The five tables

    provenance.tsv    evidence_id  paper_key  doi  locus  quote  entered_by  ts
    curation.tsv      rule_id  BioProject  predicate  field  value  evidence_id  override  note
    joins.tsv         rule_id  BioProject  dir  file  sheet  header_row  key_col  key_regex
                      sra_key  value_col  value_prefix  field  evidence_id  override  note
    decisions.tsv     human verdicts only; review.record() raises PermissionError if `who`
                      starts with "claude" and the decision is not "retracted"
    found.tsv         ts  paper_key  kind  value  note  found_by
                      kinds in use: note lead negative blocked trap prior-art paper reject
                      accession bioproject supplement species-code version

Settable fields, each with a `_source` twin written automatically:
`tissue`, `location`, `collection_date`, `setting`, `sampling_selection`,
`library_representative`.

## Precedence, so you never need `override`

    0  ""/absent/not_fetched    1  biosample    2  _import    3  _join    4  _rule

A higher tier replaces a lower one freely. A join or a hand-written rule therefore beats the
archive without any flag; `override: yes` is only for replacing an EQUAL or higher tier.
Existing `override: yes` rows mostly predate the tier system.

## Predicates and join keys

`curation.tsv predicate`: `field=value` or `field~regex`, joined by `;`, evaluated over the
BioSample attributes first and then the runinfo columns. **Empty predicate = the whole
BioProject**, which is the common case for a single-site single-season trial.

`joins.tsv sra_key`: a bare column name, or `column~regex` with the key in group 1. Resolution
is BioSample attributes first, then runinfo. So a design code living in a BioSample attribute
(`source_material_id`) is a zero-regex join key.

Joins below `MIN_JOIN_RATE = 0.60` refuse rather than merge nothing and report success.

## Resolvers

One per study, at `01a_Literature/<name>_runs.py`, writing a CSV into that study's
`studies/doi_*/` directory. Only write one when a value varies PER SAMPLE. Every resolver
carries refusal guards that `sys.exit("REFUSED: ...")` on the assumption that would silently
corrupt the result (a flattened time series, a non-unique join key, an unmapped tissue code).

`studies/**` and `data/*` are gitignored: the repo holds code and the five tables, never
papers or run data. Tracked tables are curation, joins, provenance, found, decisions,
host_summary, run_species. Verify with `git check-ignore` before staging.

## Output discipline

The expensive mistakes are printing, not thinking.

- Value-counts before raw records. Never `grep` a JSON; load it and count.
- Keyword scans print counts first; pull context only for the hits that matter.
- Cap context windows at ~200 chars and 2 hits per term.
- Do not verify a number that cannot change a decision. An off-by-one in a locality count is
  not worth a tool call.
- One worked example row is enough; never dump a column's whole vocabulary.
- Commit messages: what changed and the one surprising thing, not the full evidence chain.
  The chain is already in `provenance.tsv`.

## Where each fact comes from, in order

Try these in order and stop at the first that answers. Each line's failure mode is why the
next one exists. `sources.py` wraps the endpoints; the rule is to never retype a URL.

    setting       the paper's growth/field-setup section. ONLY the paper.
                  Never an archive field, never an author affiliation, never an abstract:
                  "naturally occurring water-limited conditions" and "25% field capacity"
                  both describe greenhouses. ~22% of everything assessed was not field.
    location      the paper's site statement, with coordinates where given
                  -> GEO characteristics
                  -> `geo_loc_name` ONLY once a paper confirms it is a site. It has meant
                     submitter country, germplasm origin, and a per-sample counter.
    collection_date  a per-sample design code in the BioSample attributes, joined to dates
                  the paper states (`project_brief.py` flags a unique code automatically)
                  -> GEO `description` codes (EPICON's are MMDDYY)
                  -> the paper's stated range, if one date per sample is unavailable
                  -> `collection_date`. Last, not first: it has been malformed (`2014-24-06`),
                     three years out, and identical across a weekly time series.
    tissue        BioSample `tissue` is usually right but has been confidently wrong
                  (`leaf` on 57 hydroponic root libraries) -> GEO `source_name_ch1`
                  -> the paper. Decide root-vs-aerial from the paper when they disagree.
    sampling_selection  `sources.py scan` for the vocabulary, then read only the non-zero
                  hits. Absence across the whole list is what licenses `unselected`, so
                  state the scan in the note. Scan case-INSENSITIVELY: a case-sensitive
                  pass reported zero `disease` hits for a paper that had one.
    host species  ENA portal per project -> the study's own supplement. `scientific_name` is
                  the submitter's subject: 929 infected wheat leaves are filed under
                  *Puccinia striiformis*.
    library chemistry  the paper's library-construction section, or GEO `molecule`.
                  `library_selection` says `cDNA` for both polyA and ribodepleted.

## Before reading any paper

`project_brief.py` then `offtarget_screen.py`. A WGS, amplicon or Hi-C library cannot be
polyA mRNA whatever the organism, but the strategy field itself can lie: PRJNA1314945 is 306
SLAF-seq genotyping runs all registered as RNA-Seq/cDNA/TRANSCRIPTOMIC, caught only by a
DNA-assay word in the paper title. Organism is never decisive on its own.

## The layout: seeds, data, gold

Three incidents on 2026-10-06 all had one cause - nothing in the layout said which files may
be rebuilt. `papers.tsv` was hand-annotated and the next `triage.py` run discarded it; an
imported accession list survived only as a shell argument; and the file written to fix that
was gitignored, because tracking depended on `git add -f` habits rather than on file kind.

    seeds/   HAND-WRITTEN, irreplaceable, TRACKED by default.
             curation.tsv  joins.tsv  provenance.tsv  decisions.tsv  found.tsv
             run_scope.tsv  backfill_accessions.tsv      (+ _exclusions.py, repo root, code)
    data/    generated or cached, IGNORED. runs.tsv, bioprojects.tsv, papers.tsv,
             paper_bioproject.tsv, registry.tsv, worklist.tsv, and the bronze caches
             runinfo/, biosample_attrs/, geo/, ena/.
    gold/    generated but small and tracked, so summary numbers have a history.
             cohort.tsv  host_summary.tsv  offtarget.tsv  barley_gap.tsv  frame_gap.tsv

**The rule: if losing a file costs a refetch it is data/; if it costs re-reading papers it is
seeds/.** Paths come from `_layout.py` - never hardcode `HERE / "data" / ...`. The repo root
has its own `_paths.py` for cross-module paths; the two must not be confused, which is why
this one is `_layout`.

Anything worth keeping goes in a seed table. A fact about a paper goes in `provenance.tsv` or
`found.tsv`, never in `papers.tsv`.

## The cohort is defined once

`cohort.py` holds the predicate (field* setting, representative, aerial tissue) and nothing
else may restate it. `from cohort import cohort, all_runs`. It was being retyped into every
ad-hoc query - four times on 2026-10-06 - and one copy drifted, giving 425 localities in one
script and 424 in another. `python 01a_Literature/cohort.py` prints the counts;
`--write` materialises `gold/cohort.tsv`.

`sampling_selection` is deliberately NOT part of the predicate. Every stratum is in the
cohort; the flag decides what counts as a finding, not whether a sample counts.

`found.tsv` is also an INPUT: a row with `kind = bioproject` feeds `found_accessions()` and
becomes a `cited` claim on the next rebuild. That is the designed way to tell the pipeline
about an accession found while reading.

**`relation = generated` is a human verdict and not ours to write.** `found_accessions()`
deliberately emits `cited` only, because who produced a deposit is Leon's call; record the
data-availability quote as evidence and leave the attribution to him.

## Rebuilding after adding accessions

    python 01a_Literature/triage.py --from-undermind --accessions PRJ... # NEVER --accessions alone
    python 01a_Literature/import_provenance.py                          # restores import-* tier
    python 01a_Literature/apply_curation.py                             # restores rule/join tiers

`--accessions` without `--from-undermind` rebuilds the universe from that list plus the Kraken
cohort only, dropping every literature-first candidate. Back up `runs.tsv` and
`bioprojects.tsv` first: they are gitignored, so there is no git safety net. Then diff
per-run values against the backup and confirm ZERO regressions before moving on.
