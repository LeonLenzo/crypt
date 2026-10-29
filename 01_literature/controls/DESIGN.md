# Negative controls: which libraries, and what they are for

Sketched 2026-10-23. Not built. Supersedes nothing; it makes the plan in
`crypt-controls-design` concrete and narrows it to pathogens and hosts. Decisions taken while
sketching it are marked in place and dated.

## The one number the controls have to settle

The saturation floor is the last provisional parameter in the method, and it is the only
floor that still moves the headline. Measured on the 1,332 runs classified so far, holding
the coverage and distinct-minimizer floors fixed and varying saturation alone:

| saturation floor | pathogen calls | runs with >=1 | co-infection, species | co-infection, genus |
|---|---|---|---|---|
| 0.00 to 0.05 | 2,970 | 1,216 | 54.0% | **46.8%** |
| 0.075 | 2,743 | 1,211 | 51.9% | 45.2% |
| 0.10 | 2,442 | 1,206 | 47.3% | 41.7% |
| 0.15 | 1,993 | 1,173 | 35.1% | 30.9% |
| 0.20 | 1,655 | 1,146 | 23.6% | 20.8% |

Two things follow. The headline is **flat below 0.05**, which is the same insensitivity
already recorded for the absolute floor, so nothing is gained by arguing about small values.
Above 0.05 it falls steeply, and the difference between a floor at 0.05 and one at 0.15 is
the difference between reporting 47% and reporting 31%.

So the controls answer one question: **does cross-mapping from a single species produce
false-positive calls that stay below 0.05 saturation, or do they reach 0.10 to 0.15?**
Everything below is in service of measuring that, and nothing else needs to be measured
to the same precision.

## What a control measures

A library containing one known organism is classified against db_v4 unchanged. Every species
called that is not that organism is a false positive produced by cross-mapping, and it comes
with the same four quantities a cohort call does: reads, distinct minimizers, coverage and
saturation. The distribution of those false positives, specifically its upper tail in
saturation, is the floor.

This is the whole reason an axenic culture is not the control set dropped on 2026-08-28. That
set was arbitrary SRA runs with unknown ground truth. Here the truth is known by construction.

## THE PRIMARY CONTROL IS A CONTROLLED SINGLE-PATHOGEN INFECTION

Leon, 2026-10-26, and it supersedes axenic culture as the first choice rather than supplementing
it. A plant inoculated with one pathogen under controlled conditions, with no other pathogen
present, is a tighter control than a pure culture for three reasons.

**It matches the condition the claim is made in.** The cohort is in planta poly(A) RNA-seq.
Saturation depends on expression breadth, and a fungus growing in culture expresses a different
gene set from one colonising a leaf. A culture calibrates the floor under the wrong regime and
the size of that error is unknown.

**It measures host cross-mapping in the same library.** A pathogen-only and a host-only control,
measured apart, cannot show what the host's reads do to the pathogen's k-mer space when both are
present at realistic ratios. A single-pathogen infection does, because it is that situation with
the answer known.

**It dissolves the obligate biotroph problem completely.** *Puccinia* and *Blumeria* cannot be
grown axenically at any price, which is why the earlier draft fell back on host-free spore
libraries and labelled them a weaker tier. Inoculated rust and mildew time courses are ordinary
experiments. Measured against the cached ENA metadata, candidate runs exist for
*P. striiformis* (20), *P. graminis* (21), *P. triticina* (14) and *B. graminis* (43). **The
spore tier is therefore dropped**, not demoted.

### What it costs, and why that is the right direction

A greenhouse plant is not sterile. Endophytes and seed-borne fungi exist, so a call of species
Y in a library inoculated with X may be a real organism rather than cross-mapping from X. That
makes this an UPPER BOUND on the false-positive rate, which is the conservative direction: if
the measured tail still sits below 0.05, the floor is safe, and the bound being loose only makes
that conclusion stronger. Axenic culture is the one thing that separates cross-mapping from
genuine contamination, which is why it is retained below rather than discarded.

### The host control comes free, and is better than the old one

A controlled infection experiment almost always sequences mock-inoculated or uninoculated
plants alongside. That arm is the host-only control, from **the same plants, the same growth
conditions, the same library preparation and usually the same flow cell**. It is strictly better
than the separately sourced lab-grown host libraries the earlier draft specified, because every
difference except the pathogen is held constant. Prefer the paired mock arm; fall back on
independent lab-grown host libraries only where a study has none.

### The cohort's own inoculated runs do NOT qualify

215 cohort runs are curated `sampling_selection = inoculated` and the funnel already excludes
them from the denominator, so they look like a free control set. **They are not.** All but six
are FIELD inoculations (`field`, `field (inoculated in situ)`): the plant carries the inoculum
and everything else the field supplied, so there is no single-pathogen ground truth. The cohort
selected field studies by construction, so it contains no controlled laboratory infections at
all and the control set has to be built fresh from the archive.

### What this changes in the screen

The discovery query inverts. The earlier cascade looked for culture and treated in planta cues
as a veto; it now looks for a controlled single-pathogen infection and the vetoes become
field collection, natural infection, co-inoculation with a second organism, and any mixed or
consortium design. The evidence standard gains two requirements: the methods must state a
controlled environment (growth chamber, greenhouse or climate room) and a single-isolate
inoculum, and must not report a second pathogen in the same plants.

## What counts as a control library

The requirement is not axenic growth, it is **single-source RNA**: reads from one organism and
nothing else db_v4 contains. Axenic culture is the usual way to get that, but it is the means
rather than the requirement, and saying it this way makes the failure modes legible. A
strictly axenic flask whose library picked up index hopping from an in planta sample on the
same flow cell is useless here, and the purity of the flask was never the question.

A run has to pass both parts.

**Evidence standard, from the methods and recorded as a quoted sentence** the same way a
cohort value is: grown in pure culture from a single-spore or single-colony isolate, on a
medium containing no material from a species in db_v4, harvested with no host tissue present,
and sequenced in a library prepared from that culture alone. No statement means not admitted.
Silence is not evidence of purity.

**Empirical screen, on the data**: classified against db_v4, the host clades and the other
pathogen clades are absent rather than low. A library curated as pure that calls wheat is not
pure, whatever the methods say. Metadata is what has misled this project most often, so the
empirical screen runs before the curated label is believed, and a failed screen sends the run
back to curation with the evidence attached.

### Three tiers, by how tightly the truth is known

1. **Pure culture, defined medium.** Clean, no caveat.
2. **Pure culture, complex medium.** PDA, PDB, V8, oatmeal and malt extract all contain
   autoclaved plant or yeast material. Almost certainly harmless, since autoclaved extract
   holds little intact nucleic acid and poly(A) selection captures less of what remains, but
   it is an assumption and the empirical screen is what tests it. Media made with wheat straw
   or leaf extract are excluded outright, and those do appear in *P. nodorum* work.
3. **Controlled single-pathogen infection.** Not axenic, and now the PRIMARY class rather than
   a fallback. See the section above. Ground truth is "one pathogen was introduced and no other
   was reported", which is weaker than a pure culture but is measured in the right condition.

The host-free spore tier that stood here is **dropped**, 2026-10-26. It existed only because the
biotrophs had no axenic option, and controlled infection gives them a better one. Its weakness
was that it measured cross-mapping and host carry-over together with no way to separate them.

Tiers 1 and 2 are retained for a job tier 3 cannot do: an axenic culture is the only library in
which a call of any other species MUST be cross-mapping, because nothing else is present. Use
culture to measure the cross-mapping spectrum cleanly where it exists, and controlled infection
to bound the false-positive rate under realistic conditions everywhere.

### The circularity to avoid

The empirical screen uses db_v4 to validate the controls that calibrate db_v4's floor. That is
circular if the screen is written in terms of the saturation floor under calibration. It is not
circular if purity is judged on the **host clade alone**, which is a large k-mer space where
contamination is unmissable and whose threshold is not the parameter being tested. The screen
is to be written that way deliberately.

## Selection rule, and why it is not the list of suspect calls

A culture of species X reveals what **X's reads do to other species' k-mer spaces**. It says
nothing about whether a call of X elsewhere is real. The control set is therefore the species
believed to be genuinely present, not the species whose calls are in doubt.

The practical consequence is that several open questions are answered by controls chosen for
other reasons. *Stemphylium lycopersici* (35 runs, 32 of them wheat, 2.1% coverage, no sibling
of the genus called) needs no *Stemphylium* culture: if it is smear, it will appear in the
wheat pathogen controls or in the wheat host control. The same holds for the *Alternaria*
cluster and for *Epicoccum* and *Cladosporium* if they are ever questioned.

## Discovery and admission: archive first, paper standard

The cohort had to be built paper first because it is a sample of the literature, and entering
from the archive would have let annotation quality decide which studies existed. **The control
set has no denominator.** It is not a sample of anything and nothing is estimated from it, so
the only methodological commitment that has to hold is the evidence standard above. The
direction of entry does not touch that chain, it only changes the order in which it is walked.
Archive first also begins from a hard requirement, since deposited data is non-negotiable for a
control, whereas paper first means reading the in vitro transcriptomics literature and then
discovering whether anything was deposited.

So: **discovery from the archive, admission on the paper.** The cascade is the same
cheap-to-expensive shape used for the STAT tail.

1. ENA portal metadata screen, one request per taxon, to get the candidate runs.
2. Sample and study title screen, to shortlist. **This is a shortlist, not a decision.**
   `invitro_SN15_000` is the submitter's word for it and carries no audit weight.
3. Recover the paper, then read the methods and take the quote. Two routes, both needed:
   the ENA study XML's `XREF_LINK` to PUBMED where it exists, and a Europe PMC full-text
   search on the accession where it does not. Tested on the two *P. nodorum* studies and both
   resolved, one by each route, so neither route alone is sufficient.
4. Purity screen on the classified output, as above.

Only steps 3 and 4 admit anything. Unpublished submissions cannot be admitted, since there are
no methods to quote, and mutant RNA-seq deposited ahead of publication is common enough that
this will discard real candidates.

**Keep the requirement of at least two independent studies per species.** Representativeness
is not the point, but strain, medium and library prep are the obvious confounders, and taking
three or four runs from a single study spans none of them.

## Screen result, 2026-10-23

Stages 1 and 2 run, by `screen_candidates.py`, over 16 taxa and 6,117 RNA-seq runs. Raw ENA
pulls are cached under `raw/`, verdicts are in `candidates_runs.tsv` with the cue and the field
it fired on, and study-level rollups with resolved papers are in `candidates_studies.tsv`.
497 studies, of which 56 carry two or more candidate runs.

**Ready for the methods read.** Two or more studies with a resolved paper, six or more
candidate runs between them:

| species | candidate runs | studies | with paper |
|---|---|---|---|
| *Fusarium graminearum* | 281 | 27 | 16 |
| *Zymoseptoria tritici* | 77 | 6 | 6 |
| *Puccinia striiformis* (spore) | 48 | 7 | 5 |
| *Fusarium verticillioides* | 26 | 3 | 2 |
| *Bipolaris sorokiniana* | 12 | 2 | 2 |
| *Blumeria graminis* (spore) | 12 | 2 | 2 |
| *Exserohilum turcicum* | 6 | 2 | 2 |

That already covers wheat on four pathogens plus *F. graminearum*, and maize on two, so the
set is usable before the remainder is resolved.

**Single study, so the strain and medium confounders are unspanned.** *B. maydis* (20
candidates across 2 studies, neither with a recoverable paper), *A. alternata* (9),
*P. nodorum* (8), *P. teres* (5), *P. graminis* (3), *P. tritici-repentis* (2).

**No candidates at all.** *P. triticina*, which matters because it carries 298 cohort runs:
245 of its 418 "RNA-Seq" runs are `library_source=GENOMIC` and the remaining 138 have no
informative sample text. *Curvularia clavata*, 27 runs and not one with usable sample text.
*A. tenuissima*, 13 runs. For these three the metadata route is exhausted and only the papers
can decide, which is the cascade working as intended rather than a dead end.

### Two things the screen taught about itself

**Study-level text must not vote.** The first version let `study_title` into the verdict and
returned zero culture candidates for *P. nodorum*, vetoing all eight `invitro_SN15_*` runs on
the phrase "in planta", which appears in the title of the study that compared the two. The same
defect suppressed *Z. tritici*, *B. maydis* and *E. turcicum*. *P. nodorum* was the one species
whose answer was already known, so it functioned as a positive control for the screen itself
and is worth keeping in that role whenever the cues change.

**`host` and `isolation_source` describe the isolate, not the flask.** Vetoing on them would
discard every culture of a wheat-derived isolate. They are carried as `isolate_origin` for the
curator and never vote. Same trap as `crypt-isolate-vs-sample-trap`.

### Paper recovery needs all three routes

Over the 56 shortlisted studies: 8 by the ENA `XREF_LINK` to PUBMED, 28 by a Europe PMC
full-text search on the BioProject, 4 by the same search on the GEO series or SRP accession,
and 16 unresolved. The GEO route is not optional, since older *Fusarium* and *Zymoseptoria*
submissions are cited by their GEO series rather than their BioProject. The 16 unresolved are
presumed unpublished and are inadmissible under the evidence standard, which includes
PRJNA976811 and PRJNA436064, two of the richest candidate pools in the screen.

## Pathogen controls

RNA-seq run counts from the ENA portal, `tax_tree(taxid) AND library_strategy="RNA-Seq"`,
species node not strain node, surveyed 2026-10-23. Counts are an upper bound on what is
usable: they mix culture with in planta, and separating those is the curation work.

| species | cohort runs | median sat. | ENA RNA-seq | control is |
|---|---|---|---|---|
| *Puccinia striiformis* | 929 | 0.452 | 1,395 | host-free spores only |
| *Zymoseptoria tritici* | 310 | 0.211 | 549 | axenic culture |
| *Puccinia triticina* | 298 | 0.227 | 418 | host-free spores only |
| *Blumeria graminis* | 159 | 0.202 | 176 | host-free spores only |
| *Puccinia graminis* | 96 | 0.354 | 203 | host-free spores only |
| *Exserohilum turcicum* | 71 | 0.290 | 138 | axenic culture |
| *Parastagonospora nodorum* | 70 | 0.118 | 16 (14 real, 2 studies) | axenic culture |
| *Bipolaris zeicola* | 66 | 0.099 | **0** | unavailable |
| *Puccinia sorghi* | 54 | 0.247 | 4 | effectively unavailable |
| *Pyrenophora tritici-repentis* | 43 | 0.143 | 6 | thin |
| *Bipolaris maydis* | 29 | 0.440 | 74 | axenic culture |
| *Curvularia clavata* | 28 | 0.188 | 27 | axenic culture |
| *Cercospora zeae-maydis* | 24 | 0.419 | **0** | unavailable |
| *Bipolaris sorokiniana* | 21 | 0.103 | 98 | axenic culture |
| *Pyrenophora teres* | 20 | 0.169 | 104 | axenic culture |
| *Fusarium graminearum* | 13 | 0.110 | 2,138 | axenic culture |
| *Fusarium verticillioides* | 14 | 0.126 | 423 | axenic culture |
| *Alternaria alternata* | 42 | 0.068 | 338 | axenic culture |

Three groups fall out of that table and they need different treatment in the text.

**Culturable, abundant, and the core of the set.** *Z. tritici*, *E. turcicum*, *B. maydis*,
*B. sorokiniana*, *P. teres*, *F. graminearum*, *F. verticillioides*, *A. alternata*,
*C. clavata*. All have enough independent studies to take replicates from more than one, which
matters because strain and library prep are the obvious confounders in a single study.

**Obligate biotrophs, where no axenic culture exists.** The four *Puccinia* species and
*Blumeria graminis* cannot be grown without a host, so there is no axenic library to find at
any count. The nearest available thing is a host-free spore library: ungerminated or
germinated urediniospores or conidia with no plant tissue. That is a different control and
must be labelled as one, because a spore transcriptome expresses a narrower gene set than
either culture or infection, and saturation depends on expression breadth. The counts are
large enough that such libraries plausibly exist for *P. striiformis* and *P. triticina*, and
*P. sorghi* has 4 runs in total so it should be written off now rather than searched for.
**Included, decided 2026-10-23**, as tier 3 above: *P. striiformis* alone carries 929 cohort
runs, and what its control buys is not a defence of its own calls, which sit at 0.45
saturation and are never at risk, but the *Puccinia*-into-*Puccinia* smear behind the 96 runs
that collapse from co-infected to single at genus level.

**Absent from the archive.** *B. zeicola* and *C. zeae-maydis* have zero RNA-seq runs
anywhere, which was true in the 2026-10-13 survey and is still true. These are precisely the
two species that head the maize congener clusters, where the real pathogen sits near 42%
coverage and its congeners trail it. **The maize clusters therefore cannot be settled by
controls at all**, and that is an argument for `04_align` rather than a gap to apologise for.
The limitation is stated, not worked around.

## Host controls

| host | ENA RNA-seq | role |
|---|---|---|
| *Zea mays* | 58,814 | per-host false-positive floor |
| *Oryza sativa* | 45,264 | per-host false-positive floor |
| *Triticum aestivum* | 23,773 | per-host false-positive floor |
| *Hordeum vulgare* | 13,695 | per-host false-positive floor |
| *Sorghum bicolor* | 10,769 | per-host false-positive floor |

**Take the paired mock-inoculated arm of a controlled infection wherever one exists**, 2026-10-26.
It holds every variable except the pathogen constant: same plants, same growth conditions, same
library preparation, usually the same flow cell. The independent lab-grown libraries below are
the fallback for studies with no mock arm, not the first choice.

Availability is not a constraint for the fallback, so selection is on growth condition only.
**Lab-grown only**, decided 2026-10-23: axenic or near-axenic seedlings, hydroponic or sterile
substrate. Healthy field-grown plants were considered as an honest comparator for a field cohort
and dropped, because the floor is what is being calibrated and a field plant's endophyte load is
exactly the thing that would blur it.

A lab-grown plant is still not a true zero. No plant is sterile and endophytes exist, so this
gives an upper bound on the false-positive rate, which is the direction that matters: an upper
bound below 0.05 is sufficient to defend the floor.

The host controls also carry the only available test for the failure mode neither pathogen
control catches, which is reads from an organism absent from the database entirely. Classifying
the same host library against db_v3 and db_v4 measures how much of what looks like a pathogen
call is simply the nearest available reference.

## *P. nodorum*, worked through as the example

The 16 ENA runs looked too thin to use and are not. Two are capillary GENOMIC runs mislabelled
as RNA-Seq and are discarded. The other 14 are real transcriptomes from two independent
studies, which satisfies the two-study requirement on its own:

| study | runs | depth | paper | recovered by |
|---|---|---|---|---|
| PRJNA632579 | 8, named `invitro_SN15_*` / `invitro_Pf2_*` | ~4.4 Gkmer | Sci Rep 2019, `10.1038/s41598-019-52444-7` | ENA `XREF_LINK` to PUBMED |
| PRJNA1044467 | 6, SN15 wildtype and *pnvea* mutant | ~25 Gkmer | BMC Microbiol 2024, `10.1186/s12866-024-03454-7` | Europe PMC full-text accession search |

Both are CCDM submissions, so the methods are in-house knowledge rather than a careful read of
someone else's paper. Neither is admitted on its title: PRJNA632579's `invitro_` prefix is the
submitter's word, and PRJNA1044467's six runs still need the culture-versus-in-planta question
settled from the PnVeA methods.

**PRJNA632579 also yields a second control of a different kind.** The study compared in vitro
against in planta, and the eight runs above are only its in vitro arm. The in planta arm is a
known-mixture library, one known pathogen on one known host, in the cohort's own condition.

## The in-house known-mixture libraries

These are controlled single-pathogen infections we already hold, so they belong to the primary
class above rather than to a category of their own.

- **`ada1` / PRJNA1509083.** *P. nodorum* SN15 and SN15 delta-*ada1* on wheat cv. Halberd,
  detached leaf, 3 and 7 dpi, four replicates, poly(A) selected. Deposited, so citable, and
  fetchable from ENA like any cohort run.
- **PRJNA632579's in planta arm**, the counterpart to the eight `invitro_` runs below. One
  study carrying both arms is worth more than either alone: it measures the same organism and
  isolate in culture and in the host, so the expression-breadth difference between the two
  control classes can be quantified rather than assumed.

## The in-house carlia libraries, not needed

`projects/carlia` holds 12 genuinely axenic *P. nodorum* SN15 libraries at 129 to 182 M read
pairs. They were considered because ENA looked nearly empty for this species, and that turned
out to be wrong. They are also not on Setonix scratch any more (592 K remains, no fastq) and
Acacia holds only `pawsey1168-llenzo-kraken-db`, so using them would mean recovering from the
research drive or re-pulling the supplier delivery, and they are undeposited, which is a
reproducibility gap a reviewer can fairly raise. **Dropped, 2026-10-23.** Reconsider only if
the two archive studies fail the methods read.

## Depth, because saturation is not scale-free

Median saturation falls 0.67, 0.23, 0.053, 0.019, 0.002 across redundancy decades, so a
control measured at one depth calibrates the floor at one depth. Each control library is
rarefied with the existing `rarefy.py` across the cohort's range, from the 0.463 Gkmer
inclusion floor to the cohort median near 1.68 Gkmer and above, and the false-positive tail is
read per decade. The floor that comes out is either a constant that holds across the range or
a function of redundancy, and which of those it is should not be assumed in advance.

## Decision rule, written before the measurement

Stated now so the floor is calibrated rather than tuned, which is the same objection that
applies to presenting any floor after the number it produces is known.

1. If the false-positive saturation tail, taken as a high percentile across controls and read
   per depth decade, sits **below 0.05**, the floor stands and the headline stands at 46.8%.
2. If it lands **between 0.05 and 0.15**, the floor moves to it and every figure and number is
   regenerated at the new floor.
3. If it **exceeds 0.15**, or overlaps the distribution of host-resolved pathogen calls so
   that no cut separates them, then the floor cannot be set this way. S1 reports the
   sensitivity range as final and the congener clusters go to `04_align`.

Outcome 3 is a real possibility and is not a failure. It is the honest reading of a method
whose discriminating power has a limit, and the TODO already allows for it.

## Scale and cost

Nine culturable pathogens and up to four biotrophs, at three or four runs each from at least
two independent studies, is 40 to 50 runs. Five hosts lab-grown, three replicates each, is 15.
Two known-mixture studies add about 25. In vitro fungal RNA-seq runs are small, order 1 to 3 GB,
so the whole set is roughly 0.2 TiB against the cohort's remaining 5.4 TiB, and classify is
about 15 minutes of an existing job plus the 41 GB database stage.

The cost is not compute. It is the methods read, which is per STUDY rather than per run and so
lands at roughly 20 to 30 papers, and the single ENA fetch stream, which is saturated by the
cohort until PRJEB83049 and PRJEB67964 are down.

## Layout, and the trap to avoid

Curation lives in `01_literature/controls/` and uses the audit chain already built there, so a
control's axenic status traces to a quoted sentence in a methods section the same way a cohort
value does. Reads and reports go in the existing bronze tree keyed by accession.

**Controls must not enter `run_list.tsv`.** They are not cohort runs and nothing may let them
into a denominator or a figure. They need their own list, and `detect.py`'s inputs should be
explicit about which list they come from rather than picking up whatever is in the reports
tree.
