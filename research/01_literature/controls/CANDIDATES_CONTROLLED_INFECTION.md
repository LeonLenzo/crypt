# Controlled-infection control candidates, verified against the archive

From the Undermind search "single pathogen controlled infection transcriptomes cereal fungi",
2026-10-26, with each accession checked to resolve in ENA before listing it here. The search
ranks by topical relevance; this file re-ranks by what a control needs, which is a retrievable
accession and a clean single-pathogen design. These are shortlisted, not admitted: admission
still needs the purity screen on the classified reads (host clade present, no OTHER pathogen
clade), per DESIGN.md.

## The complete template: tan spot

**Moolhuijzen et al. 2018** (BMC Res Notes, `10.1186/s13104-018-3993-2`), ENA **PRJEB29046**,
6 runs, all verified present. One BioProject carries every arm a control needs:

The six runs are in `control_template_runs.tsv`: two mock wheat (ERR2822756-7),
two infected wheat (ERR2822760-1), and two axenic M4 cultures (ERR2822758-9), all 3/4 dpi.

Single isolate (M4, race 1), susceptible wheat cv. Machete, mock arm time-matched at both
timepoints, and the SAME isolate in axenic culture. **ribo-depleted, not poly(A)**, un-stranded
100 bp paired. The ribo-depletion is a caveat: the cohort is poly(A), so this library retains
rRNA and ncRNA the cohort libraries drop, which changes which k-mers are present. Usable for the
cross-mapping question, read with that difference in mind for the floor.

## Obligate biotrophs, the cases that most needed covering

- **Puccinia striiformis** -- Dobon et al. 2016 (BMC Genomics, `10.1186/s12864-016-2684-4`),
  ENA **PRJEB12497**, 42 runs verified. Controlled-environment rooms, single isolate 87/66,
  poly(A), 100 bp paired, a dense time course (0 to 11 dpi) on susceptible Vuka and resistant
  Avocet-Yr5, plus germinating spores as a separate pathogen reference. **poly(A), so chemistry
  matches the cohort.** Control arm is 0 dpi uninoculated, NOT a time-matched mock; good host
  baseline, weaker as a pathogen-free control at the infection timepoints.
- **Blumeria graminis** f. sp. tritici -- Hu et al. 2018 (PLoS ONE, `10.1371/journal.pone.0198891`),
  SRA **SRP117269**, 37 runs in ENA. Growth chamber (16/8, 18 C, 80% RH), single race 15,
  poly(A), 150 bp paired, 0 h baseline then 6/12/24/48/72 hpi on wheat sister lines L658/L958.
  Chemistry and read length are the closest match to the cohort of any candidate. Control is the
  0 h uninoculated baseline only.
- **Puccinia graminis**, **Puccinia triticina** -- multiple studies surfaced (Poudel 2019 on
  barley, Bruce 2014, Rutter 2017, Henningsen 2021); accessions not yet pulled. Both species are
  well represented, so biotroph coverage is not a gap.

## Necrotrophs and hemibiotrophs with controlled-infection RNA-seq

- **Exserohilum turcicum** -- Human et al. 2020 (Front Microbiol, `10.3389/fmicb.2020.00360`),
  SRA **PRJNA560644**, 12 runs in ENA. Growth chamber, single isolate per race (13N, 23N on
  separate plants), 0/2/5/7/13 dpi on maize Va26, plus an in vitro 13N culture. Library selection
  (poly(A) vs rRNA) not stated in the paper, 90 bp strand-specific paired. Control is 0 dpi
  pre-inoculation only, no time-matched mock. Note: SIX13/SIX5-like Sanger sequences share the
  BioProject, so filter to the RNA-seq runs.
- **Bipolaris maydis** -- Meshram et al. 2022 (Front Microbiol). Accession not yet pulled.
- **Bipolaris sorokiniana** -- Basak et al. 2024 (Front Microbiol). Accession not yet pulled.
- **Pyrenophora teres** f. teres -- Moolhuijzen et al. 2021 (MPMI), same group as the tan spot
  template, so the design is likely to carry a mock arm; worth pulling for that reason.
- **Zymoseptoria tritici** -- Rudd et al. 2015 (Plant Physiol), the canonical Z. tritici dual
  RNA-seq time course. Accession not yet pulled.
- **Fusarium graminearum** -- Pan 2018, Rocher 2022, and others; abundant.

## The full candidate list lives in `control_candidates.tsv`

One row per target pathogen from `ena_availability_2026-10-23.tsv`. Every accession was checked
to resolve in both ENA and SRA on 2026-10-26; the `ena_runs`/`sra_runs` columns are what each
returns. The `tier` column ranks by what a control needs, not by topical relevance:

- **1-complete / 1-mock+live** -- live accession AND a pathogen-free host arm in the same
  deposit. The ready-to-use set: Ptr (the complete template, with mock + infection + axenic),
  *P. sorghi*, *B. sorokiniana*, *P. graminis*, *P. nodorum*.
- **2-baseline only** -- live accession, but the only pathogen-free sample is a 0 h/0 dpi
  pre-inoculation point: *P. striiformis*, *B. graminis*, *E. turcicum*.
- **3-no host floor** -- live accession, no pathogen-free host sequenced: *B. maydis*,
  *F. graminearum*. Still give the cross-mapping spectrum.
- **4-no accession** -- strong study, no repository named. *Z. tritici* (Rudd 2015) is worth
  chasing: it has the mock time course but publishes no accession.
- **5-not pulled** -- *P. teres* (no downloadable PDF), *F. verticillioides* (not yet searched).
- **6-gap** -- no cereal-host controlled infection exists: *C. clavata*, *A. tenuissima*,
  *A. alternata* (its only study used citrus).
- **7-zero archive** -- *B. zeicola*, *C. zeae-maydis*; see ABSENCE_EVIDENCE.md.

Chemistry matters for the floor: *P. striiformis*, *B. graminis*, *P. sorghi* are poly(A) and
match the cohort; Ptr and *P. triticina* are ribo-depleted, which retains rRNA/ncRNA the cohort
drops and must be read with that difference in mind.

## Gaps the search confirmed

- **No cereal-host controlled infection surfaced for *Curvularia clavata* or *Alternaria
  tenuissima*.** The *A. alternata* infection study it found used CITRUS, not a cereal.
- **The four zero-archive species** (*B. zeicola*, *C. zeae-maydis*, *B. victoriae*,
  *F. sporotrichioides*) were the first search's subject and remain uncovered; see
  [ABSENCE_EVIDENCE.md](ABSENCE_EVIDENCE.md). *C. zeae-maydis* has in planta RNA-seq published
  but behind a non-resolving accession.

## What a reviewer gets from this

The two things most likely to be challenged are now answered from retrievable data: the
obligate biotrophs have real single-pathogen infection libraries (so the saturation floor can be
calibrated for the rusts and mildew), and at least one study carries the mock arm, the infection
arm and the axenic arm of the same isolate in one deposit, which lets the expression-breadth
difference between culture and infection be measured rather than assumed.

## The archive mislabels one run of the Ptr template, and the paper settles it

Found 2026-10-26 while fetching. ENA's `sample_title` for **ERR2822761** reads "4 days post
control inoculation in planta", which is a verbatim copy of ERR2822757's title. Moolhuijzen
2018's own data-availability section names ERR2822761 **Sample P2, the 4 DPI
*P. tritici-repentis* infection** (Data set 10), and ERR2822757 Sample C2. Taking ENA at its
word would make the deposit three mock arms and one infection rather than two and two, and
would put a mock label on the only 4 dpi infection in the template.

The roles in `control_template_runs.tsv` follow the PAPER and are correct. Classification
settles it independently and without appeal to either source: a mock library carries no Ptr.
Check it when PRJEB29046 lands. The ENA capture is kept at
`02_classify/data/bronze/controls/meta/PRJEB29046.tsv` so the disagreement stays visible
rather than being corrected away.
