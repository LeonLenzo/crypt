# Which controls to fetch and classify first

Decided 2026-10-26, Leon's call: go after the pathogens that carry the most weight in the
cohort, which are the obligate biotrophs and *Zymoseptoria*, not the easiest study to fetch.
Two reasons they are highest value. They head the undeclared-detection list, so the floor
most needs to be right for them. And the rusts and mildew are obligate biotrophs for which a
controlled infection is the ONLY possible control, so there is no axenic fallback if this is
skipped.

## Priority order, by undeclared cohort load against control readiness

Undeclared-run counts are from the 1,999 classified runs (coinfection.tsv), 2026-10-26.

    Zymoseptoria tritici   343 undeclared   Child 2023    PRJEB58154  18  poly(A) 150PE  READY
    Puccinia triticina     204              Wei 2020      PRJNA613154  2  mRNA           READY (thin)
    Blumeria graminis      166              Hu 2018       SRP117269   37  poly(A) 150PE  READY
    Exserohilum turcicum    76              Human 2020    PRJNA560644 12  90PE           READY
    Parastagonospora nod.   76              Jones 2019    PRJNA632579 16  (own deposit)  READY
    Bipolaris zeicola       68              --            NONE             zero archive  BLOCKED
    Puccinia graminis       28              Henningsen21  PRJNA483957 54  mock 2/4/6dpi  READY, best
    Puccinia striiformis    26              Dobon 2016    PRJEB12497  42  poly(A)        READY

Pst and Pgt are LOW on undeclared load because they are the declared target of most wheat
projects, so their calls are mostly expected. They still matter as controls because they are
the two highest-coverage species in the cohort and the *Puccinia*-into-*Puccinia* smear behind
the 96 genus-collapse runs is theirs. Keep them in the first batch.

## The first fetch batch

Eight studies, ~160 runs total, all small in vitro / infected-leaf libraries, order 0.1 to
0.3 TiB. Chemistry note: Child (Z. tritici), Hu (Blumeria), Dobon (Pst) and Kim (P. sorghi)
are poly(A) and match the cohort; the rest vary and are read with that caveat.

    PRJEB58154   Zymoseptoria tritici    Child 2023      use IPO323 wildtype runs
    PRJNA613154  Puccinia triticina      Wei 2020        2 runs only
    SRP117269    Blumeria graminis       Hu 2018
    PRJNA483957  Puccinia graminis       Henningsen 2021 has the mock host arm
    PRJEB12497   Puccinia striiformis    Dobon 2016
    PRJNA560644  Exserohilum turcicum    Human 2020
    PRJNA632579  Parastagonospora nod.   Jones 2019      own CCDM deposit
    PRJEB29046   Pyrenophora tr-rep.     Moolhuijzen 2018 the complete-template, keep for the mock+axenic arms

## The two gaps closed today

The first pass picked the wrong paper for the two highest-value species: Rudd 2015 (Z. tritici)
and Bruce 2014 (P. triticina) are strong studies that deposited NO raw reads. Both species have
hundreds of archive runs, so a deposited controlled infection existed; found 2026-10-26:
- Z. tritici: Child 2023 PRJEB58154, ENA, poly(A) 150 PE, IPO323 wildtype plus two signalling
  mutants (all single-strain). No mock arm; Rudd 2015 has the mock time course on the same
  strain/cultivar/lab and is worth chasing for the host floor.
- P. triticina: Wei 2020 PRJNA613154, 2 runs, single-pathotype Thatcher infection. Thin but
  retrievable. Neugebauer 2018 deposited only GenBank TSA contigs.

## Not in the first batch, and why

- *B. zeicola* (68 undeclared) is BLOCKED: zero RNA-seq in any archive, no controlled infection
  published. See ABSENCE_EVIDENCE.md. It is the highest-load species with no possible control.
- *B. maydis*, *F. graminearum*, *F. verticillioides*: lower undeclared load, culturable, and
  can wait for a second batch. *B. maydis* (Yu 2024, growth chamber + axenic) is ready when
  wanted.
- *Curvularia clavata*, *Alternaria* spp.: no cereal-host controlled infection exists, see the
  gap rows in control_candidates.tsv.
