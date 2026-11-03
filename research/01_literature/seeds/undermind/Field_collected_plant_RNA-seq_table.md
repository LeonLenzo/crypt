# Field collected plant RNA-seq table

##### [**Undermind**](https://undermind.ai)

---


## Table of Contents

- [Field-collected plant RNA-seq studies](#field-collected-plant-rna-seq-studies)
  - [Eligible study designs](#eligible-study-designs)
  - [Archive and site/date metadata](#archive-and-sitedate-metadata)
  - [Pathogens and diseases reported by the authors](#pathogens-and-diseases-reported-by-the-authors)
  - [Secondary reanalysis relevant to unreported pathogen signal](#secondary-reanalysis-relevant-to-unreported-pathogen-signal)
  - [Unresolved or excluded candidates](#unresolved-or-excluded-candidates)
  - [References](#references)

# Field-collected plant RNA-seq studies

**Screening status: partial.** Thirty-two studies are currently classified as eligible from the papers and their reported raw-read deposits. This is not a complete census: the broad search returned 208 candidates and missed clearly eligible studies. Only a subset of high-likelihood field/aerial records from the earlier DOI batch has been checked; the full batch still needs reconciliation. Accessions and sample metadata have not yet all been independently verified run by run. “Not stated” is not inferred from a standard protocol.

## Eligible study designs

| Study | Host and eligible aerial tissue | Collection country/region and years | RNA-seq sample count | Library preparation and RNA selection |
|:---|:---|:---|---:|:---|
| \[Ada21\] | Bread wheat, infected leaves | 30 countries across wheat-growing regions; 2014–2018 | 538 new field libraries; the paper also reprocessed 486 older datasets | Total RNA extraction; TruSeq RNA Sample Preparation Kit; polyA/rRNA selection not specified |
| \[Hub15\] | Wheat and triticale, infected leaves | United Kingdom, 17 counties; spring–summer 2013 | 39 field libraries: 35 wheat and 4 triticale | Total RNA; TruSeq RNA Sample Preparation Kit; selection not specified |
| \[Bue17\] | Wheat, triticale, and rye, infected leaves | 16 European countries, plus current samples from Pakistan, Ethiopia, Chile, and New Zealand; December 2013–August 2014 | 133 new field RNA-seq samples (115 European, 18 additional global samples); 246 European samples were collected | Total RNA; TruSeq RNA Sample Preparation Kit; selection not specified |
| \[Dia19\] | Soybean, symptomatic leaves | Canada, Manitoba; V2/3 and R6 surveys in 2016 | 25 pooled RNA-seq libraries representing over 650 leaves from 81 fields; broader survey: 656 leaves from 103 fields | RNA-seq cDNA libraries; polyA/rRNA selection not stated in accessible methods |
| \[Mor17\] | Grapevine, woody trunk tissue at lesion margins | United States, Sonoma and Fresno Counties, California; 2014–2015 | 28 naturally infected field libraries; 12 controlled-inoculation libraries are not counted | Total RNA; Illumina TruSeq; selection not specified |
| \[Riv22\] | Tomato and weed species, leaves and fruit | Slovenia, 14 farms in six localities; summers 2019–2020 | 67 composite field libraries from 436 plants (293 tomato, 143 weeds) | Total RNA; Ribo-Zero Plant rRNA depletion; TruSeq RNA libraries |
| \[Cha22c\] | 23 weed species in rice fields, leaves and stems | China, Hangzhou, Xiaogan, Guiyang, and Sanya; 2018–2019 | 29 samples/libraries, one per sampled weed individual | Total RNA; Ribo-Zero Plant Leaf rRNA depletion; TruSeq RNA library prep |
| \[Jo20b\] | Soybean, symptomatic leaves | Republic of Korea, eight provinces; June 2016 | 172 plants represented by 12 libraries: 8 province pools and 4 single-plant libraries | PolyA-selected mRNA; NEBNext Ultra RNA Library Prep |
| \[Kam16\] | *Arabidopsis halleri* subsp. *gemmifera*, leaves | Japan, Omoide River, Hyogo Prefecture; May 27 and June 30, 2014 | 68 individual libraries | Total RNA; RNase H rRNA depletion; not polyA-selected |
| \[Elw23\] | Potato, leaves | Egypt, El Beheira Governorate, 10 farms; winter 2020 | 3 libraries selected from 20 collected plants | Total RNA; Ribo-Zero Plant Leaf rRNA depletion; TruSeq Stranded Total RNA kit |
| \[Por21\] | Bread wheat, flag leaves | Switzerland, field near Zurich; year not stated | 30 libraries: five lines, two field conditions, three biological replicates | TruSeq Stranded mRNA; mRNA selection |
| \[Sat19\] | *Arabidopsis thaliana*, leaves | Switzerland, University of Zurich Irchel field site; August 4, 2016 | 190 plants sequenced/prepared; 173 retained in final expression analysis | Total RNA; selective depletion of rRNA and highly abundant transcripts |
| \[Var19\] | Sorghum, leaves; the study also sequenced roots, which are excluded from the tissue count | United States, Parlier, California; 17-week field time series, calendar year not explicit in the article | 198 leaf libraries (198 additional root libraries) | mRNA-seq; specific selection/depletion step not stated |
| \[Ple15\] | Rice, leaves | Philippines, IRRI experimental fields, Los Baños, Laguna; dry and wet seasons, 2013 | 240 bulk libraries/samples | Total RNA; Ribo-Zero Plant Leaf rRNA depletion; ScriptSeq v2 |
| \[Mar21g\] | 24 vascular plant species, mature leaves | United States, Harvard Forest, Massachusetts; July and August 2016 | 48 samples (24 species sampled twice) | Ovation RNA-Seq System with SPIA amplification; Kapa Illumina library prep |
| \[Mey23c\] | Winter rapeseed, rosette leaf 8 | Belgium, Merelbeke; RNA-seq tissue collected November 28, 2016; field trial continued into 2017 | 62 individual plants | Total RNA and double-stranded cDNA; TruSeq DNA PCR-Free library kit; no polyA/rRNA step reported |
| \[Sez23\] | American beech, leaves | United States, Harvard Forest, Massachusetts, and SERC, Maryland; 2017–2018 | 222 samples from 40 trees | E.Z.N.A. Plant RNA Kit, DNase, TruSeq libraries; selection not stated |
| \[Cru20b\] | Maize, mature ear leaf 16 | Belgium, Zwijnaarde; August 25 and September 2, 2015 | 60 individual plants | Total RNA; purified polyA-containing mRNA; NEBNext library prep |
| \[Har23c\] | Grapevine, leaves and reproductive tissues (flower buds and berries) | United States, Mount Vernon, Missouri; 2017–2019 | 1,178 samples after sequencing/QC across three years | Bulk 3′-RNA-seq; specific polyA/rRNA step not stated in the paper summary |
| \[Sab22\] | Grapevine, ripe berries | Italy, Locorotondo, Apulia; harvest year not reported | 9 libraries: 3 cultivars × 3 replicates | Total RNA; TruSeq Stranded Total RNA kit; selection/depletion not separately specified |
| \[Trz25\] | Barley, field leaves | Poland; 2022–2023 growing season; exact source field not assigned to the pool | 1 pooled library | Total RNA; Ribo-off rRNA depletion; VAHTS Universal V6 RNA-seq library kit |
| \[Nem25c\] | Alfalfa, leaves | United States, five commercial fields in Benton and Yakima Counties, Washington; mid-July 2024 | 10 pools: one asymptomatic and one symptomatic pool per field, each representing five plants | PolyA-selected mRNA; ABclonal Fast RNA-seq Lib Prep Kit v2 |
| \[Man24b\] | Sorghum, mature leaves | United States, Havelock Research Farm, Lincoln, Nebraska; August 5, 2021 | 748 libraries (738 unique genotypes plus replicates/checks) | PolyA-selected mRNA; TruSeq Stranded mRNA |
| \[Tor23d\] | Maize, mature leaves | United States, Havelock Farm, Lincoln, Nebraska; July 8, 2020 | 693 libraries retained from 750 field plants after QC | PolyA-selected mRNA; TruSeq strand-specific RNA-seq |
| \[Li24r\] | Soybean, shoot tissue above the cotyledon node at V2 | China, Sanya; collection year not reported | 622 accessions/libraries | PolyA-selected mRNA; NEBNext Ultra RNA Library Prep |
| \[Tom26\] | *Arabidopsis thaliana*, leaves | Japan, Otsu, and Switzerland, Zurich; early summers 2017–2018 | 2,381 individual field transcriptomes | Bulk Lasy-Seq 3′-tag RNA-seq; specific selection details not stated in the paper summary |
| \[Cha13b\] | *Eutrema salsugineum*, cauline leaves | Canada, field site about 40 km northwest of Whitehorse, Yukon Territory; late June 2005 | 3 in-situ field libraries (YF1–YF3); 6 cabinet libraries are excluded | Total RNA; 3–4 rounds of oligo-dT mRNA selection; random-primed cDNA |
| \[Dan19c\] | Maize, leaves, ears, and tassels | United States, Woodland, California; 2012 | 96 libraries: 3 aerial tissues × 4 stages × 2 watering treatments × 4 replicates | mRNA-seq; exact selection step not stated in the extracted methods |
| \[Zha17e\] | Maize, leaves | China, Urumqi, Xinjiang; field-study year not stated | 16 libraries | PolyA-selected mRNA; NEBNext Ultra RNA Library Prep |
| \[Kas18\] | Rice, youngest fully expanded leaves | Japan, Takatsuki (2015) and Kizugawa (2016) | 1,026 libraries retained from 1,069 collected leaf samples | Total RNA; enzymatic depletion of abundant RNAs/rRNA, not polyA selection |
| \[Cro17\] | Douglas-fir, mature needles | United States, Central Point and Corvallis, Oregon; annual samples Oct 2010–Nov 2011, diurnal series Sep 7–8 (year not explicit) | 179 libraries from 19 trees across annual and diurnal series | Total RNA; strand-specific TruSeq mRNA-seq; polyA selection inferred from protocol, not explicitly named |
| \[Mar26e\] | *Arabidopsis thaliana*, above-ground rosette tissue | Spain, El Castillejo Botanical Garden, Cádiz; four dates in each of 2020–2021 and 2021–2022, morning and afternoon | GEO/SRA lists 80 samples; article methods report 84 per year (168 total), so the library count needs reconciliation | TruSeq RNA Sample Prep v2; selection/depletion not stated. |

## Archive and site/date metadata

| Study | Public raw-read accession reported | Accession-to-site-and-date metadata in paper or supplement? | Notes |
|:---|:---|:---|:---|
| \[Ada21\] | ENA PRJEB39201, PRJEB33109, PRJEB31334, PRJEB15280, PRJEB12497; SRA PRJNA256347, PRJNA181960, PRJNA176472, PRJNA396589 | Yes | Supplementary Table S1 includes collection location and date; project list spans new and reused data. |
| \[Hub15\] | SRA PRJNA256347 and PRJNA257181 | Partial | Supplementary Table S1 links sample to county/host; date is given as spring–summer 2013, not established per accession. |
| \[Bue17\] | ENA PRJEB15280 | Yes | Supplementary Tables S1 and S3 give sample-to-site/date data for European and additional global samples. |
| \[Dia19\] | SRA SRP105188; runs SRR5480332–SRR5480333 | Partial | Sequencing libraries pool many field leaves; the article/supplement reports survey-stage/field metadata but does not map individual leaves to a run. |
| \[Mor17\] | SRA PRJNA352065 / SRP092409 | Partial | Appendix S1 Table S3 links field samples to county, cultivar, symptoms, and year; exact collection day is not established here. |
| \[Riv22\] | SRA PRJNA772045 | Partial | Supplementary Table 2 links composite library to locality, year, plant type, and health status; no exact date per individual plant. |
| \[Cha22c\] | SRA PRJNA670260; run accessions listed in paper | Partial | Supplementary Table S1 gives host/location; exact accession-to-day linkage not confirmed. |
| \[Jo20b\] | SRA PRJNA637168 | Partial | Table 3 maps libraries to province/pool or individual plant; June 2016 is given for survey, not exact date per accession. |
| \[Kam16\] | DDBJ DRA003823 | Partial | Supporting Table S1 lists 68 individuals; paper describes two collection dates and one field site. Accession-level date crosswalk still needs checking. |
| \[Elw23\] | SRA PRJNA743866 | Partial | Table 1 maps sequenced sample IDs to farm/land; dates are given as winter 2020. |
| \[Por21\] | SRA PRJNA718488 | No | Paper gives a field near Zurich and adult-plant sampling, but no calendar date or accession-to-site/date table. |
| \[Sat19\] | SRA PRJNA488315 | Partial | Paper reports one field site and August 4, 2016; no accession-to-individual crosswalk confirmed. |
| \[Var19\] | GEO GSE128441; associated SRA project PRJNA527782 / SRP188707 | Partial | Public records include leaf/root and time-series metadata; article describes weekly sampling, but calendar-date linkage for every leaf accession was not established. |
| \[Ple15\] | GEO GSE73609; SRA BioProject PRJNA297424 | Yes | GEO/SRA sample records link leaf libraries to date, season, field type, replicate, genotype, and time point. |
| \[Mar21g\] | SRA PRJNA422719 / SRP127805 | Yes | Table 1 and Appendix 1 provide per-sample accession, locality, coordinates, and date. |
| \[Mey23c\] | ArrayExpress E-MTAB-11904 | Yes | Supporting data link plant-level metadata and field sampling. |
| \[Sez23\] | SRA PRJNA630305 | Partial | Sample naming encodes site/year/season/tree ID, but not exact collection date. |
| \[Cru20b\] | ArrayExpress E-MTAB-8944 | Yes | Field position and harvest date are recorded for sampled plants. |
| \[Har23c\] | SRA PRJNA674915 and PRJNA915033 | Partial | Sample-level metadata are in an external Figshare record; the article gives site/years but the accession-to-site/date mapping is not in paper/supplement. |
| \[Sab22\] | SRA PRJNA799026 | Partial | Paper links libraries to cultivar but does not report harvest year or a complete accession-to-date table. |
| \[Trz25\] | SRA PRJNA1232637; run SRR35566536 | Partial | One pool is described as Polish field barley from 2022–2023; exact field/date mapping is absent. |
| \[Nem25c\] | SRA PRJNA1191001 | Yes | Paper and Supplementary Table S1 link each pool to field number, health status, and the July 2024 collection window. |
| \[Man24b\] | ENA PRJEB83049 | Yes | Samples share one Nebraska field and August 5, 2021 collection date; genotype/sample metadata are reported. |
| \[Tor23d\] | ENA PRJEB67964 | Yes | Supplementary Table S2 provides sample/genotype and field-position metadata; all were collected July 8, 2020. |
| \[Li24r\] | GSA CRA009979 | Partial | Paper gives Sanya field and a shared 9–11 a.m. sampling window, but not a per-accession calendar date. |
| \[Tom26\] | SRA PRJNA1055060, PRJNA1055104, PRJNA1055317, PRJNA1055424, PRJNA1055734, PRJNA1055736, PRJNA1056126, PRJNA1055939 | Partial | Paper gives site/year/plant identity and harvest window; detailed phenotypic metadata are in external repositories, not a paper supplement. |
| \[Cha13b\] | GEO GSE49378; linked NCBI BioProject PRJNA213746 | Partial | Three field libraries are identified as YF1–YF3 from Yukon in late June 2005; individual plants are not linked to separate dates. |
| \[Dan19c\] | GEO GSE71723; SRA SRP062027 / BioProject PRJNA291919 | Yes | GEO sample records link tissue, stage, treatment, replicate, and date at the Woodland field site. |
| \[Zha17e\] | SRA SRP102142 | Partial | Paper maps all samples to Urumqi and relative drought time points; exact calendar dates are not reported. |
| \[Kas18\] | DDBJ PRJDB7234 | Yes | Supplementary Table S2 links sample IDs to transplant set and sampling time; paper gives field sites and years. |
| \[Cro17\] | SRA SRP018395; GEO GSE44058 | Partial | Supplementary files link samples to sampling intervals and environmental data; the diurnal series’ collection year still needs confirmation. |
| \[Mar26e\] | GEO GSE301981; linked SRA BioProject PRJNA1287953 | Yes, with count discrepancy | GEO lists 80 SRA-linked samples; article methods report 84 samples per year. Table 2 lists dates/developmental stages; verify which samples have raw reads. |

## Pathogens and diseases reported by the authors

| Study | Pathogen/disease focus or reported detections | Multiple pathogens reported in one sample? |
|:---|:---|:---|
| \[Ada21\] | Wheat yellow/stripe rust, *Puccinia striiformis* f. sp. *tritici* (Pst). | No multi-pathogen detection reported; host and Pst transcripts are not two pathogens. |
| \[Hub15\] | Pst field pathogenomics. | No; Pst population/genotype composition only. |
| \[Bue17\] | Pst field pathogenomics. | No other pathogens reported; samples were predominantly a single Pst genotype. |
| \[Dia19\] | Soybean field pathogen survey; bacteria, fungi, viruses, oomycetes, and incidental crop pathogens. | Multiple taxa occurred in sequenced pools; individual-leaf coinfection was not resolved by pooling. |
| \[Mor17\] | Grapevine trunk-disease fungi, including *Eutypa lata*, *Phaeomoniella chlamydospora*, *Diplodia seriata*, *Phaeoacremonium minimum*, and *Diaporthe ampelina*. | Yes; multiple fungi were detected in the same naturally infected wood samples. |
| \[Riv22\] | Tomato/weed virome; 126 viruses, including known and novel viruses. | Yes; co-detections were reported in plant samples and composite pools. |
| \[Cha22c\] | RNA viruses from weeds in rice fields; 224 viruses identified. | Yes; one *Leptochloa chinensis* sample contained up to 99 viruses. |
| \[Jo20b\] | Soybean virome; ten RNA viruses, including SMV, SYMMV, SYCMV, PeMoV, PSV, TSWV, BCMV, BCMNV, CMV, and WVMV. | Yes; individual-plant co-infection was frequent, and all 12 libraries contained multiple virus types. |
| \[Kam16\] | Natural *A. halleri* virus community: TuMV, CMV, BrYV, and novel *A. halleri* partitivirus 1. | Yes; multiple infections were common among individual plants/leaves. |
| \[Elw23\] | Potato virus survey; AMV, PLRV, PVY and additional viral/microbial sequences. | Yes; all three sequenced samples were reported to contain multiple viruses. Some additional hits may be endogenous or ambiguous. |
| \[Por21\] | Wheat leaf rust (*Puccinia triticina*) and powdery mildew (*Blumeria graminis* f. sp. *tritici*). | Yes; mixed infection by both fungi was reported in infected plots. |
| \[Har23c\] | Not pathogen-focused; rootstock/environment effects on grapevine expression. Powdery mildew is noted as affecting 2019 berry harvest, not as an RNA-seq detection. | No pathogen co-detection reported. |
| \[Sab22\] | Not pathogen-focused; grape berry phenolics, antioxidant properties, and cultivar transcriptomes. | No pathogen co-detection reported. |
| \[Trz25\] | Barley virus G targeted in a field barley pool; other viruses were also detected. | Yes; BVG, BYDV-PAV, BYDV-PAS, WDV, and WSMV were reported in the one pool. |
| \[Nem25c\] | Field alfalfa pathobiome; viruses, bacteria, and fungi, including PeSV, AMV, SRAV, BLRV, bacterial pathogens, and several fungal pathogen genera. | Yes; multiple potential pathogens in pooled samples; individual plants were not sequenced separately. |
| \[Man24b\] | Not pathogen-focused; flowering-time expression in sorghum. | No pathogen detections reported. |
| \[Tor23d\] | Not pathogen-focused; maize expression and flowering-time traits. | No pathogen detections reported. |
| \[Li24r\] | Not pathogen-focused; baseline soybean V2 shoot expression for TWAS, without biotic treatment. | No pathogen detections reported. |
| \[Tom26\] | Field herbivory and SA/JA biomarkers in Arabidopsis; studied insects, not plant pathogens. | No pathogen co-detection reported. Insects are not pathogens. |
| \[Cha13b\] | Not pathogen-focused; natural-habitat versus cabinet expression in *Eutrema*. | No pathogen detections reported. |
| \[Dan19c\] | Not pathogen-focused; field drought/heat responses in maize. | No pathogen detections reported. |
| \[Zha17e\] | Not pathogen-focused; field drought responses in maize. | No pathogen detections reported. |
| \[Kas18\] | Not pathogen-focused; field rice expression eQTLs; defense genes are discussed, but pathogens were not surveyed. | No pathogen co-detection reported. |
| \[Sat19\] | Herbivore-related field expression in Arabidopsis; mustard aphids and flea beetles were studied. | No pathogen co-detection reported. Insects are not pathogens. |
| \[Var19\] | Drought responses in sorghum; the paper discusses arbuscular mycorrhizal symbiosis, not a pathogen survey. | No pathogen co-detection reported. |
| \[Ple15\] | Rice expression under field weather and water conditions. | Not pathogen-focused; no pathogen co-detection reported. |
| \[Mar21g\] | Plant-community RNA-seq pilot; apparently healthy plants. | Not pathogen-focused; no pathogen co-detection reported. |
| \[Mey23c\] | Rapeseed yield prediction from leaf expression; bacterial/herbivore-related functional signals are not pathogen detections. | No pathogen co-detection reported. |
| \[Sez23\] | Seasonal beech leaf expression; beech bark disease is mentioned as site context, not as a sequenced pathogen target. | No pathogen co-detection reported. |
| \[Cru20b\] | Field maize gene-function/trait study; possible variable bacterial/herbivore exposure is discussed, not confirmed pathogen detection. | No pathogen co-detection reported. |
| \[Cro17\] | Not pathogen-focused; analyzed seasonal and diel expression in Douglas-fir needles. The paper notes non-host/contaminant transcripts but no pathogen survey. | No pathogen co-detection reported. |
| \[Mar26e\] | Not pathogen-focused; studied daily, seasonal, and annual expression plasticity in outdoor common gardens. Biotic-stress gene ontology terms are discussed, but no pathogens are identified. | No pathogen co-detection reported. |

## Secondary reanalysis relevant to unreported pathogen signal

\[Col25b\] reanalyzed existing polyA-enriched bulk RNA-seq from field cohorts rather than generating a new cohort. The reanalysis covered 2,194 sorghum, maize, and soybean samples across ENA PRJEB83049, ENA PRJEB67964, and GSA CRA009979. The original data-generating studies \[Man24b\], \[Tor23d\], and \[Li24r\] are listed separately above and should be counted as the sampling studies. The reanalysis recovered hundreds of fungal taxa, including taxa discussed as plant pathogens, from standard host RNA-seq and reported multi-taxon fungal communities. It is directly relevant to the unreported-pathogen-signal question but is not a separate set of field collections.

## Unresolved or excluded candidates

- \[Har24\] has a public SRA BioProject, but the paper identifies the source only as a botanic garden. It is not counted until outdoor versus glasshouse provenance is confirmed for the sequenced samples.
- \[Wam18\] reports PRJNA42371, but the NCBI project record resolves to an unrelated bacterial genome project; the correct plant raw-read accession is unresolved.
- \[Gaa20\], \[Fow21\], \[Mut18\], and \[Xia22e\] report assembled viral GenBank accessions but no raw-read archive accession in the accessible paper/data statements.
- \[Hod20\] and \[Dia23\] have relevant field RNA-seq designs, but valid raw-read archive accessions have not yet been verified.
- \[Kam18\] is a relevant natural-community RNA-seq survey, but its public raw-read accession has not yet been located.
- \[Lia21\] reports supporting information rather than a public raw-read archive; not counted.
- \[Mcc26\] sequenced wheat tissue only after field rust isolates were propagated under controlled conditions; not eligible as field-collected RNA.
- \[Ric12\] and \[San13c\] are field studies but used microarrays, not bulk RNA-seq; not eligible.

---

## References

\[Ada21\] T. Adams *et al.*, “Rust expression browser: an open source database for simultaneous analysis of host and pathogen gene expression profiles with expVIP,” *BMC Genomics*, vol. 22, Mar. 2021, doi: [10.1186/s12864-021-07488-3](https://doi.org/10.1186/s12864-021-07488-3).

\[Hub15\] A. Hubbard *et al.*, “Field pathogenomics reveals the emergence of a diverse wheat yellow rust population,” *Genome Biology*, vol. 16, Feb. 2015, doi: [10.1186/s13059-015-0590-8](https://doi.org/10.1186/s13059-015-0590-8).

\[Bue17\] V. Bueno-Sancho *et al.*, “Pathogenomic Analysis of Wheat Yellow Rust Lineages Detects Seasonal Variation and Host Specificity,” *Genome Biology and Evolution*, vol. 9, pp. 3282–3296, Nov. 2017, doi: [10.1093/gbe/evx241](https://doi.org/10.1093/gbe/evx241).

\[Dia19\] G. Díaz-Cruz, C. M. Smith, K. F. Wiebe, S. M. P. Villanueva, A. Klonowski, and B. J. Cassone, “Applications of Next-Generation Sequencing for Large-Scale Pathogen Diagnoses in Soybean.” *Plant disease*, vol. 103 6, pp. 1075–1083, Apr. 2019, doi: [10.1094/PDIS-05-18-0905-RE](https://doi.org/10.1094/PDIS-05-18-0905-RE).

\[Mor17\] A. Morales-Cruz *et al.*, “Closed-reference metatranscriptomics enables in planta profiling of putative virulence activities in the grapevine trunk-disease complex,” *bioRxiv*, Jan. 2017, doi: [10.1101/099275](https://doi.org/10.1101/099275).

\[Riv22\] M. P. Rivarez *et al.*, “In-depth study of tomato and weed viromes reveals undiscovered plant virus diversity in an agroecosystem,” Jul. 02, 2022. doi: [10.1101/2022.06.30.498278](https://doi.org/10.1101/2022.06.30.498278).

\[Cha22c\] S.-F. Chao *et al.*, “Novel RNA Viruses Discovered in Weeds in Rice Fields,” *Viruses*, vol. 14, Nov. 2022, doi: [10.3390/v14112489](https://doi.org/10.3390/v14112489).

\[Jo20b\] Y. Jo *et al.*, “Soybean Viromes in the Republic of Korea Revealed by RT-PCR and Next-Generation Sequencing,” *Microorganisms*, vol. 8, Nov. 2020, doi: [10.3390/microorganisms8111777](https://doi.org/10.3390/microorganisms8111777).

\[Kam16\] M. Kamitani, A. Nagano, M. Honjo, and H. Kudoh, “RNA-Seq reveals virus–virus and virus–plant interactions in nature,” *FEMS Microbiology Ecology*, vol. 92, Aug. 2016, doi: [10.1093/femsec/fiw176](https://doi.org/10.1093/femsec/fiw176).

\[Elw23\] E. Elwan, M. Rabie, E. A. Aleem, F. Fattouh, M. S. Kagda, and H. A. H. Zaghloul, “Exploring virus presence in field-collected potato leaf samples using RNA sequencing,” *Journal of Genetic Engineering & Biotechnology*, vol. 21, Oct. 2023, doi: [10.1186/s43141-023-00561-2](https://doi.org/10.1186/s43141-023-00561-2).

\[Por21\] M. Poretti *et al.*, “Comparative Transcriptome Analysis of Wheat Lines in the Field Reveals Multiple Essential Biochemical Pathways Suppressed by Obligate Pathogens,” *Frontiers in Plant Science*, vol. 12, Sep. 2021, doi: [10.3389/fpls.2021.720462](https://doi.org/10.3389/fpls.2021.720462).

\[Sat19\] Y. Sato *et al.*, “Transcriptional Variation in Glucosinolate Biosynthetic Genes and Inducible Responses to Aphid Herbivory on Field-Grown Arabidopsis thaliana,” *Frontiers in Genetics*, vol. 10, Feb. 2019, doi: [10.3389/fgene.2019.00787](https://doi.org/10.3389/fgene.2019.00787).

\[Var19\] N. Varoquaux *et al.*, “Transcriptomic analysis of field-droughted sorghum from seedling to maturity reveals biotic and metabolic responses,” *Proceedings of the National Academy of Sciences of the United States of America*, vol. 116, pp. 27124–27132, Dec. 2019, doi: [10.1073/pnas.1907500116](https://doi.org/10.1073/pnas.1907500116).

\[Ple15\] A. Plessis *et al.*, “Multiple abiotic stimuli are integrated in the regulation of rice gene expression under field conditions,” *eLife*, vol. 4, Nov. 2015, doi: [10.7554/eLife.08411](https://doi.org/10.7554/eLife.08411).

\[Mar21g\] H. E. Marx, S. A. Jorgensen, E. Wisely, Z. Li, K. M. Dlugosch, and M. S. Barker, “Pilot RNA‐seq data from 24 species of vascular plants at Harvard Forest,” *Applications in Plant Sciences*, vol. 9, Feb. 2021, doi: [10.1002/aps3.11409](https://doi.org/10.1002/aps3.11409).

\[Mey23c\] S. D. Meyer *et al.*, “Predicting yield of individual field-grown rapeseed plants from rosette-stage leaf gene expression,” *PLOS Computational Biology*, vol. 19, May 2023, doi: [10.1371/journal.pcbi.1011161](https://doi.org/10.1371/journal.pcbi.1011161).

\[Sez23\] U. U. Sezen, J. Shue, S. J. Worthy, S. Davies, S. McMahon, and N. Swenson, “Leaf gene expression trajectories during the growing season are consistent between sites and years in American beech,” *bioRxiv*, Nov. 2023, doi: [10.1098/rspb.2023.2338](https://doi.org/10.1098/rspb.2023.2338).

\[Cru20b\] D. Cruz *et al.*, “Using single‐plant‐omics in the field to link maize genes to functions and phenotypes,” *Molecular Systems Biology*, vol. 16, Apr. 2020, doi: [10.15252/msb.20209667](https://doi.org/10.15252/msb.20209667).

\[Har23c\] Z. N. Harris *et al.*, “Grapevine scion gene expression is driven by rootstock and environment interaction,” *BMC Plant Biology*, vol. 23, Jan. 2023, doi: [10.1186/s12870-023-04223-w](https://doi.org/10.1186/s12870-023-04223-w).

\[Sab22\] W. Sabetta *et al.*, “‘Good Wine Makes Good Blood’: An Integrated Approach to Characterize Autochthonous Apulian Grapevines as Promising Candidates for Healthy Wines,” *International Journal of Biological Sciences*, vol. 18, pp. 2851–2866, Apr. 2022, doi: [10.7150/ijbs.70287](https://doi.org/10.7150/ijbs.70287).

\[Trz25\] K. Trzmiel, A. Zarzyńska‐Nowak, and B. Hasiów‐Jaroszewska, “Occurrence and Molecular Characteristics of Polerovirus BVG Isolates from Poland,” *Pathogens*, vol. 14, Oct. 2025, doi: [10.3390/pathogens14111087](https://doi.org/10.3390/pathogens14111087).

\[Nem25c\] L. Nemchinov, B. Irish, S. Grinstead, and O. A. Postnikova, “Alfalfa transcriptomic responses to the field pathobiome,” *Plant Biology (Stuttgart, Germany)*, vol. 27, pp. 492–503, Jan. 2025, doi: [10.1111/plb.70021](https://doi.org/10.1111/plb.70021).

\[Man24b\] H. Mangal *et al.*, “Genes and pathways determining flowering time variation in temperate‐adapted sorghum,” *The Plant Journal*, vol. 122, Dec. 2024, doi: [10.1111/tpj.70250](https://doi.org/10.1111/tpj.70250).

\[Tor23d\] J. V. Torres-Rodríguez *et al.*, “Population level gene expression can repeatedly link genes to functions in maize,” *bioRxiv*, Nov. 2023, doi: [10.1101/2023.10.31.565032](https://doi.org/10.1101/2023.10.31.565032).

\[Li24r\] D.-L. Li *et al.*, “TWAS facilitates gene-scale trait genetic dissection through gene expression, structural variations, and alternative splicing in soybean,” *Plant Communications*, vol. 5, Jun. 2024, doi: [10.1016/j.xplc.2024.101010](https://doi.org/10.1016/j.xplc.2024.101010).

\[Tom26\] A. Tomita *et al.*, “Transcriptomic biomarkers reveal jasmonic and salicylic acid state under field herbivory,” Feb. 12, 2026. doi: [10.1101/2025.05.29.656841](https://doi.org/10.1101/2025.05.29.656841).

\[Cha13b\] M. Champigny *et al.*, “RNA-Seq effectively monitors gene expression in Eutrema salsugineum plants growing in an extreme natural habitat and in controlled growth cabinet conditions,” *BMC Genomics*, vol. 14, pp. 578–578, Aug. 2013, doi: [10.1186/1471-2164-14-578](https://doi.org/10.1186/1471-2164-14-578).

\[Dan19c\] O. Danilevskaya *et al.*, “Developmental and transcriptional responses of maize to drought stress under field conditions,” *Plant Direct*, vol. 3, May 2019, doi: [10.1002/pld3.129](https://doi.org/10.1002/pld3.129).

\[Zha17e\] X. Zhang *et al.*, “Genome-wide identification of gene expression in contrasting maize inbred lines under field drought conditions reveals the significance of transcription factors in drought tolerance,” *PLoS ONE*, vol. 12, Jul. 2017, doi: [10.1371/journal.pone.0179477](https://doi.org/10.1371/journal.pone.0179477).

\[Kas18\] M. Kashima *et al.*, “Genomic Basis of Transcriptome Dynamics in Rice under Field Conditions,” *Plant and Cell Physiology*, vol. 62, pp. 1436–1445, Oct. 2018, doi: [10.1093/pcp/pcab088](https://doi.org/10.1093/pcp/pcab088).

\[Cro17\] R. Cronn *et al.*, “Transcription through the eye of a needle: daily and annual cyclic gene expression variation in Douglas-fir needles,” *BMC Genomics*, vol. 18, Jul. 2017, doi: [10.1186/s12864-017-3916-y](https://doi.org/10.1186/s12864-017-3916-y).

\[Mar26e\] A. Martínez-Pérez *et al.*, “Locally adapted Arabidopsis thaliana accessions show transcriptomic plasticity in a multi‐timescale analysis of whole‐genome gene expression in a natural environment,” *Plant Biology (Stuttgart, Germany)*, vol. 28, pp. 1641–1656, Mar. 2026, doi: [10.1111/plb.70204](https://doi.org/10.1111/plb.70204).

\[Col25b\] B. Cole *et al.*, “Multi‐season analysis reveals hundreds of drought‐responsive genes in sorghum,” *The Plant Journal*, vol. 125, Jul. 2025, doi: [10.1101/2025.06.27.662006](https://doi.org/10.1101/2025.06.27.662006).

\[Har24\] V. Harju *et al.*, “First detection of orchid fleck virus on Veronica spicata and Dendrochilum magnum,” *New Disease Reports*, Oct. 2024, doi: [10.1002/ndr2.12312](https://doi.org/10.1002/ndr2.12312).

\[Wam18\] M. J. Wamaitha *et al.*, “Metagenomic analysis of viruses associated with maize lethal necrosis in Kenya,” *Virology Journal*, vol. 15, May 2018, doi: [10.1186/s12985-018-0999-2](https://doi.org/10.1186/s12985-018-0999-2).

\[Gaa20\] Y. Gaafar *et al.*, “Investigating the Pea Virome in Germany—Old Friends and New Players in the Field(s),” *Frontiers in Microbiology*, vol. 11, Nov. 2020, doi: [10.3389/fmicb.2020.583242](https://doi.org/10.3389/fmicb.2020.583242).

\[Fow21\] A. Fowkes *et al.*, “Integrating High throughput Sequencing into Survey Design Reveals Turnip Yellows Virus and Soybean Dwarf Virus in Pea (Pisum Sativum) in the United Kingdom,” *Viruses*, vol. 13, Dec. 2021, doi: [10.3390/v13122530](https://doi.org/10.3390/v13122530).

\[Mut18\] J. Mutuku *et al.*, “Metagenomic Analysis of Plant Virus Occurrence in Common Bean (Phaseolus vulgaris) in Central Kenya,” *Frontiers in Microbiology*, vol. 9, Dec. 2018, doi: [10.3389/fmicb.2018.02939](https://doi.org/10.3389/fmicb.2018.02939).

\[Xia22e\] H.-G. Xiao, W. Hao, G. Storoschuk, J. MacDonald, and H. Sanfaçon, “Characterizing the Virome of Apple Orchards Affected by Rapid Decline in the Okanagan and Similkameen Valleys of British Columbia (Canada),” *Pathogens*, vol. 11, Oct. 2022, doi: [10.3390/pathogens11111231](https://doi.org/10.3390/pathogens11111231).

\[Hod20\] B. A. Hodge, P. Paul, and L. Stewart, “Occurrence and High-Throughput Sequencing of Viruses in Ohio Wheat.” *Plant disease*, pp. PDIS08191724RE, Apr. 2020, doi: [10.1094/pdis-08-19-1724-re](https://doi.org/10.1094/pdis-08-19-1724-re).

\[Dia23\] N. Dias *et al.*, “Viromes of field-grown tomatoes and peppers in Tennessee revealed by RNA sequencing followed by bioinformatic analysis,” *Plant Health Progress*, Jan. 2023, doi: [10.1094/php-10-22-0107-rs](https://doi.org/10.1094/php-10-22-0107-rs).

\[Kam18\] M. Kamitani, A. Nagano, M. Honjo, and H. Kudoh, “A Survey on Plant Viruses in Natural Brassicaceae Communities Using RNA-Seq,” *Microbial Ecology*, vol. 78, pp. 113–121, Oct. 2018, doi: [10.1007/s00248-018-1271-4](https://doi.org/10.1007/s00248-018-1271-4).

\[Lia21\] Y. Liang, R. Tabien, L. Tarpley, A. R. Mohammed, and E. Septiningsih, “Transcriptome profiling of two rice genotypes under mild field drought stress during grain-filling stage,” *AoB Plants*, vol. 13, Jul. 2021, doi: [10.1093/aobpla/plab043](https://doi.org/10.1093/aobpla/plab043).

\[Mcc26\] B. McCallum *et al.*, “Comprehensive analysis of the wheat leaf virome associated with Puccinia triticina isolates across Canada reveals high viral diversity and novel mycoviruses,” *Phytopathology Research*, vol. 8, Jun. 2026, doi: [10.1186/s42483-026-00443-8](https://doi.org/10.1186/s42483-026-00443-8).

\[Ric12\] C. Richards, U. Rosas, J. A. Banta, N. Bhambhra, and M. Purugganan, “Genome-Wide Patterns of Arabidopsis Gene Expression in Nature,” *PLoS Genetics*, vol. 8, Apr. 2012, doi: [10.1371/journal.pgen.1002662](https://doi.org/10.1371/journal.pgen.1002662).

\[San13c\] S. D. Santo *et al.*, “The plasticity of the grapevine berry transcriptome,” *Genome Biology*, vol. 14, pp. r54–r54, Jun. 2013, doi: [10.1186/gb-2013-14-6-r54](https://doi.org/10.1186/gb-2013-14-6-r54).
