# Field RNA-seq candidate screening

##### [**Undermind**](https://undermind.ai)

---


## Table of Contents

- [Field RNA-seq candidate screening](#field-rna-seq-candidate-screening)
  - [Accession-list screening batch](#accession-list-screening-batch)
  - [Search-result screening batch](#search-result-screening-batch)
  - [Accession-list screening batch two](#accession-list-screening-batch-two)
  - [Search-result screening batch two](#search-result-screening-batch-two)
  - [Search-result screening batch three](#search-result-screening-batch-three)
  - [Search-result screening batch four](#search-result-screening-batch-four)
  - [Search-result screening batch five](#search-result-screening-batch-five)
  - [Accession-list screening batch three](#accession-list-screening-batch-three)
  - [Rules applied in this audit](#rules-applied-in-this-audit)
  - [References](#references)

# Field RNA-seq candidate screening

**Purpose.** Classify candidate papers against the field-collected bulk RNA-seq criteria: primary research; total RNA or mRNA sequencing; aerial plant tissue; tissue collected outdoors in a field, farm, orchard, vineyard, or forest (including field trials; excluding glasshouse/chamber and detached-leaf assays); and public raw reads in SRA, ENA, DDBJ, or GSA. A paper is marked **meets** only when all required elements are supported. **Unresolved** means at least one required element still needs verification. Reanalyses of an already-counted dataset are linked rather than counted as new collections.

**Coverage.** As of October 5, 2026, this partial audit has screened 103 of the 208 papers in the search results and 54 distinct accession rows from the inventory supplied in chat. It is not a complete census; unresolved candidates and the remaining records still need review. The existing eligible-studies table remains provisional until both candidate pools have been reconciled.

| Paper | Decision | Decisive evidence and reason | Accession or follow-up |
|:---|:---|:---|:---|
| \[Wam18\] | Unresolved | The paper reports bulk, non-polyA-selected total RNA from leaves of maize, sorghum, and Napier grass sampled on Kenyan farms. However, its stated BioProject does not resolve to plant reads. | Paper reports PRJNA42371, but NCBI identifies that project as *Amphibacillus xylanus* genome sequencing. Need to locate a corrected raw-read accession before counting it. citeturn0search0 |
| \[Elm22\] | Meets | Bulk ribodepleted total-RNA sequencing of symptomatic field-collected soybean leaves, stems, and pods; samples were collected from soybean fields across eight U.S. states. | Raw reads reported under PRJNA730888, PRJNA730918, and PRJNA731151. |
| \[Gaa20\] | Does not meet | Field pea samples and ribodepleted total-RNA sequencing are supported, but the paper deposits assembled viral sequences rather than public raw reads in the specified archives. | No qualifying raw-read accession located in the paper’s data statement. |
| \[San13c\] | Does not meet | Grapevine berries were field-collected, but expression was measured by NimbleGen microarray, not RNA sequencing. | GEO microarray series GSE41633 is not a raw-read archive. |
| \[Nag12\] | Does not meet | Rice leaves were sampled in a paddy field, but the study’s expression measurements were DNA microarrays, not bulk RNA-seq. | GEO microarray series; no qualifying raw reads for this study. |
| \[Nag19\] | Meets | Seasonal and diurnal bulk transcriptomes were collected from leaves of wild *Arabidopsis halleri* in a natural habitat over two years; the study also includes separate controlled-environment experiments, which should not be counted as field samples. | DDBJ SRA accessions DRA005871–DRA005876. Confirm field-only sample count during final record reconciliation. |
| \[Sch25d\] | Meets | Bulk ribodepleted total-RNA sequencing of leaves and shoot tips from blackberry and wild *Rubus* collected at South Carolina farms and nearby sites. | SRA run accessions reported in the paper. |
| \[Han16\] | Unresolved | The paper describes field-collected soybean leaves and bulk total-RNA sequencing, pooled first by region and then into two yearly sequencing pools. | A qualifying public raw-read accession has not yet been verified; do not count pending archive check. |
| \[Kam18\] | Unresolved | The abstract supports RNA-seq with rRNA depletion in wild Brassicaceae communities, but tissue-level sampling details and the raw-read accession have not yet been checked against full text/archive records. | Also check whether its reads overlap a previously counted *A. halleri* dataset. |
| \[Red21b\] | Meets | Bulk ribodepleted total-RNA sequencing of field-collected wheat leaves; public raw reads are reported. | PRJNA722004. Check for related papers using the same project so the dataset is not double-counted. |
| \[Mut18\] | Does not meet | Bulk total-RNA sequencing of common-bean leaves collected from Kenyan farms is supported, but the paper reports GenBank accessions for assembled viral sequences, not public raw reads in a specified archive. | No qualifying raw-read accession located in the paper’s data statement. |
| \[Mat17b\] | Does not meet | The paper includes bulk transcriptome sequencing of citrus leaves (and roots) from field groves, but its reported accessions are assembled viral sequences, not deposited raw reads. | No qualifying raw-read accession located in the paper. |
| \[Hon19b\] | Same dataset, not a new collection | The paper analyzes field-collected *A. halleri* leaf transcriptomes and discusses uncontrolled pathogen status, but its RNA-seq data are from an earlier seasonal transcriptome dataset. Link to the data-generating paper rather than count another collection. | DRA005871–DRA005876 overlap the \[Nag19\] seasonal transcriptome dataset; verify exact sample overlap when consolidating. |
| \[Mul22\] | Unresolved | Full-text methods and archive records have not yet been verified. | Pending screening. |
| \[Obo24\] | Unresolved | Full-text methods and archive records have not yet been verified. | Pending screening. |

## Accession-list screening batch

This batch checks the first set of higher-count accessions from the inventory supplied in chat. A project can qualify even when its row has the wrong title or DOI, but those metadata mismatches are recorded and must be corrected before the final table is merged. Projects are not counted as separate studies when they reuse the same field collection.

| Supplied accession and candidate record | Decision | Decisive evidence and reason | Crosswalk note |
|:---|:---|:---|:---|
| PRJEB39201 — 539 wheat field samples | Meets | \[Ada21\] generated 538 new bulk total-RNA datasets from field-collected infected wheat leaves across 30 countries. | Accession is the new field dataset in \[Ada21\]. The DOI/title in the supplied row points to a different paper; correct the bibliographic link. |
| PRJNA1217477 — Botrytis multiplant atlas | Does not meet | \[Sin26b\] generated bulk host–pathogen RNA-seq from inoculated leaves in controlled experiments, not field-collected tissue. | Accession and paper match; controlled inoculation is disqualifying. |
| PRJEB31334 — listed as a Rust expression browser record | Meets, metadata mismatch | \[Rad19\] reports bulk total-RNA sequencing of naturally infected wheat leaves collected from fields in nine countries; reads are deposited under ERP113880/PRJEB31334. | The supplied title/DOI fields are mismatched: DOI 10.1186/s12915-019-0684-y is the MARPLE paper, not the Rust expression browser. |
| PRJEB15280 — 117 wheat samples | Meets, reused dataset | \[Bue17\] reports field-collected wheat, triticale, and rye leaf RNA-seq, with raw reads under PRJEB15280. | This is an existing field-pathogenomics dataset reused in later resources; do not count it as new data generated by the later database paper. |
| PRJEB84390 — wheat rust study | Meets | \[Tsu25\] includes bulk RNA-seq of field-collected *Puccinia triticina*-infected wheat leaves, alongside separate pathogen-isolate genomic work; the RNA data are deposited in ENA. | Count only the field plant RNA-seq component, not the pathogen-isolate DNA sequencing. |
| PRJNA513863 — grape deacclimation/budbreak | Unresolved | The field-versus-controlled provenance and exact sequenced tissue have not yet been verified from full text. | The DOI/title may refer to a related preprint/published version; confirm accession-to-paper linkage. |
| PRJNA526829 — *Magnaporthe* MoAa91 | Does not meet | \[Li20e\] sequenced controlled rice seedling infection samples collected at defined post-inoculation times, not field tissue. | Public expression data do not overcome the controlled-growth exclusion. |
| PRJEB65589 — wheat stem rust in western Europe | Unresolved | Full methods and the contents of the listed project have not yet been reconciled. | Do not infer host-plant RNA-seq from a field disease report or an accession count alone. |
| PRJNA256347 — field pathogenomics of stripe rust | Meets | \[Hub15\] generated bulk RNA-seq from field-collected infected wheat and triticale leaves in the UK and deposited raw reads. | One of the already established field pathogenomics collections. |
| PRJNA746402 — maize senescence-associated traits | Meets | \[Cai21b\] reports ear-leaf RNA-seq from field-grown maize at Pontevedra, Spain, with raw reads under this BioProject. | The paper cites the underlying thesis dataset; treat this as one field collection. |
| PRJEB75530 — *Fusarium*–wheat interaction | Unresolved | The DOI supplied with this row maps to a different paper \[Dil19\], whose RNA-seq used controlled inoculations. The listed accession/title pairing has not been verified. | Keep separate from \[Dil19\] until the project-to-paper link is established. |
| PRJEB36485 — South African stripe-rust population | Unresolved | Need to establish whether the project contains bulk infected wheat tissue or pathogen-only isolate sequencing. | Field origin of a rust isolate does not by itself qualify as field-collected plant RNA. |
| PRJNA328045 — *Leptosphaeria maculans*–canola interaction | Does not meet | \[Spa16\] used controlled seedling inoculations and sampled leaves/cotyledons at defined post-inoculation stages. | Not field-collected plant tissue. |
| PRJNA950118 — Canadian stripe-rust population structure | Meets, subset only | \[Hol23c\] reports 18 libraries made from RNA extracted directly from single-lesion infected wheat leaves; raw reads include host and pathogen RNA. Other project libraries are purified pathogen urediniospores and do not count as plant tissue. | Count only the 18 field leaf libraries; do not count pathogen-only libraries or reused sequences. |
| PRJEB46918 — Asian soybean rust genome/transcriptome | Does not meet | \[Gup23\] used pathogen isolate/genome material and controlled inoculated soybean leaves for transcriptomics; it is not field-collected host RNA-seq. | Field origin of the original isolate is not the collection setting for the transcriptome libraries. |
| PRJEB47693 — wheat stem rust in Ireland | Meets | \[Tsu22\] generated bulk RNA-seq from field-collected infected wheat leaves and stems, including Irish field plots, and deposited reads in ENA. | The project also includes other European samples; use the paper’s sample-level breakdown. |
| PRJEB59062 — wheat susceptibility to *Fusarium* | Does not meet | \[Roc24\] used controlled growth-chamber inoculations of wheat spikelets, not field-collected tissue. | Raw reads are public, but field provenance is absent. |
| PRJEB48949 — Furmint berries and noble rot | Meets | \[Pog22\] used poly(A) RNA-seq on grape berries collected during natural noble-rot development in a Hungarian vineyard and deposited raw reads in ENA. | Eligible aerial fruit tissue; keep library/sample count separate from the number of SRA runs. |
| PRJNA1013336 — lettuce during Botrytis infection | Unresolved | The paper’s collection setting and whether RNA-seq came from field samples or controlled infection have not yet been verified. | Do not assume field provenance from the host or accession label. |
| PRJNA1279063 — finger-millet-infecting *Magnaporthe* | Does not meet | The linked study profiles fungal strains/pathogen genomes, not bulk RNA from field-collected aerial plant tissue. | Pathogen collected from a crop is not the same as sequencing plant tissue. |
| PRJNA802726 — UK pea virus survey | Unresolved | \[Fow21\] confirms bulk ribodepleted RNA-seq on field pea-leaf pools, but its paper data statement lists GenBank viral-sequence records rather than raw-read accessions. | The supplied BioProject needs direct verification that it contains the study’s plant RNA-seq reads before this candidate can count. |
| PRJNA540228 — barley transcriptome-wide association study | Unresolved | Tissue provenance, assay design, and raw-read linkage have not yet been checked against the paper and archive. | Pending screening. |

## Search-result screening batch

| Paper | Decision | Decisive evidence and reason | Accession or follow-up |
|:---|:---|:---|:---|
| \[Cho26c\] | Meets | Pepper-leaf libraries used bulk total-RNA depletion and poly(A) selection on commercial field-collected plants. Garlic bulbs are excluded, but the pepper leaves independently satisfy the tissue criterion. | PRJNA1420407. |
| \[Cho20d\] | Meets | Bulk total-RNA sequencing of tomato and pepper leaves collected from farms/open fields in Vietnam. | PRJNA636575. |
| \[Dan24\] | Unresolved | Full text was not available for checking whether the HTS libraries were qualifying bulk RNA-seq and whether raw reads were deposited. | Verify assay and archive before counting. |
| \[Rum19\] | Does not meet | Bulk leaf RNA-seq is described, but the paper reports assembled viral sequences and no raw-read archive for its own libraries. | A cited SRX accession belongs to another study. |
| \[Ram22b\] | Meets | Poly(A)-selected bulk RNA-seq of naturally infected citrus leaves from commercial orchards in Veracruz, Mexico. | PRJNA748945. |
| \[Dia23\] | Unresolved | Full text and raw-read archive not verified. | Check field collection, library type, and project contents. |
| \[Xia22e\] | Does not meet | Field-orchard apple leaves were sequenced after rRNA depletion, but the paper reports viral contig accessions rather than raw reads in a specified archive. | No qualifying raw-read accession found in the paper. |
| \[Kim21g\] | Meets | Bulk rRNA-depleted total-RNA sequencing of pear leaves collected from 35 Korean orchards. | SRA run SRR8466619. |
| \[Hu21c\] | Unresolved | PDF and accession linkage not verified. | Check that the study generated bulk plant RNA-seq from outdoor grape tissue. |
| \[Hod20\] | Unresolved | PDF and raw-read accession not verified; the title alone does not establish that the assay was bulk plant RNA-seq rather than another HTS workflow. | Check sample tissue, library method, and SRA records. |
| \[Jo18\] | Meets | Poly(A)-selected mRNA-seq of field-collected barley leaves sampled across 17 Korean regions. | SRA runs SRR6706097–SRR6706102. |
| \[Nhl18\] | Unresolved | Full text not available to determine whether sequencing was bulk total RNA/mRNA or small-RNA/targeted virome sequencing. | Check assay and raw-read archive. |
| \[Esc26\] | Meets | Bulk rRNA-depleted RNA-seq of cotton leaves and petioles collected from field sites across six southern U.S. states. | PRJNA1015735. |
| \[Wri20\] | Meets | Ribodepleted bulk total-RNA sequencing included apple leaves from commercial orchards in Washington; root libraries are not counted. | PRJNA562540. |
| \[Jo18b\] | Meets | Poly(A)-selected mRNA-seq of peach leaves collected from orchard-grown trees in Korea. | SRA run accessions listed in the paper. |
| \[Von21\] | Meets | Bulk berry-pericarp RNA-seq from grapevines planted in an outdoor experimental vineyard; the plants were field-sampled after greenhouse acclimation. | RNA BioProject PRJNA701940. |
| \[Pag20\] | Meets | Bulk RNA-seq of grapevine leaf veins from naturally infected and recovering field vines in an Italian vineyard. | PRJNA587882. |
| \[Bel23c\] | Meets | Ribo-depleted bulk total-RNA sequencing of grapevine leaf petioles and veins from an outdoor ampelographic collection in Russia. | PRJNA1043183. |
| \[Jo20\] | Meets | Poly(A)-selected bulk mRNA-seq of leaves from field-grown sweet-potato plants sampled across Korean regions. | PRJNA517178. |
| \[Lov18\] | Meets | Poly(A)-selected bulk mRNA-seq of *Panicum hallii* leaves harvested outdoors at a field site during a drought experiment. The user’s criteria include field-grown studies; a field experiment is not excluded merely because it has a treatment. | PRJNA251785 and PRJNA250527. |

## Accession-list screening batch two

| Supplied accession and candidate record | Decision | Decisive evidence and reason | Crosswalk note |
|:---|:---|:---|:---|
| PRJNA151285 — *Colletotrichum graminicola* on maize | Does not meet | \[Oco12\] sequenced excised maize leaf sheaths after controlled fungal inoculation, not field-collected tissue. | Host tissue is present, but the controlled assay setting is disqualifying. |
| PRJNA1231252 — PstS10 virulence-gain study | Does not meet as an independent field collection | \[Car26\] analyzes previously published field datasets but the reads in this project are from newly generated greenhouse-inoculated wheat leaves. | Link the field-data reanalysis to \[Ada21\]; do not count this project as a new field collection. |
| PRJNA853492 — tobacco–*Rhizoctonia* time course | Does not meet | \[Li22g\] used controlled inoculations of tobacco leaves harvested at set post-inoculation times. | Field origin of the fungal isolate does not make the sequenced leaves field-collected. |
| PRJNA505891 — maize senescence | Meets | \[Sek19\] sequenced bulk leaf and internode tissue collected from field trials at the Clemson University research farm. | Field trial; count aerial leaf/internode libraries. |
| PRJEB8798 — wheat–*Zymoseptoria* infection cycle | Does not meet | \[Rud15\] used greenhouse-inoculated wheat and pathogen culture samples, not field-collected plant tissue. | Controlled infection time course. |
| PRJNA1113144 — banana TLP16/Fusarium wilt | Does not meet | The transcriptome data are from controlled inoculations and root tissue, outside the allowed aerial-tissue set. | The accession’s RNA-seq is not eligible even though the host is a crop. |
| PRJNA1055939 — Arabidopsis field RNA-seq | Meets, same study as other projects | \[Tom26\] sampled leaf tissue from outdoor sites in Japan and Switzerland and deposited the field reads across several BioProjects, including this one. | Link all eight projects to the same field study; do not count each project as a separate study. |
| PRJEB34186 — barley STARR-seq | Does not meet | STARR-seq is a reporter/enhancer assay using plasmid-derived transcripts, not bulk endogenous plant RNA-seq from field tissue. | Assay criterion fails. |
| PRJEB22223 — *Puccinia graminis* isolate sequencing | Does not meet | The sequenced material is pathogen isolates/genomes, not bulk host plant aerial tissue. | Pathogen field provenance alone is insufficient. |
| PRJNA1183537 — wheat Fusarium-head-blight landraces | Meets | \[Wu25c\] sampled bulk wheat spikelets from plants growing in an outdoor breeding field and deposited raw reads. Plants were point-inoculated in the field; this is not a greenhouse or detached-leaf assay. | Mark as field-collected, artificially inoculated—not natural infection surveillance. |
| PRJNA540228 — barley rpg4/stem-rust transcriptomics | Does not meet | \[Pou19\] used greenhouse-inoculated barley seedlings and sampled leaves at a fixed post-inoculation time. | RNA-seq data are plant–pathogen mixtures, but collection was controlled. |
| PRJNA281236 — Sémillon noble rot | Meets | Bulk RNA-seq of ripe Sémillon berries naturally infected in a commercial vineyard over three years; raw reads are linked to this BioProject. | PRJNA281236; 16 SRA experiments. |

## Search-result screening batch two

| Paper | Decision | Decisive evidence and reason | Accession or follow-up |
|:---|:---|:---|:---|
| \[Mje26\] | Unresolved | Full text and raw-read linkage not verified. | Check exact outdoor samples, aerial tissue, and archive. |
| \[Dia22b\] | Unresolved | Full text and raw-read linkage not verified. | Check whether it used bulk RNA-seq rather than a targeted virus assay. |
| \[Kha25b\] | Does not meet | The study used bulk ribodepleted total RNA from field soybean leaves, but the paper lists assembled viral-sequence accessions and no qualifying raw-read archive. | No SRA/ENA/DDBJ/GSA accession established in the paper. |
| \[Ber26\] | Meets | Bulk grapevine leaf RNA sequencing from symptomatic, asymptomatic, and control vines sampled in an outdoor vineyard. | ENA PRJEB94866. |
| \[Nab22\] | Does not meet | Bulk ribodepleted RNA-seq of field apple leaves is described, but raw reads are not deposited in a qualifying archive; the data statement offers data on request and lists viral assemblies. | No qualifying raw-read accession found. |
| \[Rod26\] | Unresolved | Full text and raw-read archive not verified. | Check orchard tissue, library type, and raw-read project. |
| \[Xu19\] | Does not meet | The study used small-RNA sequencing from greenhouse nectarine trees and reports assembled viral sequences rather than public raw reads. | Fails assay, setting, and archive criteria. |
| \[Kim24c\] | Meets | Bulk rRNA-depleted total-RNA sequencing of symptomatic citrus leaves sampled from Korean production regions. | PRJNA928104; runs SRR23239772–SRR23239774. |
| \[Sur26\] | Meets | Bulk rRNA-depleted RNA-seq of symptomatic pepper leaves collected from field sites in West Java, Indonesia. | PRJNA1492756. |
| \[Kus20\] | Unresolved | Full text unavailable; setting, library, and raw-read archive not yet checked. | Pending screening. |
| \[Lia21\] | Does not meet | Field-grown rice leaves were sequenced by bulk RNA-seq, but the paper reports supporting information rather than a qualifying public raw-read accession. | No SRA/ENA/DDBJ/GSA reads established. |
| \[Wu23d\] | Does not meet | The study sequenced dsRNA-enriched virome material from field grapevine tissue, not bulk total RNA or mRNA as required. | Viral assemblies are reported; no qualifying bulk-RNA archive accession established. |
| \[Jo20c\] | Meets | Poly(A)-selected bulk RNA-seq of leaf samples from orchard-grown plum cultivars in Korea. | PRJNA295439. Some viral reads were attributed to index misassignment; this does not alter sequencing-method eligibility. |
| \[San16c\] | Does not meet | Vineyard-collected grape berries were analyzed by NimbleGen microarray, not RNA sequencing. | GEO GSE75565 is microarray expression data. |
| \[San25b\] | Does not meet | Orchard/nursery bark was analyzed using HTS, but only assembled virus accessions are reported; no qualifying public raw reads were found. | GenBank accessions are viral assemblies, not raw reads. |
| \[Wil16d\] | Same field dataset, not a new collection | This paper analyzes field rice-leaf RNA-seq among broader controlled and field datasets. The field data overlap the \[Ple15\] collection already listed in the main table. | Link to \[Ple15\]; do not count a second field collection. |
| \[Lee23h\] | Meets | Bulk ribodepleted total-RNA sequencing of symptomatic wheat leaves sampled in South Korean wheat-growing regions. | SRA runs SRR16774428–SRR16774430 and SRR22371364–SRR22371366. |
| \[Rat22b\] | Meets | Bulk RNA-seq of leaves collected over 12 seasonal time points from a natural alpine habitat in India. | PRJNA751462. |
| \[Jo21\] | Meets, field subset only | Poly(A)-enriched mRNA-seq included overwintering pepper fruit from two open-field regions; a third sample came from a greenhouse and is excluded. | SRA runs SRR13319506–SRR13319511; identify the field runs in final sample-level crosswalk. |
| \[Maa21\] | Meets, library subset | Field-collected tomato leaves were sequenced in two ribodepleted total-RNA libraries and one dsRNA-enriched library. The bulk total-RNA libraries satisfy the assay criterion; the dsRNA-enriched library alone does not. | Runs SRR14066946–SRR14066948; distinguish the total-RNA libraries from the dsRNA library. |

## Search-result screening batch three

| Paper | Decision | Decisive evidence and reason | Accession or follow-up |
|:---|:---|:---|:---|
| \[Zha22k\] | Does not meet | The RNA-seq leaves and stems were sampled from plants grown and drought-treated in a greenhouse; field work was limited to phenotyping. | PRJEB43702. |
| \[Xia19\] | Does not meet | Bulk rRNA-depleted RNA-seq of phloem scrapings from vineyard grape cuttings is supported, but only assembled viral genomes are deposited; no raw-read archive is reported. | No qualifying raw-read accession located. |
| \[Cho16b\] | Does not meet | Bulk rRNA-depleted total-RNA sequencing of orchard apple leaves is supported, but the paper provides no raw-read accession. | Only viral reference accessions are listed. |
| \[Goz25\] | Meets | rRNA-depleted bulk RNA sequencing of symptomatic leaf pools from peach/nectarine trees sampled in Bulgarian orchards in 2022–2023. The paper describes the libraries as lncRNA sequencing, but the workflow starts from total RNA and uses rRNA depletion rather than targeted capture. | SRA PRJNA1219823. |
| \[Zha20r\] | Does not meet | Bulk rRNA-depleted RNA-seq of field-collected Camellia leaves is supported, but only assembled virus sequences are reported, not raw reads. | No qualifying raw-read accession found. |
| \[Ye23\] | Meets | Bulk RNA-seq of leaf and flower buds collected outdoors at four elevations in the Gaoligong Mountains. | PRJNA916369. |
| \[Pat25c\] | Meets | Field-collected grapevine spur/cortical tissue was sequenced using stranded mRNA libraries; the GEO series links to raw SRA data and BioProject PRJNA1192272. | GEO GSE283169; SRA BioProject PRJNA1192272. |
| \[Die23\] | Meets, plant samples only | Bulk rRNA-depleted total-RNA sequencing includes field-collected tree bark, wood, and leaves from Mediterranean forests. | PRJNA1032577; exclude arthropod and cultured-fungus samples from the plant-tissue count. |
| \[Amo26\] | Meets | Bulk rRNA-depleted total-RNA sequencing of field-collected leaves from tomato, pepper, bean, and a weed in Angola. | Raw SRA runs reported under PRJNA1399677 and PRJNA1470426. |
| \[Mur25\] | Same dataset, not a new collection | This paper reanalyzes the field leaf transcriptomes from \[Nag19\]. | DRA005871–DRA005876; link to \[Nag19\]. |
| \[Pao21\] | Meets | Bulk metatranscriptomic RNA from wood of mature field-grown grapevines; this is the data-generating study underlying the later viral reanalysis. | ENA PRJEB31098; link the later \[Deb25b\] analysis to this collection. |
| \[Ras22b\] | Does not meet | Bulk rRNA-depleted total RNA from field bean leaves was sequenced, but the paper reports assembled viral genomes and no qualifying raw-read deposit. | No SRA/ENA/DDBJ/GSA raw-read accession established. |
| \[Har21e\] | Meets | Bulk 3′-RNA-seq of grapevine leaves sampled in an outdoor experimental vineyard. | PRJNA674915; cross-check for overlap with the later \[Har23c\] paper before counting studies. |
| \[Nan25c\] | Does not meet | Bulk total-RNA sequencing of field cereal tillers is described, but the listed GenBank accessions are assembled virus genomes, not deposited raw reads. | No qualifying raw-read accession found. |
| \[Ren22g\] | Meets | Bulk RNA-seq of developing rice grains collected from field-grown panicles in a paddy. | SRA runs SRR17706479, SRR17709283, SRR17709378, SRR17709379, and SRR17709726. |
| \[Lov16\] | Meets, field samples only | Bulk tag-based mRNA-seq includes aerial leaves sampled from outdoor field plots/shelters; the study also includes greenhouse samples. | PRJNA322529; count only the outdoor field libraries. |
| \[Bon19\] | Unresolved | Full text and raw-read archive not available for verification. | Pending screening. |
| \[Rog23c\] | Meets | Poly(A)-selected bulk RNA-seq of stem xylem/cambium from trees grown outdoors in a common-garden field experiment. | SRP188754. |
| \[Man25e\] | Unresolved | Full text and archive linkage not verified. | Pending screening. |
| \[Has21c\] | Meets, field samples only | The primary study generated bulk RNA-seq from field-grown rice leaves as well as controlled-environment samples. | PRJNA726716; count only field samples. |

## Search-result screening batch four

| Paper | Decision | Decisive evidence and reason | Accession or follow-up |
|:---|:---|:---|:---|
| \[Gla20\] | Does not meet | Bulk rRNA-depleted RNA-seq of a field-grown grapevine leaf is supported, but raw reads are not deposited; only assembled viral/viroid sequences are reported. | No qualifying raw-read accession found. |
| \[Zag26\] | Meets | Bulk total-RNA sequencing of field-collected pepper leaves in Iraq. | SRA run SRR21032823. |
| \[Min25\] | Meets, field subset only | Bulk total-RNA sequencing includes field-collected lily leaves and flowers; other libraries include greenhouse-grown or bulb samples and are excluded. | PRJNA1313793; count only field aerial-tissue libraries. |
| \[Fas12b\] | Does not meet as an independent RNA-seq study | The grapevine atlas uses newly generated microarrays and reanalyzes older RNA-seq datasets; it does not generate new bulk RNA-seq reads. | Link any qualifying source datasets to their original studies. |
| \[Bel24\] | Meets | Bulk rRNA-depleted total-RNA sequencing of Cnidium leaves collected from South Korean farms. | PRJNA1074493. |
| \[Mar16c\] | Unresolved | Full text and a valid raw-read accession have not been verified. | Pending screening. |
| \[Van22c\] | Does not meet | Although the paper includes field switchgrass leaf mRNA-seq, the public raw-read projects identified in the paper correspond to microbial amplicon or host genetic data; host RNA-seq is available as expression counts, not qualifying raw reads. | Do not mistake ITS/LSU amplicon reads for plant RNA-seq reads. |
| \[Cho23c\] | Same dataset, not a new collection | The study reanalyzes public poly(A)-selected soybean transcriptomes from field-grown buds, flowers, and pods. | SRA PRJNA481793; link to the data-generating study and do not count a new collection here. |
| \[Pad23\] | Meets | Bulk rRNA-depleted total-RNA sequencing of symptomatic citrus leaves collected in field surveys across Colombia and the United States. | PRJNA895362. |
| \[Yu26b\] | Meets, sorghum field subset | The sorghum component used bulk RNA-seq of leaves sampled in an outdoor Arizona field trial; other plant components of the paper are not field RNA-seq. | PRJNA1119650. |
| \[Dik25c\] | Unresolved | Full text and raw-read linkage have not yet been verified. | Pending screening. |
| \[Nak22\] | Unresolved | Full text and exact tissue/library/archive details were not available for verification. | Confirm that eligible aerial tissue—not storage roots alone—was sequenced and that raw reads are public. |
| \[Nis24\] | Same dataset, not a new collection | Field measurements in this paper are RT-qPCR; its RNA-seq component reuses the \[Nag19\] dataset. | DRA005871–DRA005876; link to \[Nag19\]. |
| \[Deb25b\] | Same dataset, not a new collection | This paper reanalyzes the field grapevine wood metatranscriptomes generated by \[Pao21\]. | ENA PRJEB31098; do not count a second sampling study. |
| \[Col26\] | Reanalysis, not a new collection | The paper reanalyzes field leaf mRNA-seq from already identified maize, sorghum, and soybean studies to profile fungal signal. | PRJEB83049, PRJEB67964, and CRA009979; link to the data-generating studies. |
| \[Wei25c\] | Does not meet | Bulk leaf RNA-seq and public reads are reported, but plants were grown and sampled in a controlled glasshouse. | PRJNA1094407. |
| \[Paa26\] | Meets | Bulk BrAD-seq of leaves from a natural *Arabidopsis halleri* population sampled outdoors across seasons. | ENA PRJEB83648. |
| \[Mwa20\] | Unresolved | The paper reports field surveillance and NGS, but the exact plant library type and public raw-read deposit have not been confirmed. | Check whether the listed project contains raw plant RNA-seq rather than only viral assemblies. |
| \[Kha23c\] | Unresolved | Field pear leaves and total-RNA sequencing are supported by the paper record, but a qualifying raw-read accession has not yet been verified. | Search the linked project/SRA records before counting. |
| \[Dan21\] | Does not meet | Field-grown sugarcane tissue was analyzed by oligoarray microarrays, not RNA-seq. | GEO GSE129543 and GSE171222 are expression-array records. |

## Search-result screening batch five

| Paper | Decision | Decisive evidence and reason | Accession or follow-up |
|:---|:---|:---|:---|
| \[Alb20\] | Unresolved | The paper surveys field wheat and mite-transmitted viruses, but the specific sequenced material, assay, and raw-read archive have not been verified. | Check whether any libraries are bulk plant RNA-seq, rather than mite/virus sequence analyses. |
| \[Cho18b\] | Search-record mismatch | The cite-key record resolves to an unrelated hepatitis C paper, not the apparent apple-virus study. Do not treat this search hit as an eligible plant study until the correct bibliographic record is identified. | Repair the record mapping; the apple-virus paper is a separate candidate. |
| \[Mat16\] | Does not meet | The thesis describes citrus virus sequencing associated with the same dataset as \[Mat17b\], but no qualifying public raw-read accession was located; only viral sequence records are reported. | Link to \[Mat17b\] and avoid counting the same dataset twice. |
| \[Zyk25\] | Unresolved, likely reanalysis | The article analyzes five in-house rye datasets plus 50 public SRA datasets to screen non-rye reads. Whether its in-house accessions are new field-collected samples and whether their reads are publicly deposited has not been established. | Do not count the 50 reused datasets as new collections; verify the five in-house records. |
| \[Yu25e\] | Same study/version as \[Yu26b\] | This preprint and the later article describe the same outdoor sorghum field RNA-seq study. | One study, one accession: PRJNA1119650. |
| \[Hod17\] | Does not meet | Bulk rRNA-depleted RNA-seq was performed on pooled field-collected wheat leaves, but the short report deposits the assembled viral sequence, not the raw reads in a qualifying archive. | GenBank viral sequence MF621330 is not a raw-read accession. |
| \[Kaw20\] | Unresolved | The RNA-seq tissue, growth setting, and public raw-read accession could not be verified from the available full text. | Pending verification of field versus controlled-grown samples. |
| \[Rie19\] | Unresolved | The transcriptome methods, field provenance, and raw-read accession have not been verified. | Pending full-text/archive check. |

## Accession-list screening batch three

| Supplied accession and candidate record | Decision | Decisive evidence and reason | Crosswalk note |
|:---|:---|:---|:---|
| PRJEB44222 — *Phakopsora pachyrhizi* isolate PpUFV02 | Does not meet | \[Gup23\] includes pathogen-isolate sequencing and controlled soybean inoculations, not field-collected host RNA-seq. | Pathogen isolate reads are not eligible plant-tissue reads. |
| PRJNA577523 — BIOCOMES RNASeq | Unresolved | The accession is associated with RNA-seq on the wheat powdery-mildew antagonist *Golubevia*, but the supplied DOI points to an unrelated fungal genome note; the project’s sample material and field provenance need checking. | Do not count until the BioSamples establish field-collected aerial plant RNA. |
| PRJEB59238 — wheat Fusarium-head-blight RNA-seq | Does not meet | \[Roc24\] used controlled growth-chamber inoculations of wheat spikelets. | Public raw reads do not satisfy the field-setting criterion. |
| PRJNA1245489 — Xinluzao84/J8031 cotton RNA-seq | Does not meet | \[Li25f\] sequenced cotton roots after controlled *Verticillium* inoculation, not field-collected aerial tissue. | Accession/title in the supplied row is misleading about tissue and setting. |
| PRJNA407369 — maize–*Ustilago maydis* time course | Does not meet | \[Lan18\] used inoculated maize seedlings grown in growth chambers. | Controlled infection assay. |
| PRJNA415866 — Australian *Puccinia graminis* pathotypes | Does not meet | The project concerns pathogen pathotype/genome sequencing, not bulk host plant RNA from field tissue. | Pathogen-only dataset. |
| PRJNA1055939 — Arabidopsis field transcriptome | Meets, duplicate project within one study | \[Tom26\] includes field leaf libraries from outdoor sites; this is one of several BioProjects for the same field study. | Link to \[Tom26\]; do not count each project separately. |
| PRJEB34186 — barley STARR-seq | Does not meet | The assay is a reporter-enhancer screen, not endogenous bulk plant RNA-seq. | Assay criterion fails. |
| PRJEB22223 — *Puccinia graminis* isolate sequencing | Does not meet | The libraries are pathogen isolate/genome data, not plant-tissue RNA-seq. | Pathogen-only. |
| PRJNA994854 — “Tomato raw sequence reads” | Unresolved, metadata mismatch | The supplied DOI resolves to an unrelated yak gut virome paper, not a tomato transcriptome study. | Confirm the correct paper and inspect BioSamples before deciding. |
| PRJNA1183537 — wheat FHB landrace study | Meets | \[Wu25c\] used bulk RNA-seq from field-grown wheat spikelets; plants were point-inoculated in the outdoor breeding field. | Field-inoculated sample, not natural infection surveillance. |
| PRJNA281236 — Sémillon noble rot | Meets | Bulk RNA-seq of naturally infected field berries from a commercial Napa vineyard, with raw reads in SRA. | PRJNA281236 has 16 SRA experiments; the supplied table lists 12, so reconcile the count at sample level. |
| PRJNA274853 — canola–*Sclerotinia* transcriptome | Unresolved, metadata mismatch | The supplied DOI resolves to an in-silico Brassica gene-family paper, not the listed infection transcriptome. | Verify which paper and biological samples the BioProject actually contains. |
| PRJNA1428298 — Botrytis pilot co-transcriptome | Does not meet | \[Sin26b\] used controlled inoculated leaves, not field-collected plants. | Controlled pilot experiment. |
| PRJEB33109 — rust expression browser record | Unresolved, metadata mismatch | The supplied DOI maps to a different wheat-rust population paper; the project’s source-study and field sample links are unclear. | Check BioProject/BioSample records before merging with the \[Ada21\] datasets. |
| PRJNA1405621 — *Xanthomonas perforans* strain reads | Does not meet | The record describes bacterial strain sequencing, not bulk plant-tissue RNA. | Pathogen-only. |
| PRJNA925193 — cacao *Moniliophthora roreri* population genomics | Does not meet | The project sequences fungal isolates/genomes, not bulk cacao tissue. | Pathogen-only. |
| PRJNA674985 — ThatcherLr14a/Puccinia infection | Does not meet | \[Kol21\] generated RNA-seq from controlled inoculations of wheat leaves; its field trials were phenotyping only. | Controlled inoculation libraries. |
| PRJNA1414147 — tomato–*Phytophthora infestans* transcriptome | Unresolved | The preprint/project linkage and sample collection setting have not been verified from full text. | Check whether the RNA-seq is field-collected or controlled infection. |
| PRJNA895362 — citrus leprosis virus surveillance | Meets | \[Pad23\] used bulk rRNA-depleted total RNA from field-collected symptomatic leaves. | Public raw reads under PRJNA895362. |
| PRJNA293991 — *Colletotrichum*–maize leaf sheath interaction | Does not meet | The dataset is from controlled infection/developmental assays, not field-collected tissue. | No eligible field sample series established. |
| PRJNA554133 — grass–*Epichloë* transcriptomes | Meets | \[Ber22\] sequenced bulk mRNA from outdoor-grown grass reproductive tissues under cover in Lexington, Kentucky; reads include host and fungal transcripts. | PRJNA554133; count the plant-containing tissue libraries, not a fungus-only interpretation. |
| PRJNA671670 — *Ralstonia solanacearum* stress response | Does not meet | The RNA-seq is from bacteria exposed to nitrosative/oxidative stress, not plant tissue. | Pathogen-only in vitro assay. |
| PRJNA1166635 — *Erwinia amylovora* apple-flower colonization | Does not meet | Although flower samples were inoculated in a research orchard, the libraries were depleted of host poly(A) RNA and prepared to profile bacterial mRNA, not bulk plant RNA-seq. | Pathogen-focused transcriptome. |
| PRJNA1263482 — TSWV D-RB | Unresolved | The accession-to-paper link, source tissue, and raw-read library type have not been verified. | Pending archive and methods check. |
| PRJNA590139 — Polish plant-virus HTS survey | Meets, field subset | \[Min20b\] used bulk rRNA-depleted total RNA from field-collected symptomatic plant leaves; some additional libraries came from mechanically inoculated test plants. | PRJNA590139; count only original field-plant libraries. |

## Rules applied in this audit

- A paper’s claim that sequencing occurred is not enough: the raw-read accession must resolve to the relevant plant dataset.
- GenBank accessions for assembled pathogen genomes do not satisfy the raw-read archive criterion.
- Microarrays, targeted assays, and pathogen-only sequencing do not qualify as bulk RNA-seq from plant tissue.
- Multiple papers using the same field RNA-seq dataset are linked to avoid inflating the number of independent field collections.

---

## References

\[Wam18\] M. J. Wamaitha *et al.*, “Metagenomic analysis of viruses associated with maize lethal necrosis in Kenya,” *Virology Journal*, vol. 15, May 2018, doi: [10.1186/s12985-018-0999-2](https://doi.org/10.1186/s12985-018-0999-2).

\[Elm22\] M. G. Elmore *et al.*, “Detection and discovery of plant viruses in soybean by metagenomic sequencing,” *Virology Journal*, vol. 19, Sep. 2022, doi: [10.1186/s12985-022-01872-5](https://doi.org/10.1186/s12985-022-01872-5).

\[Gaa20\] Y. Gaafar *et al.*, “Investigating the Pea Virome in Germany—Old Friends and New Players in the Field(s),” *Frontiers in Microbiology*, vol. 11, Nov. 2020, doi: [10.3389/fmicb.2020.583242](https://doi.org/10.3389/fmicb.2020.583242).

\[San13c\] S. D. Santo *et al.*, “The plasticity of the grapevine berry transcriptome,” *Genome Biology*, vol. 14, pp. r54–r54, Jun. 2013, doi: [10.1186/gb-2013-14-6-r54](https://doi.org/10.1186/gb-2013-14-6-r54).

\[Nag12\] A. Nagano *et al.*, “Deciphering and prediction of transcriptome dynamics under fluctuating field conditions.” *Cell*, vol. 151 6, pp. 1358–69, Dec. 2012, doi: [10.1016/j.cell.2012.10.048](https://doi.org/10.1016/j.cell.2012.10.048).

\[Nag19\] A. Nagano, T. Kawagoe, J. Sugisaka, M. Honjo, K. Iwayama, and H. Kudoh, “Annual transcriptome dynamics in natural environments reveals plant seasonal adaptation,” *Nature Plants*, vol. 5, pp. 74–83, Jan. 2019, doi: [10.1038/s41477-018-0338-z](https://doi.org/10.1038/s41477-018-0338-z).

\[Sch25d\] E. Schnabel, C. A. D. Xavier, A. Whitfield, Z. Dubrow, G. M. Pham, and E. Cieniewicz, “Exploring the Virome of Blackberry and Wild Rubus spp. in South Carolina,” *Phytobiomes journal*, vol. 9, pp. 80–94, Feb. 2025, doi: [10.1094/pbiomes-11-24-0106-r](https://doi.org/10.1094/pbiomes-11-24-0106-r).

\[Han16\] J. Han, L. L. Domier, B. J. Cassone, A. Dorrance, and F. Qu, “Assessment of common soybean-infecting viruses in Ohio, USA, through multi-site sampling and high-throughput sequencing,” *Plant Health Progress*, vol. 17, pp. 133–140, 2016, doi: [10.1094/PHP-RS-16-0018](https://doi.org/10.1094/PHP-RS-16-0018).

\[Kam18\] M. Kamitani, A. Nagano, M. Honjo, and H. Kudoh, “A Survey on Plant Viruses in Natural Brassicaceae Communities Using RNA-Seq,” *Microbial Ecology*, vol. 78, pp. 113–121, Oct. 2018, doi: [10.1007/s00248-018-1271-4](https://doi.org/10.1007/s00248-018-1271-4).

\[Red21b\] C. D. Redila, V. Prakash, and S. Nouri, “Metagenomics Analysis of the Wheat Virome Identifies Novel Plant and Fungal-Associated Viral Sequences,” *Viruses*, vol. 13, Dec. 2021, doi: [10.3390/v13122457](https://doi.org/10.3390/v13122457).

\[Mut18\] J. Mutuku *et al.*, “Metagenomic Analysis of Plant Virus Occurrence in Common Bean (Phaseolus vulgaris) in Central Kenya,” *Frontiers in Microbiology*, vol. 9, Dec. 2018, doi: [10.3389/fmicb.2018.02939](https://doi.org/10.3389/fmicb.2018.02939).

\[Mat17b\] E. Matsumura *et al.*, “Deep Sequencing Analysis of RNAs from Citrus Plants Grown in a Citrus Sudden Death-Affected Area Reveals Diverse Known and Putative Novel Viruses,” *Viruses*, vol. 9, Apr. 2017, doi: [10.3390/v9040092](https://doi.org/10.3390/v9040092).

\[Hon19b\] M. Honjo *et al.*, “Seasonality of interactions between a plant virus and its host during persistent infection in a natural environment,” *The ISME Journal*, vol. 14, pp. 506–518, Oct. 2019, doi: [10.1038/s41396-019-0519-4](https://doi.org/10.1038/s41396-019-0519-4).

\[Mul22\] R. M. Mulenga *et al.*, “Survey for virus diversity in common bean (Phaseolus vulgaris L.) fields and the detection of a novel strain of cowpea polerovirus 1 in Zambia.” *Plant disease*, Feb. 2022, doi: [10.1094/PDIS-11-21-2533-RE](https://doi.org/10.1094/PDIS-11-21-2533-RE).

\[Obo24\] D. Obonyo, G. Ouma, R. Ikawa, and D. Odeny, “Meta-transcriptomic identification of groundnut RNA viruses in western Kenya and the novel detection of groundnut as a host for Cauliflower mosaic virus.” *Virology*, vol. 593, pp. 110011, Feb. 2024, doi: [10.1016/j.virol.2024.110011](https://doi.org/10.1016/j.virol.2024.110011).

\[Ada21\] T. Adams *et al.*, “Rust expression browser: an open source database for simultaneous analysis of host and pathogen gene expression profiles with expVIP,” *BMC Genomics*, vol. 22, Mar. 2021, doi: [10.1186/s12864-021-07488-3](https://doi.org/10.1186/s12864-021-07488-3).

\[Sin26b\] R. Singh, A. J. Muhich, C. Tom, C. Caseys, and D. Kliebenstein, “A multiplant transcriptomic atlas reveals conserved and lineage-specific defense architectures in response to Botrytis cinerea,” *Proceedings of the National Academy of Sciences of the United States of America*, vol. 123, Jan. 2026, doi: [10.64898/2026.01.14.699558](https://doi.org/10.64898/2026.01.14.699558).

\[Rad19\] G. V. Radhakrishnan *et al.*, “MARPLE, a point-of-care, strain-level disease diagnostics and surveillance tool for complex fungal pathogens,” *BMC Biology*, vol. 17, Aug. 2019, doi: [10.1186/s12915-019-0684-y](https://doi.org/10.1186/s12915-019-0684-y).

\[Bue17\] V. Bueno-Sancho *et al.*, “Pathogenomic Analysis of Wheat Yellow Rust Lineages Detects Seasonal Variation and Host Specificity,” *Genome Biology and Evolution*, vol. 9, pp. 3282–3296, Nov. 2017, doi: [10.1093/gbe/evx241](https://doi.org/10.1093/gbe/evx241).

\[Tsu25\] A. Tsushima *et al.*, “k-mer-based GWAS reveals a candidate avirulence gene and structural variation in Puccinia triticina linked to gain of Lr20 virulence,” *BMC Genomics*, vol. 26, Nov. 2025, doi: [10.1186/s12864-025-12230-4](https://doi.org/10.1186/s12864-025-12230-4).

\[Li20e\] Y. Li *et al.*, “Magnaporthe oryzae Auxiliary Activity Protein MoAa91 Functions as Chitin-Binding Protein To Induce Appressorium Formation on Artificial Inductive Surfaces and Suppress Plant Immunity,” *mBio*, vol. 11, Mar. 2020, doi: [10.1128/mBio.03304-19](https://doi.org/10.1128/mBio.03304-19).

\[Hub15\] A. Hubbard *et al.*, “Field pathogenomics reveals the emergence of a diverse wheat yellow rust population,” *Genome Biology*, vol. 16, Feb. 2015, doi: [10.1186/s13059-015-0590-8](https://doi.org/10.1186/s13059-015-0590-8).

\[Cai21b\] M. Caicedo, E. D. Munaiz, R. Malvar, J. C. Jimenez, and B. Ordás, “Precision Mapping of a Maize MAGIC Population Identified a Candidate Gene for the Senescence-Associated Physiological Traits,” *Frontiers in Genetics*, vol. 12, Oct. 2021, doi: [10.3389/fgene.2021.716821](https://doi.org/10.3389/fgene.2021.716821).

\[Dil19\] T. Dilks, K. Halsey, R. P. D. Vos, K. Hammond-Kosack, and N. Brown, “Non-canonical fungal G-protein coupled receptors promote Fusarium head blight on wheat,” *PLoS Pathogens*, vol. 15, Apr. 2019, doi: [10.1371/journal.ppat.1007666](https://doi.org/10.1371/journal.ppat.1007666).

\[Spa16\] H. Sonah, X. Zhang, R. Deshmukh, M. H. Borhan, W. D. D. Fernando, and R. Bélanger, “Comparative Transcriptomic Analysis of Virulence Factors in Leptosphaeria maculans during Compatible and Incompatible Interactions with Canola,” *Frontiers in Plant Science*, vol. 7, Dec. 2016, doi: [10.3389/fpls.2016.01784](https://doi.org/10.3389/fpls.2016.01784).

\[Hol23c\] S. Holden *et al.*, “Uncovering the history of recombination and population structure in western Canadian stripe rust populations through mating type alleles,” *BMC Biology*, vol. 21, Apr. 2023, doi: [10.1186/s12915-023-01717-9](https://doi.org/10.1186/s12915-023-01717-9).

\[Gup23\] Y. Gupta *et al.*, “Major proliferation of transposable elements shaped the genome of the soybean rust pathogen Phakopsora pachyrhizi,” *Nature Communications*, vol. 14, Apr. 2023, doi: [10.1038/s41467-023-37551-4](https://doi.org/10.1038/s41467-023-37551-4).

\[Tsu22\] A. Tsushima, C. Lewis, K. Flath, S. Kildea, and D. G. O. Saunders, “Wheat stem rust recorded for the first time in decades in Ireland,” *Plant Pathology*, vol. 71, pp. 890–900, Jan. 2022, doi: [10.1111/ppa.13532](https://doi.org/10.1111/ppa.13532).

\[Roc24\] F. Rocher *et al.*, “Integrative systems biology of wheat susceptibility to Fusarium graminearum uncovers a conserved gene regulatory network and identifies master regulators targeted by fungal core effectors,” *BMC Biology*, vol. 22, Mar. 2024, doi: [10.1186/s12915-024-01852-x](https://doi.org/10.1186/s12915-024-01852-x).

\[Pog22\] M. Pogány *et al.*, “Redox and Hormonal Changes in the Transcriptome of Grape (Vitis vinifera) Berries during Natural Noble Rot Development,” *Plants*, vol. 11, Mar. 2022, doi: [10.3390/plants11070864](https://doi.org/10.3390/plants11070864).

\[Fow21\] A. Fowkes *et al.*, “Integrating High throughput Sequencing into Survey Design Reveals Turnip Yellows Virus and Soybean Dwarf Virus in Pea (Pisum Sativum) in the United Kingdom,” *Viruses*, vol. 13, Dec. 2021, doi: [10.3390/v13122530](https://doi.org/10.3390/v13122530).

\[Cho26c\] H. Choi *et al.*, “Library Preparation Biases Plant Virome Detection: Poly(A) mRNA Enrichment vs. rRNA Depletion in Pepper and Garlic,” *International Journal of Molecular Sciences*, vol. 27, Feb. 2026, doi: [10.3390/ijms27052300](https://doi.org/10.3390/ijms27052300).

\[Cho20d\] H. Choi *et al.*, “Identification of Viruses and Viroids Infecting Tomato and Pepper Plants in Vietnam by Metatranscriptomics,” *International Journal of Molecular Sciences*, vol. 21, Oct. 2020, doi: [10.3390/ijms21207565](https://doi.org/10.3390/ijms21207565).

\[Dan24\] W. Dantes, L. Boatwright, and E. Cieniewicz, “Comparing RT-PCR of individual samples with high throughput sequencing of pooled plant samples for field-level surveillance of viruses in blackberry and wild Rubus.” *Plant disease*, Apr. 2024, doi: [10.1094/PDIS-11-23-2428-RE](https://doi.org/10.1094/PDIS-11-23-2428-RE).

\[Rum19\] A. Rumbou *et al.*, “Unravelling the virome in birch: RNA-Seq reveals a complex of known and novel viruses,” *PLoS ONE*, vol. 15, Aug. 2019, doi: [10.1371/journal.pone.0221834](https://doi.org/10.1371/journal.pone.0221834).

\[Ram22b\] J.-A. Ramírez-Pool *et al.*, “Transcriptomic Analysis of the Host Response to Mild and Severe CTV Strains in Naturally Infected Citrus sinensis Orchards,” *International Journal of Molecular Sciences*, vol. 23, Feb. 2022, doi: [10.3390/ijms23052435](https://doi.org/10.3390/ijms23052435).

\[Dia23\] N. Dias *et al.*, “Viromes of field-grown tomatoes and peppers in Tennessee revealed by RNA sequencing followed by bioinformatic analysis,” *Plant Health Progress*, Jan. 2023, doi: [10.1094/php-10-22-0107-rs](https://doi.org/10.1094/php-10-22-0107-rs).

\[Xia22e\] H.-G. Xiao, W. Hao, G. Storoschuk, J. MacDonald, and H. Sanfaçon, “Characterizing the Virome of Apple Orchards Affected by Rapid Decline in the Okanagan and Similkameen Valleys of British Columbia (Canada),” *Pathogens*, vol. 11, Oct. 2022, doi: [10.3390/pathogens11111231](https://doi.org/10.3390/pathogens11111231).

\[Kim21g\] N.-Y. Kim, H.-J. Lee, H. Kim, S.-H. Lee, J. Moon, and R. Jeong, “Identification of Plant Viruses Infecting Pear Using RNA Sequencing,” *The Plant Pathology Journal*, vol. 37, pp. 258–267, Jun. 2021, doi: [10.5423/PPJ.OA.01.2021.0009](https://doi.org/10.5423/PPJ.OA.01.2021.0009).

\[Hu21c\] R.-B. Hu *et al.*, “Cultivated and wild grapevines in Tennessee possess overlapping but distinct virus populations.” *Plant disease*, Feb. 2021, doi: [10.1094/PDIS-11-20-2483-SC](https://doi.org/10.1094/PDIS-11-20-2483-SC).

\[Hod20\] B. A. Hodge, P. Paul, and L. Stewart, “Occurrence and High-Throughput Sequencing of Viruses in Ohio Wheat.” *Plant disease*, pp. PDIS08191724RE, Apr. 2020, doi: [10.1094/pdis-08-19-1724-re](https://doi.org/10.1094/pdis-08-19-1724-re).

\[Jo18\] Y. Jo, J. Bae, S.-M. Kim, H. Choi, B.-C. Lee, and W. Cho, “Barley RNA viromes in six different geographical regions in Korea,” *Scientific Reports*, vol. 8, Sep. 2018, doi: [10.1038/s41598-018-31671-4](https://doi.org/10.1038/s41598-018-31671-4).

\[Nhl18\] T. F. Nhlapo, D. J. G. Rees, D. Odeny, J. Mulabisana, and M. Rey, “Viral metagenomics reveals sweet potato virus diversity in the Eastern and Western Cape provinces of South Africa,” Jul. 01, 2018. doi: [10.1016/J.SAJB.2018.05.024](https://doi.org/10.1016/J.SAJB.2018.05.024).

\[Esc26\] C. Escalante *et al.*, “Metatranscriptomics analysis reveals the cotton virome in the southern United States,” *Scientific Reports*, vol. 16, Feb. 2026, doi: [10.1038/s41598-026-40828-5](https://doi.org/10.1038/s41598-026-40828-5).

\[Wri20\] A. Wright, A. Cross, and S. Harper, “A bushel of viruses: Identification of seventeen novel putative viruses by RNA-seq in six apple trees,” *PLoS ONE*, vol. 15, Jan. 2020, doi: [10.1371/journal.pone.0227669](https://doi.org/10.1371/journal.pone.0227669).

\[Jo18b\] Y. Jo *et al.*, “Peach RNA viromes in six different peach cultivars,” *Scientific Reports*, vol. 8, Jan. 2018, doi: [10.1038/s41598-018-20256-w](https://doi.org/10.1038/s41598-018-20256-w).

\[Von21\] A. M. Vondras *et al.*, “Rootstock influences the effect of grapevine leafroll‐associated viruses on berry development and metabolism via abscisic acid signalling,” *Molecular Plant Pathology*, vol. 22, pp. 984–1005, Jun. 2021, doi: [10.1111/mpp.13077](https://doi.org/10.1111/mpp.13077).

\[Pag20\] C. Pagliarani *et al.*, “Molecular memory of Flavescence dorée phytoplasma in recovering grapevines,” *Horticulture Research*, vol. 7, Aug. 2020, doi: [10.1038/s41438-020-00348-3](https://doi.org/10.1038/s41438-020-00348-3).

\[Bel23c\] D. Belkina, D. Karpova, E. Porotikova, I. Lifanov, and S. Vinogradova, “Grapevine Virome of the Don Ampelographic Collection in Russia Has Concealed Five Novel Viruses,” *Viruses*, vol. 15, Dec. 2023, doi: [10.3390/v15122429](https://doi.org/10.3390/v15122429).

\[Jo20\] Y. Jo, S.-M. Kim, H. Choi, J.-W. Yang, B.-C. Lee, and W. Cho, “Sweet potato viromes in eight different geographical regions in Korea and two different cultivars,” *Scientific Reports*, vol. 10, Feb. 2020, doi: [10.1038/s41598-020-59518-x](https://doi.org/10.1038/s41598-020-59518-x).

\[Lov18\] J. Lovell *et al.*, “The genomic landscape of molecular responses to natural drought stress in Panicum hallii,” *Nature Communications*, vol. 9, Dec. 2018, doi: [10.1038/s41467-018-07669-x](https://doi.org/10.1038/s41467-018-07669-x).

\[Oco12\] R. O’Connell *et al.*, “Lifestyle transitions in plant pathogenic Colletotrichum fungi deciphered by genome and transcriptome analyses,” *Nature Genetics*, vol. 44, pp. 1060–1065, Aug. 2012, doi: [10.1038/ng.2372](https://doi.org/10.1038/ng.2372).

\[Car26\] R. D. Carvalho *et al.*, “Virulence gains in the Puccinia striiformis f. sp. tritici PstS10 lineage correlate with expression polymorphism in a candidate Avr effector,” *Communications Biology*, vol. 9, Apr. 2026, doi: [10.1038/s42003-026-10018-0](https://doi.org/10.1038/s42003-026-10018-0).

\[Li22g\] X. Li *et al.*, “Integrative transcriptome analysis revealed the pathogenic molecular basis of Rhizoctonia solani AG-3 TB at three progressive stages of infection,” *Frontiers in Microbiology*, vol. 13, Oct. 2022, doi: [10.3389/fmicb.2022.1001327](https://doi.org/10.3389/fmicb.2022.1001327).

\[Sek19\] R. Sekhon *et al.*, “Integrated Genome-Scale Analysis Identifies Novel Genes and Networks Underlying Senescence in Maize\[OPEN\],” *Plant Cell*, vol. 31, pp. 1968–1989, Jun. 2019, doi: [10.1105/tpc.18.00930](https://doi.org/10.1105/tpc.18.00930).

\[Rud15\] J. Rudd *et al.*, “Transcriptome and Metabolite Profiling of the Infection Cycle of Zymoseptoria tritici on Wheat Reveals a Biphasic Interaction with Plant Immunity Involving Differential Pathogen Chromosomal Contributions and a Variation on the Hemibiotrophic Lifestyle Definition1\[OPEN\],” *Plant Physiology*, vol. 167, pp. 1158–1185, Jan. 2015, doi: [10.1104/pp.114.255927](https://doi.org/10.1104/pp.114.255927).

\[Tom26\] A. Tomita *et al.*, “Transcriptomic biomarkers reveal jasmonic and salicylic acid state under field herbivory,” Feb. 12, 2026. doi: [10.1101/2025.05.29.656841](https://doi.org/10.1101/2025.05.29.656841).

\[Wu25c\] L.-J. Wu, J. Wang, S. Shen, Z. Yang, and X. Hu, “Transcriptomic analysis of two Chinese wheat landraces with contrasting Fusarium head blight resistance reveals miRNA-mediated defense mechanisms,” *Frontiers in Plant Science*, vol. 16, Feb. 2025, doi: [10.3389/fpls.2025.1537605](https://doi.org/10.3389/fpls.2025.1537605).

\[Pou19\] R. S. Poudel, J. Richards, S. Shrestha, S. Solanki, and R. Brueggeman, “Transcriptome-wide association study identifies putative elicitors/suppressor of Puccinia graminis f. sp. tritici that modulate barley rpg4-mediated stem rust resistance,” *BMC Genomics*, vol. 20, Dec. 2019, doi: [10.1186/s12864-019-6369-7](https://doi.org/10.1186/s12864-019-6369-7).

\[Mje26\] E. Y. Mjema, M. L. Bonatelli, and S. Laubinger, “Molecular and phenotypic footprints of climate in native Arabidopsis thaliana,” Mar. 04, 2026. doi: [10.64898/2026.03.02.709013](https://doi.org/10.64898/2026.03.02.709013).

\[Dia22b\] N. Dias, R.-B. Hu, D. Hensley, Z. Hansen, L. Domier, and M. R. Hajimorad, “A Survey for Viruses and Viroids of Peach in Tennessee Orchards by RNA Sequencing,” *Plant Health Progress*, Apr. 2022, doi: [10.1094/php-01-22-0007-sc](https://doi.org/10.1094/php-01-22-0007-sc).

\[Kha25b\] M. F. Khatun, M. Kwak, M. Kwon, M. M. Hossain, and E.-J. Kil, “New insights into viral threats in soybean (Glycine max) crops from Bangladesh, including a novel crinivirus,” *Frontiers in Microbiology*, vol. 16, Feb. 2025, doi: [10.3389/fmicb.2025.1523767](https://doi.org/10.3389/fmicb.2025.1523767).

\[Ber26\] M. M. J. Berger *et al.*, “Esca disease triggers local transcriptomic response and DNA methylation changes in grapevine,” *Journal of Experimental Botany*, vol. 77, pp. 5244–5266, Jun. 2026, doi: [10.1093/jxb/erag192](https://doi.org/10.1093/jxb/erag192).

\[Nab22\] S. Nabi *et al.*, “High-Throughput RNA Sequencing of Mosaic Infected and Non-Infected Apple (Malus × domestica Borkh.) Cultivars: From Detection to the Reconstruction of Whole Genome of Viruses and Viroid,” *Plants*, vol. 11, Mar. 2022, doi: [10.3390/plants11050675](https://doi.org/10.3390/plants11050675).

\[Rod26\] L. Rodríguez-Robles, P. J. Martínez-García, P. Martínez-Gómez, and M. Rubio, “Detection and Characterization of Plum Pox Virus (Potyvirus plumpoxi) Marcus Strains in Spanish Apricot and Peach Orchards Through RNA-Seq Analysis,” *Agronomy*, Mar. 2026, doi: [10.3390/agronomy16060608](https://doi.org/10.3390/agronomy16060608).

\[Xu19\] Y. Xu, S. Li, C.-Y. Na, L. Yang, and M.-G. Lu, “Analyses of virus/viroid communities in nectarine trees by next-generation sequencing and insight into viral synergisms implication in host disease symptoms,” *Scientific Reports*, vol. 9, Aug. 2019, doi: [10.1038/s41598-019-48714-z](https://doi.org/10.1038/s41598-019-48714-z).

\[Kim24c\] H.-J. Kim, S.-R. Choi, I.-S. Cho, and R. Jeong, “Viral Metatranscriptomic Analysis to Reveal the Diversity of Viruses Infecting Satsuma Mandarin (Citrus unshiu) in Korea,” *The Plant Pathology Journal*, vol. 40, pp. 115–124, Apr. 2024, doi: [10.5423/PPJ.OA.01.2024.0009](https://doi.org/10.5423/PPJ.OA.01.2024.0009).

\[Sur26\] A. S. Suryaningsih *et al.*, “Metatranscriptomic profiling uncovers a distinctive pepper virome in tropical Indonesia,” *Frontiers in Plant Science*, vol. 17, Sep. 2026, doi: [10.3389/fpls.2026.1933649](https://doi.org/10.3389/fpls.2026.1933649).

\[Kus20\] S. K. Kushwaha, R. Vetukuri, F. Odilbekov, N. Pareek, T. Henriksson, and A. Chawade, “Differential Gene Expression Analysis of Wheat Breeding Lines Reveal Molecular Insights in Yellow Rust Resistance under Field Conditions,” Nov. 29, 2020. doi: [10.3390/agronomy10121888](https://doi.org/10.3390/agronomy10121888).

\[Lia21\] Y. Liang, R. Tabien, L. Tarpley, A. R. Mohammed, and E. Septiningsih, “Transcriptome profiling of two rice genotypes under mild field drought stress during grain-filling stage,” *AoB Plants*, vol. 13, Jul. 2021, doi: [10.1093/aobpla/plab043](https://doi.org/10.1093/aobpla/plab043).

\[Wu23d\] Q. Wu *et al.*, “A Metagenomic Investigation of the Viruses Associated with Shiraz Disease in Australia,” *Viruses*, vol. 15, Mar. 2023, doi: [10.3390/v15030774](https://doi.org/10.3390/v15030774).

\[Jo20c\] Y. Jo, H. Choi, S. Lian, J. K. Cho, H. Chu, and W. Cho, “Identification of viruses infecting six plum cultivars in Korea by RNA-sequencing,” *PeerJ*, vol. 8, Jul. 2020, doi: [10.7717/peerj.9588](https://doi.org/10.7717/peerj.9588).

\[San16c\] S. D. Santo *et al.*, “Plasticity of the Berry Ripening Program in a White Grape Variety,” *Frontiers in Plant Science*, vol. 7, Jul. 2016, doi: [10.3389/fpls.2016.00970](https://doi.org/10.3389/fpls.2016.00970).

\[San25b\] V. M. Santi *et al.*, “Molecular characterization and survey of viruses and viroids infecting apple and pear trees in southern Brazil,” *Tropical Plant Pathology*, vol. 50, Jul. 2025, doi: [10.1007/s40858-025-00747-8](https://doi.org/10.1007/s40858-025-00747-8).

\[Wil16d\] O. Wilkins *et al.*, “EGRINs (Environmental Gene Regulatory Influence Networks) in Rice That Function in the Response to Water Deficit, High Temperature, and Agricultural Environments\[OPEN\],” *Plant Cell*, vol. 28, pp. 2365–2384, Sep. 2016, doi: [10.1105/tpc.16.00158](https://doi.org/10.1105/tpc.16.00158).

\[Ple15\] A. Plessis *et al.*, “Multiple abiotic stimuli are integrated in the regulation of rice gene expression under field conditions,” *eLife*, vol. 4, Nov. 2015, doi: [10.7554/eLife.08411](https://doi.org/10.7554/eLife.08411).

\[Lee23h\] H.-J. Lee, S.-M. Kim, and R. Jeong, “Analysis of Wheat Virome in Korea Using Illumina and Oxford Nanopore Sequencing Platforms,” *Plants*, vol. 12, Jun. 2023, doi: [10.3390/plants12122374](https://doi.org/10.3390/plants12122374).

\[Rat22b\] N. Rathore, P. Kumar, N. Mehta, M. Swarnkar, R. Shankar, and A. Chawla, “Time-series RNA-Seq transcriptome profiling reveals novel insights about cold acclimation and de-acclimation processes in an evergreen shrub of high altitude,” *Scientific Reports*, vol. 12, Sep. 2022, doi: [10.1038/s41598-022-19834-w](https://doi.org/10.1038/s41598-022-19834-w).

\[Jo21\] Y. Jo *et al.*, “Comparative Study of Metagenomics and Metatranscriptomics to Reveal Microbiomes in Overwintering Pepper Fruits,” *International Journal of Molecular Sciences*, vol. 22, Jun. 2021, doi: [10.3390/ijms22126202](https://doi.org/10.3390/ijms22126202).

\[Maa21\] A. Maachi, C. Torre, R. Sempere, Y. Hernando, M. A. Aranda, and L. Donaire, “Use of High-Throughput Sequencing and Two RNA Input Methods to Identify Viruses Infecting Tomato Crops,” *Microorganisms*, vol. 9, May 2021, doi: [10.3390/microorganisms9051043](https://doi.org/10.3390/microorganisms9051043).

\[Zha22k\] R. Zhang *et al.*, “Response of Multiple Tissues to Drought Revealed by a Weighted Gene Co-Expression Network Analysis in Foxtail Millet \[Setaria italica (L.) P. Beauv.\],” *Frontiers in Plant Science*, vol. 12, Jan. 2022, doi: [10.3389/fpls.2021.746166](https://doi.org/10.3389/fpls.2021.746166).

\[Xia19\] H.-G. Xiao, C. Li, M. A. Rwahnih, V. Dolja, and B. Meng, “Metagenomic Analysis of Riesling Grapevine Reveals a Complex Virome Including Two New and Divergent Variants of Grapevine leafroll-associated virus 3.” *Plant disease*, vol. 103 6, pp. 1275–1285, Apr. 2019, doi: [10.1094/pdis-09-18-1503-re](https://doi.org/10.1094/pdis-09-18-1503-re).

\[Cho16b\] I. Cho *et al.*, “Deep Sequencing Analysis of Apple Infecting Viruses in Korea,” *The Plant Pathology Journal*, vol. 32, pp. 441–451, Oct. 2016, doi: [10.5423/PPJ.OA.04.2016.0104](https://doi.org/10.5423/PPJ.OA.04.2016.0104).

\[Goz25\] M. Gozmanova *et al.*, “Viral and Viroid Communities in Peach Cultivars Grown in Bulgaria,” *Horticulturae*, May 2025, doi: [10.3390/horticulturae11050503](https://doi.org/10.3390/horticulturae11050503).

\[Zha20r\] S. Zhang *et al.*, “Virome of Camellia japonica: Discovery of and Molecular Characterization of New Viruses of Different Taxa in Camellias,” *Frontiers in Microbiology*, vol. 11, May 2020, doi: [10.3389/fmicb.2020.00945](https://doi.org/10.3389/fmicb.2020.00945).

\[Ye23\] L.-J. Ye *et al.*, “Variation in gene expression along an elevation gradient of Rhododendron sanguineum var. haemaleum assessed in a comparative transcriptomic analysis,” *Frontiers in Plant Science*, vol. 14, Mar. 2023, doi: [10.3389/fpls.2023.1133065](https://doi.org/10.3389/fpls.2023.1133065).

\[Pat25c\] M. Patanita *et al.*, “Transcriptome profiling of symptomatic vs. asymptomatic grapevine plants reveals candidate genes for plant improvement against trunk diseases,” *BMC Plant Biology*, vol. 25, Jul. 2025, doi: [10.1186/s12870-025-06763-9](https://doi.org/10.1186/s12870-025-06763-9).

\[Die23\] S. Diez-Hermano *et al.*, “Diversity of RNA Viruses in Declining Mediterranean Forests,” *Microorganisms*, vol. 14, Oct. 2023, doi: [10.3390/microorganisms14071445](https://doi.org/10.3390/microorganisms14071445).

\[Amo26\] S. S. Amoia, A. Giampetruzzi, F. F. de Sousa, L. F. António, A. T. P. da Cunha, and A. Minafra, “Plant Viral Metagenomic Analysis from a Preliminary Field Survey in Angola Reveals Complex Mixed Infections in Vegetable Crops,” *Viruses*, vol. 18, Jul. 2026, doi: [10.3390/v18080822](https://doi.org/10.3390/v18080822).

\[Mur25\] T. Muranaka, G. Yumoto, M. Honjo, A. J. Nagano, J. Zhou, and H. Kudoh, “Coincidence of the threshold temperature of seasonal switching for diel transcriptomic oscillations and growth,” *Plant and Cell Physiology*, vol. 66, pp. 1412–1425, Aug. 2025, doi: [10.1093/pcp/pcaf092](https://doi.org/10.1093/pcp/pcaf092).

\[Pao21\] M. Paolinelli, G. Escoriaza, C. Césari, S. García-Lampasona, and R. Hernández-Martínez, “Characterization of Grapevine Wood Microbiome Through a Metatranscriptomic Approach,” *Microbial Ecology*, vol. 83, pp. 658–668, Jun. 2021, doi: [10.1007/s00248-021-01801-z](https://doi.org/10.1007/s00248-021-01801-z).

\[Deb25b\] H. Debat, M. Paolinelli, G. Escoriaza, S. García-Lampasona, S. Gomez-Talquenca, and N. Bejerman, “Grapevine holobiome metatranscriptomics provides a glimpse into the wood mycovirome,” *bioRxiv*, Mar. 2025, doi: [10.1101/2025.03.21.644598](https://doi.org/10.1101/2025.03.21.644598).

\[Ras22b\] S. Rashid, F. Wani, G. Ali, T. Sofi, Z. Dar, and A. Hamid, “Viral metatranscriptomic approach to study the diversity of virus(es) associated with Common Bean (Phaseolus vulgaris L.) in the North-Western Himalayan region of India,” *Frontiers in Microbiology*, vol. 13, Sep. 2022, doi: [10.3389/fmicb.2022.943382](https://doi.org/10.3389/fmicb.2022.943382).

\[Har21e\] Z. N. Harris *et al.*, “Multi-dimensional leaf phenotypes reflect root system genotype in grafted grapevine over the growing season,” *GigaScience*, vol. 10, Dec. 2021, doi: [10.1093/gigascience/giab087](https://doi.org/10.1093/gigascience/giab087).

\[Har23c\] Z. N. Harris *et al.*, “Grapevine scion gene expression is driven by rootstock and environment interaction,” *BMC Plant Biology*, vol. 23, Jan. 2023, doi: [10.1186/s12870-023-04223-w](https://doi.org/10.1186/s12870-023-04223-w).

\[Nan25c\] N. Nancarrow, B. Rodoni, S.-K. Lam, W. M. Kinoti, and P. Trebicki, “Hidden diversity of yellow dwarf viruses in Australian cereals,” *Archives of Virology*, vol. 171, Nov. 2025, doi: [10.1007/s00705-025-06451-x](https://doi.org/10.1007/s00705-025-06451-x).

\[Ren22g\] Z.-W. Ren, P. Kopittke, F.-J. Zhao, and P. Wang, “Nutrient accumulation and transcriptome patterns during grain development in rice,” *Journal of Experimental Botany*, vol. 74, pp. 909–930, Oct. 2022, doi: [10.1093/jxb/erac426](https://doi.org/10.1093/jxb/erac426).

\[Lov16\] J. Lovell *et al.*, “Promises and Challenges of Eco-Physiological Genomics in the Field: Tests of Drought Responses in Switchgrass1\[OPEN\],” *Plant Physiology*, vol. 172, pp. 734–748, May 2016, doi: [10.1104/pp.16.00545](https://doi.org/10.1104/pp.16.00545).

\[Bon19\] L. Bono *et al.*, “Spatiotemporal dynamics of RNA virus diversity in a phyllosphere microbial community,” Sep. 18, 2019. doi: [10.1101/772475](https://doi.org/10.1101/772475).

\[Rog23c\] O. Rogier *et al.*, “RNAseq based variant dataset in a black poplar association panel,” *BMC Research Notes*, vol. 16, Oct. 2023, doi: [10.1186/s13104-023-06521-w](https://doi.org/10.1186/s13104-023-06521-w).

\[Man25e\] Z. Mangral *et al.*, “Transcriptomic Insights Into Elevation‐Dependent Gene Expression in Rhododendron anthopogon D.Don: Implications for Climate Resilience,” *Physiologia Plantarum*, vol. 177, Jul. 2025, doi: [10.1111/ppl.70419](https://doi.org/10.1111/ppl.70419).

\[Has21c\] Y. Hashida *et al.*, “Fillable and unfillable gaps in plant transcriptome under field and controlled environments,” *Plant, Cell & Environment*, vol. 45, pp. 2410–2427, Aug. 2021, doi: [10.1111/pce.14367](https://doi.org/10.1111/pce.14367).

\[Gla20\] M. Glasa, L. Predajňa, N. Sihelská, K. Šoltys, and A. Ruiz-García, “Analysis of Virome by High-Throughput Sequencing Revealed Multiple Infection and Intra-Virus Diversity in a Single Grapevine Plant,” May 01, 2020. doi: [10.2478/ahr-2020-0009](https://doi.org/10.2478/ahr-2020-0009).

\[Zag26\] S. Zagier, O. Alisawi, W. Aljuaifari, and F. A. Al-fadhal, “Metatranscriptomic analysis reveals two new viruses and highly expressed endogenous elements in pepper in Iraq,” *Jurnal Hama dan Penyakit Tumbuhan Tropika*, Jul. 2026, doi: [10.23960/jhptt.226463-472](https://doi.org/10.23960/jhptt.226463-472).

\[Min25\] D. Min, Y. Jo, J. Park, G.-G. Min, J.-S. Hong, and W. Cho, “Comprehensive Virome Analysis of Commercial Lilies in South Korea by RT-PCR, High-Throughput Sequencing, and Phylogenetic Analyses,” *International Journal of Molecular Sciences*, vol. 26, Oct. 2025, doi: [10.3390/ijms26199598](https://doi.org/10.3390/ijms26199598).

\[Fas12b\] M. Fasoli *et al.*, “The Grapevine Expression Atlas Reveals a Deep Transcriptome Shift Driving the Entire Plant into a Maturation Program\[W\]\[OA\],” *Plant Cell*, vol. 24, pp. 3489–3505, Sep. 2012, doi: [10.1105/tpc.112.100230](https://doi.org/10.1105/tpc.112.100230).

\[Bel24\] M. T. Belete *et al.*, “Deciphering the virome of Chunkung (Cnidium officinale) showing dwarfism-like symptoms via a high-throughput sequencing analysis,” *Virology Journal*, vol. 21, Apr. 2024, doi: [10.1186/s12985-024-02361-7](https://doi.org/10.1186/s12985-024-02361-7).

\[Mar16c\] S. Marzano and L. Domier, “Novel mycoviruses discovered from metatranscriptomics survey of soybean phyllosphere phytobiomes.” *Virus research*, vol. 213, pp. 332–342, Feb. 2016, doi: [10.1016/j.virusres.2015.11.002](https://doi.org/10.1016/j.virusres.2015.11.002).

\[Van22c\] A. VanWallendael *et al.*, “Host genotype controls ecological change in the leaf fungal microbiome,” *PLoS Biology*, vol. 20, Aug. 2022, doi: [10.1371/journal.pbio.3001681](https://doi.org/10.1371/journal.pbio.3001681).

\[Cho23c\] H. Choi *et al.*, “Investigating Variability in Viral Presence and Abundance across Soybean Seed Development Stages Using Transcriptome Analysis,” *Plants*, vol. 12, Sep. 2023, doi: [10.3390/plants12183257](https://doi.org/10.3390/plants12183257).

\[Pad23\] C. Padmanabhan *et al.*, “High-throughput sequencing application in the detection and discovery of viruses associated with the regulated citrus leprosis disease complex,” *Frontiers in Plant Science*, vol. 13, Jan. 2023, doi: [10.3389/fpls.2022.1058847](https://doi.org/10.3389/fpls.2022.1058847).

\[Yu26b\] L. Yu *et al.*, “The mRNA covalent modification dihydrouridine regulates transcript turnover and photosynthetic capacity during plant abiotic stress,” *The Plant Cell*, vol. 38, Jun. 2026, doi: [10.1093/plcell/koag196](https://doi.org/10.1093/plcell/koag196).

\[Dik25c\] D. Diksha *et al.*, “Peach viromics unveils it to be new host of grapevine red globe virus, citrus sudden death-associated virus, grapevine asteroid mosaic-associated virus and a novel marafivirus,” *Journal of Plant Pathology*, vol. 107, pp. 2041–2054, Jul. 2025, doi: [10.1007/s42161-025-01960-9](https://doi.org/10.1007/s42161-025-01960-9).

\[Nak22\] E. Nakasu, G. Silva, S. Montes, and A. F. S. Mello, “Virome analysis of sweetpotato in three Brazilian regions using high-throughput sequencing,” Oct. 11, 2022. doi: [10.1007/s40858-022-00532-x](https://doi.org/10.1007/s40858-022-00532-x).

\[Nis24\] H. Nishio *et al.*, “Circadian and environmental signal integration in a natural population of Arabidopsis,” *Proceedings of the National Academy of Sciences of the United States of America*, vol. 121, Jan. 2024, doi: [10.1073/pnas.2402697121](https://doi.org/10.1073/pnas.2402697121).

\[Col26\] C. Colvin and S. Chopra, “Population-scale transcriptomics reveals host genetic control of phyllosphere fungal communities,” May 26, 2026. doi: [10.64898/2026.05.21.727028](https://doi.org/10.64898/2026.05.21.727028).

\[Wei25c\] E. Weinheimer, S. Cory, N. Kortessis, T. M. Anderson, and J. B. Pease, “Differential gene reactions reveal drought response strategies in African acacias,” *The Plant Journal*, vol. 123, Aug. 2025, doi: [10.1111/tpj.70385](https://doi.org/10.1111/tpj.70385).

\[Paa26\] P. Paajanen *et al.*, “Circadian entrainment and gating in a natural plant population,” Aug. 30, 2026. doi: [10.64898/2026.01.23.701304](https://doi.org/10.64898/2026.01.23.701304).

\[Mwa20\] F. Mwatuni, A. Nyende, J. Njuguna, X. Zhonguo, E. M. Machuka, and F. Stomeo, “Occurrence, genetic diversity, and recombination of maize lethal necrosis disease-causing viruses in Kenya,” *Virus Research*, vol. 286, Jul. 2020, doi: [10.1016/j.virusres.2020.198081](https://doi.org/10.1016/j.virusres.2020.198081).

\[Kha23c\] Z. A. Khan *et al.*, “Genome analysis of viruses of Phenuiviridae, Betaflexiviridae and Bromoviridae, and apple scar skin viroid in pear by high-throughput sequencing revealing host expansion of a rubodvirus and an ilarvirus,” *Physiological and Molecular Plant Pathology*, Nov. 2023, doi: [10.1016/j.pmpp.2023.102196](https://doi.org/10.1016/j.pmpp.2023.102196).

\[Dan21\] L. Dantas *et al.*, “Field microenvironments regulate crop diel transcript and metabolite rhythms,” *bioRxiv*, Apr. 2021, doi: [10.1101/2021.04.08.439063](https://doi.org/10.1101/2021.04.08.439063).

\[Alb20\] T. Albrecht, S. White, M. L. Layton, M. Stenglein, S. Haley, and P. Nachappa, “Ecology and Epidemiology of Wheat Curl Mite and Mite-Transmissible Viruses in Colorado and Insights into the Wheat Virome,” Aug. 10, 2020. doi: [10.1101/2020.08.10.244806](https://doi.org/10.1101/2020.08.10.244806).

\[Cho18b\] C. K. J. D. Ibaba and A. Gubba, “The Conceivable Influence of Persistent Genotype-4 Hepatitis C Virus Infection on Cellular Immune Subsets before, during and after Pegylated Interferon-α and Ribavirin Therapy,” Jul. 12, 2018. doi: [10.4172/2161-0517-C2-027](https://doi.org/10.4172/2161-0517-C2-027).

\[Mat16\] E. Matsumura, “Estudo das populações de vírus presentes em plantas de citros cultivadas em uma região afetada pela morte súbita dos citros,” Dec. 12, 2016.

\[Zyk25\] P. A. Zykin, E. A. Andreeva, N. Tsvetkova, A. N. Bulanov, and A. V. Voylokov, “Rna-seq contamination as a metatranscriptomic data for screening of plant pests and symbionts,” *Ecological genetics*, Sep. 2025, doi: [10.17816/ecogen642484](https://doi.org/10.17816/ecogen642484).

\[Yu25e\] L. Yu *et al.*, “The mRNA covalent modification dihydrouridine regulates transcript turnover and photosynthetic capacity during plant abiotic stress,” *bioRxiv*, vol. 38, Nov. 2025, doi: [10.1101/2025.01.17.633510](https://doi.org/10.1101/2025.01.17.633510).

\[Hod17\] B. Hodge, P. Paul, and L. Stewart, “First report of Cocksfoot mottle virus infecting wheat (Triticum aestivum) in Ohio,” *Plant Disease*, vol. 102, pp. 464–464, Dec. 2017, doi: [10.1094/PDIS-08-17-1224-PDN](https://doi.org/10.1094/PDIS-08-17-1224-PDN).

\[Kaw20\] T. Kawakatsu *et al.*, “The transcriptomic landscapes of diverse rice cultivars grown under mild drought conditions,” Dec. 11, 2020. doi: [10.1101/2020.12.11.421685](https://doi.org/10.1101/2020.12.11.421685).

\[Rie19\] M. Rienth, S. Ghaffari, and Jean, “IMPACT OF GRAPEVINE LEAFROLL VIRUS INFECTIONS ON VINE PHYSIOLOGY AND THE BERRY TRANSCRIPTOME,” 2019.

\[Li25f\] C. Li, G. Zhang, G. Cheng, and Q. Wang, “Phenylalanine Ammonia-Lyase GhPAL9 Confers Resistance to Verticillium Wilt in Cotton,” *International Journal of Molecular Sciences*, vol. 26, May 2025, doi: [10.3390/ijms26114983](https://doi.org/10.3390/ijms26114983).

\[Lan18\] D. Lanver *et al.*, “The Biotrophic Development of Ustilago maydis Studied by RNA-Seq Analysis\[OPEN\],” *Plant Cell*, vol. 30, pp. 300–323, Jan. 2018, doi: [10.1105/tpc.17.00764](https://doi.org/10.1105/tpc.17.00764).

\[Kol21\] M. C. Kolodziej *et al.*, “A membrane-bound ankyrin repeat protein confers race-specific leaf rust disease resistance in wheat,” *Nature Communications*, vol. 12, Feb. 2021, doi: [10.1038/s41467-020-20777-x](https://doi.org/10.1038/s41467-020-20777-x).

\[Ber22\] D. Berry *et al.*, “Cross-species transcriptomics identifies core regulatory changes differentiating the asymptomatic asexual and virulent sexual life cycles of grass-symbiotic Epichloë fungi,” *G3: Genes\|Genomes\|Genetics*, vol. 12, Feb. 2022, doi: [10.1093/g3journal/jkac043](https://doi.org/10.1093/g3journal/jkac043).

\[Min20b\] J. Minicka, A. Zarzyńska‐Nowak, D. Budzyńska, N. Borodynko-Filas, and B. Hasiów‐Jaroszewska, “High-Throughput Sequencing Facilitates Discovery of New Plant Viruses in Poland,” *Plants*, vol. 9, Jun. 2020, doi: [10.3390/plants9070820](https://doi.org/10.3390/plants9070820).
