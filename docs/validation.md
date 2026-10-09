# What is known about the data

This document summarises what is established about the input data, cell
sets, metadata and shared methods of the PlaViSca atlas, and what remains
open. The evidence (download records, comparison tables, reports) is kept in
the archive `archive/pre_process_data-audit.tar` of the PlaViSca project
folder; paths below of the form `audit/...` refer to files in that archive.

## How the pipeline outputs are verified

- **Inputs.** Every input file has a checksum in `metadata/input_files.tsv`,
  `metadata/starsolo_counts.sha256` and `metadata/study_count_objects.sha256`
  (`./run.sh check-inputs`).
- **Checks while running.** The stages stop if a cell cannot be matched to
  its run or to the authors' data, if a cell gets a stage impossible for its
  tissue, or if cell numbers differ from `metadata/expected_counts.tsv`.
- **Count import.** The counts of the four STARsolo studies are produced in
  the current run from per-study count objects (the alternative entry point,
  `R/tools/import_count_objects.R`). For Sa2020 the full import from the
  STARsolo matrices (`R/10_import_counts.R`) was also run: it gives exactly
  the same cells, genes and counts. For Hazzard2022, Hazzard2024 and
  Ruberto2022_2 the import from the STARsolo matrices has not been run with
  the current code. For Hazzard2024, the STARsolo outputs of 4 of its 24 runs
  (SRR27021991 to SRR27021994) are not in the local inputs, so this import
  cannot be run until they are added.
- **Author matrices.** The import of Ruberto2022_1 and Mancio-Silva2022 from
  the authors' objects gives exactly the same cells, genes and counts as the
  per-study objects of data release v1.0-rc2.
- **Release comparison.** The outputs were compared, cell by cell, with
  data release v1.0-rc2:
  - identical: cells and their order, genes, counts, normalised
    expression, all metadata and stages, variable genes, the clusters at
    every resolution, the integrated UMAP, and all integration metrics;
  - 21 of the 50 principal components (and of the Harmony components) have
    the opposite sign. The sign of a component is arbitrary and changes no
    distance;
  - the unintegrated UMAP and the t-SNE differ slightly, as expected for
    approximate methods. UMAP: median shift 0.07 on a range of 21.5,
    correlation 0.9996 per axis, 43 cells shifted by more than 1 unit.
    t-SNE: median shift 0.0002, correlation 1.000.

The current atlas contains 105,949 cells from six studies:

| Study | Cells | Cell set taken from |
|---|---:|---|
| Mancio-Silva2022 (liver) | 1,480 | authors' processed Seurat object, minus 14 duplicate records |
| Sa2020 (blood) | 9,766 | STARsolo filtered cells (own reprocessing of ENA reads) |
| Hazzard2022 (sporozoites + blood) | 3,294 | own emptyDrops call on STARsolo raw matrices |
| Hazzard2024 (blood) | 80,024 | authors' deposited processed-cell list |
| Ruberto2022_1 (liver) | 1,438 | authors' processed Seurat object (cells and counts) |
| Ruberto2022_2 (sporozoites) | 9,947 | locally held combined object, `Stage == "Sporozoite"` |


## Mancio-Silva2022 (PMID 35443155)

**Inputs.**
- Authors' Zenodo record 6280956 (version 2.0), linked to GEO GSE197409.
- `PvData_CHM_Final.RDS`, MD5 `b58ba79595bee10121dc374aa160bc30`.
- `TableS2_Cell_Metadata.xls`, MD5 `5450e12b6f2a51c0e5197c4f30fd91c8`.
- Both MD5 values match the checksums advertised on Zenodo.
- The object has 1,494 cells, 78 metadata fields and 5,934 RNA features. Its counts are non-negative integers.
- Counts are taken directly from this object. PlaViSca does not realign this study.

**Accessions.**
- GEO GSE197409 is linked to PMID 35443155 and to PRJNA810174.
- ENA returns 58 runs, SRR18134227–SRR18134284 (SRP361380).
- The 15 retained source groups map to 15 GEO samples and 23 runs. Each group matches a unique GEO sample, and GEO day and treatment agree for every cell.
- Infection2 groups have two runs each. A cell is therefore linked to a set of runs, not to a single run.

**Cell set and duplicates.**
- Each source cell maps to exactly one record in the PlaViSca files. Raw counts are identical to the source for all 8,865,396 values.
- Normalised values equal `log1p(count / cell_total * 10000)` exactly.
- Table S2 has 1,494 rows but only 1,480 distinct identifiers.
- 14 records carry a terminal `.1`. For each one, both the `.1` record and its partner match the same single GEO count column across all 5,934 genes and all assay layers. They are duplicate representations of one upstream cell. These 14 `.1` records are removed. Their partners are kept.
- The mechanism that created these duplicates on the author side is unknown.
- `D5Seq1_CCCCGATTGACG` and `D5Seq2_CCCCGATTGACG` come from different GEO samples (SEQ5r1, SEQ5r2) and different runs. Their count profiles are still identical over all 85,578 host and parasite features in the GEO combined matrix (612 counts, 581 of them parasite).
- Both cells are kept as ordinary cells, with no flag. This is a policy decision: no evidence proves duplication, and none proves independence.
- No other exact duplicate profiles exist among the 1,494 source records.

**Metadata verified.**
- Days d1/d4/d5/d8/d11/d14 (207/80/306/65/603/233 cells, before duplicate removal) come from the source `Day` field and agree with GEO.
- Treatment: source CTRL (1,341) / PI4K (153). The article identifies PI4K as MMV390048 at 1 µM, days 5–8. CTRL is kept as "study control", not "no treatment".
- Technology: Seq-Well S3. Host: primary human hepatocytes procured in the USA.
- Parasite origin: the publication supports Thailand at country level only. It names three provinces (Tak, Songkla, Ubon-Ratchathani) but does not say which infection came from which province.

**Discrepancies.**
- Sequencer. GEO/ENA register Infection1 as NextSeq 500 and Infection2 as NovaSeq 6000. The article instead says NovaSeq for uncaptured material and NextSeq 500 for captured material in both infections. The Infection2 ENA experiment record also describes a NextSeq capture protocol under a NovaSeq label. Only the registered instrument is stored, and the field for the instrument stated in the publication is empty.
- Stage labels. Source `State` is Hypnozoite or Replicative. In the current atlas, all 1,208 retained Replicative cells are labelled "Schizont (Liver stage)" (272 are Hypnozoite). The study audit judged that converting every Replicative cell to schizont is not justified by `State` alone. No document records that this was resolved.
- Origin location. The current data give the parasite origin as `Thailand_Ubon-Ratchathani` for every cell. This is more specific than the evidence, which supports Thailand only.

## Sa2020 (PMID 32365102)

**Inputs.**
- ENA PRJNA603327 contains exactly the ten runs used, SRR11008269–SRR11008278. All are registered as Illumina HiSeq 4000.
- Counts are PlaViSca's own STARsolo `GeneFull/filtered` output: all filtered barcodes, with no author list and no further QC.
- Publisher supplements S1, S5 and S7 Data, the article XML and the ENA snapshot are kept with SHA-256 values.
- The local S5 Data workbook is byte-identical to the publisher download (MD5 `03ad55660ba413d7bd8b6140e98c6397`).
- No FASTQ checksums are recorded for this study.

**Cell set compared with the authors.**
- The authors' main analysis has 9,215 cells (5,000–75,000 unique reads after PCR deduplication, HISAT2 pipeline). Their relaxed analysis reports 13,503 cells (1,000-read cutoff), but S7 Data contains 13,502 unique IDs.
- Author prefixes were matched to runs using Table 1 cell totals, ENA aliases and barcode overlap. The labels are historical and misleading: for example, `NIH_Ao` is the Indonesia-I Aotus library and `PB_MACS` is the second NIH-1993 Saimiri library.
- The 9,766 PlaViSca cells split into three groups:
  - 9,018 cells in the authors' main analysis;
  - 414 cells only in the relaxed list;
  - 334 cells found only by the reprocessing.
- 197 main-analysis cells are not in PlaViSca.
- Per run, the shortfall is 0 to 57 cells, and the excess is 9 to 233 (AMRU_Sa: 22 author cells vs 255 PlaViSca cells).
- All 9,766 cells are kept. A per-cell flag (`sa2020_main_analysis`) marks the 9,018 cells in the main analysis.
- The extra cells are not shown to be invalid droplets. The two methods (STARsolo UMIs and author read windows) count differently. Why individual cells were gained or lost is not established.

**Expression.**
- The counts in the study object equal the STARsolo filtered matrices exactly, with zero unequal entries in all ten runs.
- The 50 rRNA genes are removed, and 1,167 genes that are zero in every cell are dropped, leaving 5,644 genes.
- No profile is an exact duplicate. 183 barcode tokens recur across runs, which is expected; cell IDs include the run.
- There is no comparison with the authors' own counts, because no author count matrix was deposited.

**Metadata verified against publication Table 1 and ENA.**
- The following match for all 9,766 cells: host species (Saimiri boliviensis / Aotus nancymaae), animal ID, strain (NIH-1993, Indonesia-I, Chesson, AMRU-I), collection day and chloroquine status.
- The three chloroquine samples were collected 16 h after one oral dose: 5 mg/kg for Indonesia-I, and 10 mg/kg for Chesson and AMRU-I. They come from the same Aotus animals as the matching untreated libraries.
- 10x Chromium 3′ is supported, but the chemistry version is not. The data record no version.
- "USA_Rockville" is the registry location, not where the parasites come from. The origin of NIH-1993 is recorded as unresolved: it is closely related to Salvador-I but not the same isolate.

**Author sex labels.**
- S5 Data has 9,158 IDs (7,587 Asexual, 1,477 Female, 94 Male). Of the 1,571 sex-labelled rows, 1,532 are present in PlaViSca.
- These 1,532 labels are imported with an exact key and take precedence over the PlaViSca cluster call.
- In 124 of these cells the cluster call disagrees: 20 have the opposite sex and 104 have no gametocyte call.

## Hazzard2022 (PMID 36525464)

**Inputs.**
- ENA PRJNA863611 has 8 runs: 4 PacBio (not used) and 4 Illumina, SRR20710498–501.
- Sample aliases: 498 = illumina_spz_sample_2 and 499 = illumina_spz_sample_1 (sporozoites from An. freeborni and An. stephensi); 500/501 = Saimiri boliviensis blood.
- Counts are PlaViSca's own STARsolo `GeneFull/raw` matrices. Cells are called with seeded emptyDrops (lower = 1, FDR ≤ 0.001).

**Cell set.**
- Retained cells per run: 1,433 / 1,425 / 362 / 74 (2,858 sporozoites, 436 blood cells).
- Every retained barcode exists in the raw STARsolo matrix of its run.
- The cell lineage is exact, and counts in the study object equal the exported raw values.
- Normalised values follow the stated formula exactly (17,517,492 values, zero mismatches).
- No exact duplicate profiles were found.
- Count values were not compared element by element with the raw STARsolo matrix; only barcode presence was checked.
- Sporozoite cells have low depth. Median UMI is 103 (run 498) and 114 (run 499), against 1,545 and 1,688 for the blood runs.

**Not reconciled with the publication.**
- The publication uses read/mapping QC and reports 2,609 + 2,363 sporozoite transcriptomes. PlaViSca keeps 2,858 sporozoites under its own QC.
- How the two sets correspond (cell IDs, count definitions) was never checked. The size of the blood population in the publication is not established.
- No complete audit of this study's publication, accessions and QC exists. Unlike the other five studies, there is no dedicated study-level report.

**Metadata.**
- Sequencer: ENA registers NovaSeq 6000 for 498/499 and HiSeq 2500 for 500/501. The publication text says NovaSeq. The cells carry the registered instrument; `metadata/runs.tsv` also has the published one (`sequencer_in_publication`). The conflict is not resolved.
- Strain: the data use Sal1/Chesson (mosquito runs) and Sal1 (blood), as in the registry. The publication text describes NIH/Chesson mosquito feeds. This is not reconciled.
- Day post infection and host ID are empty. The publication states sporozoites were collected 21 days after the feed, but this was not mapped to runs.
- The sporozoites are labelled by their source (all 2,858 are "Sporozoite"). The 48 gametocyte calls (46 female, 2 male) are all among the 436 blood cells.

## Hazzard2024 (PMID 39223117)

**Inputs.**
- Processed-cell list `Proccessed_Data.txt` from Zenodo 12775216, MD5 `6ede1f7da4e5d57841c9f90883b2b2c1`, which matches the advertised value.
- BioProject PRJNA1047651 (SRP475376). Author code: GitHub `bhazzard11/Coinfection_Analysis`, Zenodo 12775129.
- Counts are PlaViSca's own STARsolo processing, subset to the author list.
- The authors' count table (`Single_Cell_Count_Table.txt.zip`, advertised MD5 `276aa89b96636e9091a05cd43166465e`) could not be downloaded (HTTP 504). PlaViSca's counts have therefore never been compared with the authors' values.

**Runs.**
- ENA lists 43 runs, all NovaSeq 6000: 24 blood scRNA-seq, 6 salivary-gland sporozoite, 11 auxiliary and 2 WGS runs.
- The 24 blood runs contribute cells (`metadata/runs.tsv`). The sporozoite and auxiliary runs are not used.
- The publication describes 8 sporozoite libraries but ENA has 6 sporozoite runs. How they correspond is unresolved; this does not affect the atlas.

**Cell set.**
- The publication reports 80,024 cells from 24 libraries and 10 female Saimiri boliviensis. The table below gives the cells per design.

| Design | Animals | Libraries per animal | Cells |
|---|---:|---:|---:|
| Chesson mono-infection | 3 | 1–2 | 22,847 |
| NIH-1993-F3 mono-infection | 3 | 1–2 | 16,292 |
| Consecutive | 2 | 3 | 18,286 |
| Simultaneous | 2 | 2 | 11,715 |
| Sporozoite-initiated | 2 | 3 | 10,884 |

- The author list, the study object and the exports agree exactly (80,024 unique IDs).
- There are 6,583 genes and no exact duplicate profiles.
- Counts in the study object equal the exported raw and normalised values exactly. Comparison with raw STARsolo values was not performed.

**Metadata verified against the publication, ENA aliases and the author table.**
- Run 986 is library 5537_2 (3,489 cells); run 993 is library 5708_1 (16,321 cells, Chesson).
- Runs 970/971/972 were sampled on days 30/24/29 (9,644 cells).
- Animal is the replicate unit. Libraries from the same animal are repeated measures.
- Library, animal, lineage, day and design fields now match the reference run table for all 80,024 cells.
- NCBI BioSample records register three libraries as host "Homo sapiens" in "USA: Bethesda". This conflicts with the publication and is kept only as registry information.
- The origin of NIH-1993-F3 is recorded as unresolved. The data do not assign it Salvador-I's El Salvador origin.

**Author stages compared with PlaViSca stages.**
- The author table has PseudoGroup labels: GroupA 39,444, GroupB 12,422, GroupC 7,149, Sexual 21,009. These are kept as a separate field and are not used as stage labels.
- PlaViSca calls 20,769 female and 880 male gametocytes in this study.
- Of the 21,009 author Sexual cells, 20,392 receive a gametocyte call and 617 do not. 1,257 GroupA–C cells receive a gametocyte call.
- Hazzard2024 is about 75.5% of the atlas, so it dominates the shared embedding and clustering.

## Ruberto2022_1 (liver; PMID 36093191)

**Inputs.**
- BioProject PRJNA843856. Runs SRR19573609–612 (infected hepatocytes) are used; SRR19573613/614 (uninfected controls) are not.
- All runs are registered as HiSeq X Ten, in agreement with the publication.
- Counts come from the authors' processed object `Hep59.1.2.seu_20aug2025.rds` (MD5 `203013218936548744516954cb73c39e`), `RNA` counts layer, for all 1,438 cells.
- This object matches, by name, merge structure and cluster-renaming logic, the `Hep59.1.2.seu` built in the authors' GitHub workflow (`AnthonyRuberto/Pv_LS_singleCell`, additional file 3).
- It was never compared byte for byte with the deposited copy in Zenodo
  6463338. That file sits inside a 14.5 GB archive which was not downloaded.

**Why the author matrix is used instead of PlaViSca's own counts.**
- PlaViSca's STARsolo counts for the 1,438 author barcodes were near-empty for 538 cells: 517 with 0 UMI and 21 with 1 UMI.
- These cells were 89/90 cells of SRR19573609 and 449/457 of SRR19573610, which form biological replicate 2. No replicate-1 cell was affected.
- In the authors' object the same cells have real signal.
- The local FASTQ files match ENA read counts exactly.
- A controlled kallisto|bustools reprocessing, following the authors' own method, found plenty of parasite-mappable reads in runs 609/610: 2.6% of reads, against 0.6–2.0% in the control runs.
- That signal is spread over about 170,000–195,000 barcodes. At the 538 known barcodes, 98% have zero UMI and the maximum is 1.
- The same method recovers strong signal at the known cells of the control runs 611/612.
- Conclusion: the public replicate-2 data cannot reproduce the authors' cells. The problem lies upstream of alignment, in the library or the deposit, not in PlaViSca.
- Reasons to trust the author matrix:
  - the counts are non-negative integers;
  - the cell crosswalk (51inf→612, 52inf→610, 91inf→611, 92inf→609) is exactly one-to-one;
  - all 4,722 author genes match PlaViSca gene IDs exactly;
  - for the 891 unaffected cells of runs 611/612, each author cell's best-correlated PlaViSca cell is its true counterpart 100% (611) and 85% (612) of the time, against 0.39% and 0.16% by chance.
- Using the author matrix for all 1,438 cells avoids mixing two counting pipelines within one study.
- The adopted counts equal the author matrix exactly (6,790,236 entries).
- The author matrix comes from a different pipeline: kallisto 0.46.0 / bustools 0.39.3, a combined GRCh38 + PlasmoDB-51 index, and emptyDrops (lower = 1000). The other studies use STARsolo on PlasmoDB-68.
- 2,089 PlaViSca genes are absent from the author matrix and are filled with zero:
  - 288 ncRNA/pseudogenes excluded by the authors' protein-coding filter;
  - 1,801 protein-coding genes inferred to be missing from PlasmoDB-51. This was not checked gene by gene against the PlasmoDB-51 GFF.
- The cells' `count_source` field says that the counts are the authors' (kallisto|bustools, PlasmoDB-51).

**Metadata verified against ENA sample titles and the Methods.**

| Run | Day | Replicate | Treatment | Cells |
|---|---:|---:|---|---:|
| SRR19573609 | 9 | 2 | MMV390048 | 90 |
| SRR19573610 | 5 | 2 | none | 457 |
| SRR19573611 | 9 | 1 | MMV390048 | 255 |
| SRR19573612 | 5 | 1 | none | 636 |

- 1,438 cells matches the publication's stated total.
- Liver form comes from the source: 1,147 Hypnozoite and 291 Schizont. All 345 day-9 cells are hypnozoites, as the authors' workflow implies.
- Hepatocyte donor: BioIVT lot BGW. Parasites: Cambodian field isolates from Mondulkiri.
- 10x chemistry: the cells' `sc_technology` is "10x_Chromium" with no version. The authors' scripts use the 10x v3 whitelist and `kallisto bus -x 10xv3`, and the STARsolo runs of this study used the v3 whitelist; v3 is not stated in the publication text.

## Ruberto2022_2 (sporozoites; PMID 35926062)

**Inputs.**
- ENA PRJEB42435 has exactly three runs, ERR5087438–440, with samples SAMEA7804580–82 (rep1–3). All are HiSeq X Ten.
- Counts are PlaViSca's own STARsolo `GeneFull/raw` matrices (10x v3 whitelist, 6,794,880 barcodes per run).
- Cells are selected from the locally held object `6_PvSPZ.BS.combined.rds`.
- That object was integrated locally (Seurat 4.1.0, CCA integration, scDblFinder and scater fields). It holds 9,947 Sporozoite and 8,304 BloodStage cells.
- It is not a minimally modified author download. No build script for it exists, and no checksum against an author deposit is recorded.
- Whether its `Stage` label came from the authors or was added locally is unknown. Author deposits are on Zenodo 6474355.

**Cell set.**
- The cell set is defined entirely by the source list. No other QC is applied.
- Sample-to-run mapping: Case_909→438, Case_922→439, Case_923→440, with 3,375 / 4,002 / 2,570 cells.
- These counts are identical in the source, the STARsolo barcodes, the study object and the exports.
- 9,947 matches the publication. Lineage is exact, with no missing, extra or duplicate IDs.
- STARsolo counts = study object = exported raw values (67,749,017 values, zero differences). Normalised values follow the formula exactly.
- There are no exact duplicate profiles. Twelve barcode tokens recur across runs, which is harmless.
- PlaViSca counts were not compared with the source object's own counts, which cover a smaller gene set (4,170–4,254 features).

**Metadata verified.**
- Anopheles dirus salivary glands.
- Days 16/18/17 after the mosquito feed (ENA `dev_stage`).
- Three independent patient isolates from Mondulkiri (Kaev Seima), Cambodia.
- 10x Chromium v3, confirmed in the Methods.
- Sequencing at Macrogen, Seoul. Sample preparation at Institut Pasteur du Cambodge, Phnom Penh (the ENA coordinates).
- All cells are labelled Sporozoite because of how they were selected.

## Cross-study methods, annotation and integration

**Features and normalisation.**
- STARsolo studies use the PlasmoDB-68 P. vivax P01 annotation (6,861 genes). The 50 GFF rRNA genes are removed, leaving 6,811.
- Mancio-Silva2022 carried 49 rRNA genes from the author object (206,690 UMIs). These are now removed, so the atlas has 6,811 genes.
- Integration uses only the 4,709 genes measured in every study. From these, 2,000 variable genes are used for PCA and Harmony.
- Genes a study did not measure are recorded per study, separately from true zeros.
- Normalisation is `log1p(count / cell_total * 10000)` for all studies.

**Life-cycle stage.**
- The harmonised stage uses this order of precedence:
  1. a stage defined by the source (Ruberto2022_2 and Hazzard2022 sporozoites, Ruberto2022_1 liver form, Mancio-Silva2022 state, Sa2020 author sex labels);
  2. a cluster-level gametocyte call;
  3. a blood-stage bin from the reference.
- Author labels and PlaViSca predictions are kept in separate fields.
- Gametocyte calls are made per cluster, from marker scores (19 female and 6 male markers) at Louvain resolution 0.05.
- Marker discrimination: female AUC 0.862 (22,278 cells), male AUC 0.979 (988 cells). All calls lie in clusters that pass the marker threshold.
- The calls are stable: 23,137 of 23,266 are unchanged at resolution 0.08, and every resolution from 0.02 to 0.1 agrees with 0.05 on at least 99.5% of cells.
- In total, 23,370 cells carry a gametocyte stage.
- Blood-stage time labels come from SingleR against the Zhu et al. 2016 SMRU1 RNA-seq time course (4,908 genes with a PvP01 ortholog x 7 time points).
- They measure similarity to a reference time point. They are not measured hours after invasion. The stage of each time point (6-12 h ring, 24 h trophozoite, 32-42 h schizont, 48 h late schizont/merozoite; `metadata/idc_timepoints.tsv`) is a PlaViSca choice, not taken from the literature.
- These labels are computed for host-blood cells with at least 10 UMIs, including the Sa2020 cells labelled as gametocytes by the authors (whose final stage stays the authors' label).
- The pruned SingleR label and the score delta are stored for each cell.
- Checks on the stage labels:
  - no liver or mosquito cell has a blood-stage label or a gametocyte call;
  - no liver cell is a gametocyte;
  - seven Hazzard2024 cells with 5–9 UMI have no stage at all;
  - no impossible combinations of stage and source were found.

**Integration.**
- Harmony runs per study. The six input layers are checked to map one-to-one to the six studies.
- Seeds: PCA and UMAP 42, t-SNE 123. Harmony uses the random number state left by the unintegrated UMAP, so the order of the steps in `R/30_integrate.R` is part of the method.
- Mixing was assessed only within populations that are biologically comparable. The table gives the mean fraction of each cell's 30 nearest neighbours that come from its own study.

| Broad stage | Cells | Studies | Before | After |
|---|---:|---:|---:|---:|
| Sporozoite | 12,805 | 2 | 0.864 | 0.661 |
| Blood | 90,219 | 3 | 0.996 | 0.919 |
| Liver | 2,918 | 2 | 0.972 | 0.830 |

- Every study keeps neighbourhoods enriched for itself, well above what random mixing would give.
- Whether one global integration is the right design was tested only through these metrics. Study and biology are strongly confounded.

## Open points

1. **Ruberto2022_1 author object, byte identity.** The local `Hep59.1.2.seu_20aug2025.rds` was never checksum-compared with the Zenodo 6463338 deposit. The authors were not asked whether it was built from the public FASTQ files. See `audit/ruberto2022_1/author_matrix_validation/REPORT.md` (section 12).
2. **Ruberto2022_1 replicate 2.** Why the public FASTQ files for SRR19573609/610 do not carry signal at the real cell barcodes is unknown. Possible causes are ambient contamination or barcode corruption at deposit. See `audit/ruberto2022_1/D036_resolution.md`.
3. **Ruberto2022_1 gene gap.** That the 1,801 protein-coding genes are absent from PlasmoDB-51 is inferred, not checked against that GFF. See `audit/ruberto2022_1/phase2_feature_coverage.tsv`.
4. **Hazzard2022 QC.** How PlaViSca's 2,858 sporozoites relate to the publication's 2,609 + 2,363 was never checked. The publication's blood population is unknown. See `audit/mancio_silva2022/AUDIT.md`, section "Hazzard2022: confirmed stage contradiction, different QC".
5. **Hazzard2022 sequencer and strain.** Sequencer: HiSeq 2500 in ENA vs NovaSeq in the publication for runs 500/501. Strain: Sal1/Chesson in the registry vs NIH/Chesson in the text. Neither is resolved. See `audit/mancio_silva2022/cross_study_ena_hazzard2022.tsv`.
6. **Hazzard2022 day and host.** These are empty. The 21-day collection time is not mapped to runs. See the same section of `audit/mancio_silva2022/AUDIT.md`.
7. **Hazzard2022 count values.** STARsolo raw counts were not compared element by element with the study object; only barcode presence was checked. See `audit/mancio_silva2022/hazzard2022_raw_inventory.tsv`.
8. **Hazzard2024 author counts.** The deposited count table could not be downloaded, so PlaViSca counts are not verified against it. See `audit/hazzard2024/AUDIT.md`.
9. **Hazzard2024 sporozoite runs and NIH-1993-F3.** How 8 published sporozoite libraries relate to 6 ENA runs is unresolved. So is the geographic origin of NIH-1993-F3. See `audit/hazzard2024/all_deposited_run_membership.tsv`.
10. **Sa2020 population.** The 748 cells outside the main analysis and the 197 missing author cells are not explained cell by cell. S7 Data has 13,502 IDs against 13,503 in the article. See `audit/mancio_silva2022/sa2020_inventory_summary.tsv`.
11. **Sa2020 and Hazzard2024 author counts.** No author raw count matrix was compared for either study. See `audit/mancio_silva2022/sa2020_feature_disposition.tsv`.
12. **Mancio-Silva2022 sequencer.** The registered and published instruments for Infection2 disagree, and which run is pre- or post-capture is unknown. See `audit/mancio_silva2022/source_group_accession_map.tsv`.
13. **Mancio-Silva2022 province.** The data assert Ubon-Ratchathani although only Thailand is supported. See `audit/phase1_code_review/defect_fix_verification.tsv`.
14. **Mancio-Silva2022 stages.** All 1,208 Replicative cells are labelled "Schizont (Liver stage)". The study audit considered this unsupported. See `audit/phase1_code_review/REVIEW.md`.
15. **Mancio-Silva2022 duplicates.** The author-side mechanism behind the 14 removed duplicates is unknown. Why the two retained day-5 cells are identical is unexplained. See `audit/mancio_silva2022/exact_duplicate_trace.tsv` and `audit/mancio_silva2022/DUPLICATE_TABLE.md`.
16. **Ruberto2022_2 source object.** `6_PvSPZ.BS.combined.rds` has no build script and no checksum against Zenodo 6474355. The origin of its `Stage` label is unknown. See `audit/ruberto2022_2/source_object_provenance_summary.txt`.
17. **Integration design.** Global Harmony with strongly confounded studies is supported only by neighbourhood metrics. No stratified alternative was tested. See `audit/phase3_preflight/batch_mixing_within_comparable_populations.tsv`.
18. **Stage biology.** Gametocyte calls and blood-stage stages are checked against markers and robustness tests, not against independent biological ground truth. The SingleR uncertainty tables in `audit/phase3_preflight/` were computed on an atlas with the same cells but before the current gametocyte precedence; they were not recomputed.
19. **STARsolo import.** Not run with the current code for Hazzard2022, Hazzard2024 and Ruberto2022_2; four Hazzard2024 runs are missing from the local STARsolo outputs (see the top of this document).
20. **Input sources.** The source of several input files is not documented (`to be confirmed` in `metadata/input_files.tsv`): the Zhu et al. 2016 table, the ortholog table, Mancio-Silva2022 Table S2 and the two Ruberto objects.
