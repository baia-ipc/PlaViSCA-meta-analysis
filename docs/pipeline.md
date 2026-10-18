# The pipeline, stage by stage

All parameters are in `config/params.yaml`. All curated tables are in
`metadata/`. The numbers obtained (cells, genes, metrics) are written by the
stages to the run's `results/` and are not repeated here.

## 00 Reference

- Reads the PlasmoDB-68 PvivaxP01 GFF. Every `protein_coding_gene`,
  `ncRNA_gene` and `pseudogene` record is a gene.
- rRNA genes are the `ncRNA_gene` records with `ebi_biotype=rRNA`.
  Protein-coding genes whose description mentions rRNA (e.g.
  methyltransferases) are not rRNA.
- Reads the gene list of every STARsolo run (`features.tsv`) and checks that
  all runs were counted against the same genes.

Outputs: `genes.tsv` (gene ID, name, description, type, rRNA flag),
`starsolo_genes.tsv`, `gene_annotation.rds` (for the app).

## 10 Import counts

Each study is imported as one Seurat object with the raw UMI counts. Cells
are named `<last 3 digits of the run accession>_<barcode>`. Mancio-Silva2022
keeps the authors' cell names. `author_cell_id` is the cell's name in the
authors' data.

STARsolo studies (GeneFull counts; settings under `import` in params):

| Study | Matrix | Cells |
|---|---|---|
| Sa2020 | raw | the cells of the authors' main analysis (`metadata/sa2020_main_analysis_cells.tsv`) |
| Hazzard2022 | raw | the authors' rule: at least 5,000 (blood) or 250 (salivary gland) unique reads over the whole genome |
| Hazzard2024 | raw | the cells in the authors' processed table (Zenodo 12775216) |
| Ruberto2022_2 | raw | the cells with Stage = Sporozoite in the authors' object |

For these studies:
- rRNA genes are removed before the objects are made.
- Each run is made into an object with `min.cells`/`min.features` as in
  params, and the runs are merged.
- A barcode-rank (knee) plot is drawn for each run.
- The runs of a study are read in parallel, on `PLAVISCA_THREADS` cores
  (set in `config/local.env`).
- Hazzard2022: the authors define a cell by its unique reads over the whole
  genome (reads with the same barcode, UMI, strand and position counted
  once), not only on genes. `starsolo/run_starsolo.sh --genome-reads` writes
  these counts per barcode (`genome_reads_per_barcode.tsv` in each run's
  folder); the thresholds are `import.Hazzard2022.min_genome_reads`. The
  deposited blood runs give fewer cells than the publication reports (879
  against 2,756), also when processed with the authors' published method.
- To match a cell with the authors' data, the authors' sample names are
  translated to runs using `metadata/runs.tsv` (`author_sample`).

`metadata/sa2020_main_analysis_cells.tsv` lists the authors' names of the
9,215 cells of their main analysis. The authors did not deposit a cell list;
the names come from the source data of their Figure 1B (S1 Data,
`author_data/pbio.3000711.s030.xlsx`, sheet `data 1B`), which gives every
cell of the figure by sample name and barcode. The P. berghei cells of the
same sheet (`Pberg_`) are left out; the remaining 9,215 match the parasite
numbers per sample of the authors' Table 1. Sample names are matched to runs
with `metadata/runs.tsv` (`author_sample`). Stage 10 keeps exactly these
cells (all are in the STARsolo raw matrices) and gives each its
`author_cell_id`.

Studies whose counts come from the authors:

- **Ruberto2022_1**: the raw UMI counts (`counts` layer) of the authors'
  object. The authors counted with kallisto|bustools against PlasmoDB-51 and
  kept protein-coding genes only. Their matrix is placed on our gene list:
  genes they did not count are zero rows.
  `metadata/ruberto2022_1_genes.tsv` lists, for every gene, whether the
  authors counted it and, if not, why. The authors' liver forms
  (`LiverForm`) are kept for stage 20.
- **Mancio-Silva2022**: the raw UMI counts of the authors' Zenodo object
  (Seq-Well). Gene names are cut to the PVP01 ID. 14 cells appear twice: a
  second copy with a `.1` suffix and identical counts. The copies are
  removed (`metadata/mancio_silva2022_duplicate_cells.tsv`), and the stage
  checks that their counts are identical.

No further filtering on UMIs or genes per cell, and no doublet detection, is
applied. Each study keeps the cells defined by its source.

## 20 Cell metadata

Every cell gets:
- its study (`metadata/studies.tsv`);
- the run-level information of its run (`metadata/runs.tsv`): host, tissue,
  strain, origin, infection design, day post infection, treatment,
  replicate, technology, sequencer;
- the authors' own annotations:

| Study | Annotation | Source |
|---|---|---|
| Sa2020 | female/male gametocytes → `author_stage`, `author_sex` | S5 Data |
| Ruberto2022_1 | Hypnozoites / Schizonts → `author_stage` | authors' object |
| Mancio-Silva2022 | State, sex.state, AP2G; State gives `author_stage` (Replicative = liver schizont) | Table S2 |
| Hazzard2024 | sexual/asexual group → `author_pseudogroup` | processed table |
| Ruberto2022_2, Hazzard2022 (runs SRR20710498/9) | the whole sample is sporozoites → `author_stage` | `metadata/runs.tsv`, `sample_stage` |

## 30 Integrate

1. **rRNA removal.** The rRNA genes are removed from every study. Only the
   Mancio-Silva2022 matrix still contains them at this point. The UMIs and
   genes per cell are then recomputed.
2. **Measured genes.** For each study, the genes it was counted against:
   - STARsolo studies: the STARsolo gene list;
   - Ruberto2022_1: the genes in the authors' matrix;
   - Mancio-Silva2022: the genes of its matrix.

   Detection is not required, so stage-specific genes stay measured.
3. **Merge and normalise.** The studies are merged in the order
   `integration$merge_order`, with one count layer per study, and
   log-normalised (`NormalizeData`, scale factor 10,000).
4. **Variable genes and PCA.** The variable genes (2,000, `vst`) are chosen
   among the genes measured in all six studies. A gene missing from one
   study's annotation therefore cannot look variable. Then come `ScaleData`
   on these genes and `RunPCA` (50 components, seed 42).
5. **Embeddings and clusters.**
   - Unintegrated: UMAP on PCs 1–30 (3 components, seed 42).
   - Harmony (`IntegrateLayers`, default settings): with one layer per study,
     the study is the batch variable.
   - On the Harmony PCs 1–30: UMAP (3 components); a shared nearest-neighbour
     graph (k = 20); Louvain clusters at resolutions 0.01 to 0.10; t-SNE (3
     dimensions, seed 123).

Harmony draws random numbers from the state left by the unintegrated UMAP.
The order of the steps in the script is therefore part of the method.

Outputs: `atlas_integrated.rds`, `rrna_genes_removed.tsv`,
`measured_genes.tsv` (which genes each study measured).

## 35 Blood embedding

The cells eligible for blood-stage inference (host blood, at least 20
UMIs) are analysed on their own, so that the asexual cycle and the sexual
branch shape their own embedding, without the liver and sporozoite cells:
- 2,000 variable genes, PCA, and Harmony with one batch per study;
- a 2D UMAP on the Harmony PCs 1–15, with 100 neighbours and minimum
  distance 0.5 (`umap_b_1`, `umap_b_2` in the exported data);
- Louvain clusters at resolution 1 on a 50-neighbour graph of the Harmony
  PCs 1–15 (`blood_cluster`), used by stage 40 to call gametocytes.

## 40 Annotate stages

1. **Eligible cells.** Cells of host blood samples with at least 20 UMIs
   are eligible for blood-stage inference.
2. **Asexual stage.** SingleR (default settings) compares each eligible
   cell with the SMRU1 time course of Zhu et al. 2016:
   - the reference is RNA-seq with 7 time points (6, 12, 24, 32, 36, 42 and
     48 h);
   - its genes are mapped to PvP01 by PlasmoDB's one-to-one orthologs, and
     genes without an ortholog are dropped.

   `idc_timepoint` is SingleR's label: the time point the cell most
   resembles. It is a similarity, not a measured age. Every eligible cell
   keeps its label; `idc_pruned` is TRUE where SingleR's pruning (default
   settings) would discard it as unreliable. `idc_delta` is SingleR's
   `delta.next`. `metadata/idc_timepoints.tsv` turns time points
   into stages (`idc_stage`): 6–12 h ring/early trophozoite, 24 h trophozoite, 32–42 h
   schizont, 48 h late schizont/merozoite.
3. **Gametocytes.** Gametocytes are called from the
   marker sets of `metadata/gametocyte_markers.tsv`: 19 female and 4 male
   genes from Sa et al. 2020, and a panel of 8 known male gametocyte
   markers (Lasonder et al. 2016, Table 1: kinesin-8B, three dynein heavy
   chains, CDPK4, MAPK2, alpha-tubulin 2, HAP2). The S4 dot plot also shows
   5 known female markers from the same table (P25, P47, NEK2, NEK4, DOZI).
   The P. falciparum ortholog of each panel gene is given in the table.
   The Sa et al. 2020 sets are the genes highlighted as specifically
   expressed in gametocytes in their S3 Table. Two of its six male genes are
   left out: fructose-bisphosphate aldolase (PVP01_1262200), also detected
   in most asexual parasites, and allantoicase (PVP01_1266500), also detected
   in most late female gametocytes.
   For each blood-stage cluster of stage 35 and each marker set, the
   enrichment is the share of the cluster's cells detecting the set's genes
   minus the share detecting control genes of similar expression (100 per
   marker gene), in percentage points. Each cluster's median UMIs are
   divided by the median over the clusters (relative UMIs). Neither depends
   on sequencing depth. Each cluster is assigned by this rule (thresholds in
   `config/params.yaml`):

   | Assignment | Condition (per cluster) |
   |---|---|
   | female branch | female enrichment ≥ 10 and > male enrichment |
   | male branch | male enrichment ≥ 10 and > female enrichment |
   | late male | male branch with mature-male panel enrichment ≥ 20 |
   | early male | other male-branch clusters |
   | developing gametocyte | female branch with relative UMIs < 0.4 (sex not called) |
   | late female | female branch with relative UMIs ≥ 1.5 |
   | early female | other female-branch clusters |

   Module scores of the female and male sets (`AddModuleScore`) are kept per
   cell for display (`female_gam_score`, `male_gam_score`).

   The UMI cut-offs use the growth of gametocytes as they mature (more RNA
   per cell). Developing gametocytes get the stage "Developing gametocyte",
   not a sex. Eligible cells take the call (`cluster_gametocyte_call`) and
   the state (`gametocyte_state`) of their cluster.
4. **Final stage.**
   - `stage` of eligible (host-blood) cells: the cluster gametocyte call,
     else `idc_stage`. All other cells (sporozoites, liver stages) take the
     authors' stage, as these are cells we do not classify ourselves.
     `stage_source` records which one was used.
   - `broad_stage` (blood, liver or sporozoite) comes from
     `metadata/stages.tsv`.
   - The authors' gametocyte labels of Sa2020 (`author_sex`) and the
     sexual/asexual groups of Hazzard2024 (`author_pseudogroup`) do not
     enter `stage` of host-blood cells; they are compared with the cluster
     call in `author_labels_vs_cluster_call.tsv`.

Checks:
- only eligible cells have blood-stage labels;
- gametocytes occur only in host blood;
- liver cells have liver stages.

Results (`results/stage_annotation/`):
- `eligible_cells.tsv`;
- `idc_reference.tsv` (genes of the reference and how many were mapped);
- `cluster_gametocyte_scores.tsv` (mean scores, median UMIs, call and
  state per blood-stage cluster);
- `gametocyte_marker_panel.tsv` (fraction of cells detecting each curated
  sexual-stage gene, and its mean expression, per gametocyte state);
- `author_labels_vs_cluster_call.tsv`.

## 50 Export of the app data

`results/app_data/` contains:
- `cleaned_dataset.rds`, a list of:
  - `mr_data`: the per-cell columns of `metadata/columns.tsv`, in that order,
    cells sorted by name;
  - `pca_df`: per-study PCA standard deviations;
  - `hiv_data`: per-study variable genes, recomputed on each study alone;
- `normalize_df.rds`: `mr_data` plus log-normalised expression;
- `raw_df.rds`: `mr_data` plus counts;
- `data_source.csv`;
- `gene_annotation.rds`.

Gene columns use PlasmoDB IDs (`PVP01_...`).

## 60 Assess integration

Compares the PCA before and after Harmony (PCs 1–30, Euclidean distances, 30
nearest neighbours):
- the fraction of neighbours from another study and from another stage, and
  the silhouette width by study and by stage, on 8,000 random cells (seed 1);
- for each stage present in more than one study, on up to 3,000 cells of the
  stage: the fraction of same-stage neighbours from another study, divided
  by its value under complete mixing, Σ pᵢ (n − nᵢ)/(n − 1) (n cells used,
  nᵢ of them from study i, pᵢ = nᵢ/n): the neighbours of a cell are then a
  random draw of the other n − 1 cells.

Results in `results/integration/`: count tables, `integration_metrics.tsv`,
`per_stage_sharing.tsv` and `per_stage_pairwise_sharing.tsv`.

## 70 Tables and values

Collects the numbers of the earlier stages into `results/tables_and_values/`,
ready to be cited:
- `values.tsv` (one row per number, with its source file), and
  `values.tex` with one LaTeX command per value;
- the dataset and stage-count tables (`.tsv`, and the table rows as `.tex`);
- the STARsolo options read from the run logs;
- the parameters and the software versions;
- the column dictionary.

The manuscript is not changed by the pipeline: it is edited separately,
using these files.

## 80 Figures

Draws the figures of the manuscript into `results/figures/`, one PDF (for
the manuscript) and one PNG (for a quick look) each, named by their number
in the manuscript:

| File | Content |
|---|---|
| `fig1_overview` | workflow; cells per study; stage composition of each study |
| `fig2_blood_stages` | blood-stage UMAP: stage and gametocyte state, UMIs, female and male marker scores; curated marker panel; our calls in the cells labelled by the source studies |
| `fig3_integration` | UMAP before and after Harmony, by study and stage; mixing metrics; cross-study sharing within stage |
| `fig4_invasion_genes` | DBP, CSP and AMA1 on the integrated UMAP, and the fraction of cells detecting them per stage and study |
| `s1_fig_qc` | genes and UMIs per cell, per study |
| `s2_fig_stage_call` | SingleR delta per reference time point; share of labels flagged by pruning |
| `s3_fig_blood_clusters` | blood-stage clusters and their stage |
| `s4_fig_gametocyte_calls` | how clusters are called: marker enrichment, relative UMIs, thresholds |
| `s5_fig_technical` | platform and UMIs on the embeddings |

The numbers printed in the figures (cells, genes) come from
`results/tables_and_values/values.tsv`, so stage 70 runs first.

The style is shared (`R/lib/plots.R`): one colour and one display name per
study and per stage, Helvetica, at most 7.5 x 8.75 in (PLOS). The stage
colours are chosen so that every pair stays distinguishable with
protanopia and deuteranopia. Each figure subject has its own module in
`R/figures/`; to add a figure, write a function there that returns a
ggplot/patchwork object and save it with `save_figure()` in
`R/80_figures.R`.

Fig 5 (the app, side-by-side coloring) is a screenshot, taken with the
screenshot tool of the app repository (`tools/screenshots/manuscript.yml`).

