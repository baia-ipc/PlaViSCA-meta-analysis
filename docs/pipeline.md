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
| Sa2020 | filtered | all cells called by STARsolo |
| Hazzard2022 | raw | `emptyDrops` (lower = 1, FDR ≤ 0.001, seed 123456) |
| Hazzard2024 | raw | the cells in the authors' processed table (Zenodo 12775216) |
| Ruberto2022_2 | raw | the cells with Stage = Sporozoite in the authors' object |

For these studies:
- rRNA genes are removed before the objects are made.
- Each run is made into an object with `min.cells`/`min.features` as in
  params, and the runs are merged.
- A barcode-rank (knee) plot is drawn for each run.
- To match a cell with the authors' data, the authors' sample names are
  translated to runs using `metadata/runs.tsv` (`author_sample`).

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
| Sa2020 | female/male gametocytes → `author_stage` | S5 Data |
| Sa2020 | member of the authors' main analysis → `sa2020_main_analysis` | `metadata/sa2020_main_analysis_cells.tsv` |
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

## 40 Annotate stages

1. **Eligible cells.** Cells of host blood samples with at least 10 UMIs
   are eligible for blood-stage inference.
2. **Asexual stage.** SingleR (default settings) compares each eligible
   cell with the SMRU1 time course of Zhu et al. 2016:
   - the reference is RNA-seq with 7 time points (6, 12, 24, 32, 36, 42 and
     48 h);
   - its genes are mapped to PvP01 by PlasmoDB's one-to-one orthologs, and
     genes without an ortholog are dropped.

   `idc_timepoint` is the pruned label: the time point the cell most
   resembles. It is a similarity, not a measured age. `idc_delta` is
   SingleR's `delta.next`. `metadata/idc_timepoints.tsv` turns time points
   into stages (`idc_stage`): 6–12 h ring, 24 h trophozoite, 32–42 h
   schizont, 48 h late schizont/merozoite.
3. **Gametocytes.**
   - Module scores (`AddModuleScore`) are computed for 19 female and 6 male
     marker genes (`metadata/gametocyte_markers.tsv`, from Sa et al. 2020).
   - For each cluster (resolution 0.05), the scores are averaged over its
     eligible cells.
   - A cluster is called female (male) gametocytes if its mean female (male)
     score is at least 0.10 and exceeds the other mean by at least 0.05.
   - Eligible cells take the call of their cluster
     (`cluster_gametocyte_call`).
4. **Final stage.**
   - `stage` takes the first available of: the authors' stage, the cluster
     gametocyte call, `idc_stage`. `stage_source` records which one was
     used.
   - `broad_stage` (blood, liver or sporozoite) comes from
     `metadata/stages.tsv`.
   - The authors' sexual/asexual groups of Hazzard2024 do not enter `stage`;
     they are compared with the cluster call in
     `author_labels_vs_cluster_call.tsv`.

Checks:
- only eligible cells have blood-stage labels;
- gametocytes occur only in host blood;
- liver cells have liver stages;
- cell numbers per study and stage equal `metadata/expected_counts.tsv`.

Results (`results/stage_annotation/`):
- `eligible_cells.tsv`;
- `idc_reference.tsv` (genes of the reference and how many were mapped);
- `idc_delta_distribution.pdf`;
- `cluster_gametocyte_scores.tsv` (mean scores and call per cluster);
- `author_labels_vs_cluster_call.tsv`.

## 50 Export for the app

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
  by its value under complete mixing, 1 − Σ pᵢ² (pᵢ = share of study i in
  the stage).

Results in `results/integration/`: count tables, `integration_metrics.tsv`,
`per_stage_sharing.tsv`, `per_stage_pairwise_sharing.tsv` and the figure
`integration.pdf`.

## 70 Manuscript

Collects the numbers of the earlier stages into `results/manuscript/`:
- `values.tsv`, and `values.tex` with one LaTeX command per value;
- the dataset and stage-count tables (`.tsv` and `.tex`);
- the STARsolo options read from the run logs;
- the parameters and the software versions;
- the column dictionary;
- the supplementary figures.

`--copy-to <dir>` copies the `.tex` and `.pdf` files into the manuscript.
