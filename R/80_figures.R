# Stage 80: figures of the manuscript.
#
# Reads the atlas (stage 40), the blood-stage embedding (stage 35) and the
# result tables of stages 40 and 60, and draws every figure with the shared
# style of R/lib/plots.R. One module per subject in R/figures/:
#   overview.R        workflow and cells per study and stage
#   blood_stages.R    blood-stage embedding, gametocyte calls and markers
#   integration.R     embeddings before and after Harmony, mixing metrics
#   invasion_genes.R  DBP, CSP and AMA1 across stages and studies
#   quality.R         UMIs and genes per cell, confidence of the stage call
#
# Writes to results/figures/ one PDF (for the manuscript) and one PNG (for a
# quick look) per figure, named by their number in the manuscript:
#   fig1_overview          fig2_blood_stages      fig3_integration
#   fig4_invasion_genes    s1_fig_qc              s2_fig_stage_call
#   s3_fig_blood_clusters  s4_fig_gametocyte_calls  s5_fig_technical
# Fig 5 (screenshots of the app) is made in the app repository
# (tools/screenshots/).

STAGE <- "80_figures"
source("R/lib/setup.R")
source("R/lib/plots.R")
suppressPackageStartupMessages(library(Seurat))
for (module in list.files("R/figures", pattern = "\\.R$", full.names = TRUE)) source(module)

FIGURES_DIR <- results_dir("figures")
read_result <- function(...) utils::read.delim(file.path(RESULTS, ...), stringsAsFactors = FALSE, check.names = FALSE)

# --- data -------------------------------------------------------------------

genes <- utils::read.delim(file.path(stage_dir("00_reference"), "genes.tsv"), stringsAsFactors = FALSE)
case_genes <- c("DBP", "CSP", "AMA1")
check(all(case_genes %in% genes$name) && !anyDuplicated(genes$name[genes$name %in% case_genes]),
  "each case-study gene must have exactly one gene of that name in genes.tsv")
case_gene_ids <- setNames(genes$seurat_id[match(case_genes, genes$name)], case_genes)

atlas <- readRDS(file.path(stage_dir("40_annotate_stages"), "atlas.rds"))
meta <- atlas@meta.data
meta$cell <- rownames(meta)
umap <- list(
  unintegrated = Embeddings(atlas, "umap_unintegrated")[, 1:2],
  integrated = Embeddings(atlas, "umap_integrated")[, 1:2])
check_same_order(rownames(umap$integrated), meta$cell, "cells of the UMAP and of the metadata")
expression <- as.matrix(t(GetAssayData(atlas, layer = "data")[case_gene_ids, ]))
colnames(expression) <- case_genes
rm(atlas)
invisible(gc())

check(all(meta$study_label %in% names(STUDY_COLOURS)), "a study has no figure colour (R/lib/plots.R)")
check(all(meta$stage %in% c(NA, names(STAGE_COLOURS))), "a stage has no figure colour (R/lib/plots.R)")
check_same_set(read_meta("stages.tsv")$stage, names(STAGE_COLOURS), "stages of metadata/stages.tsv and of the figure colours")

blood_cells <- utils::read.delim(file.path(stage_dir("35_blood_embedding"), "blood_cells.tsv"), stringsAsFactors = FALSE)
blood <- meta[match(blood_cells$cell, meta$cell), ]
blood$x <- blood_cells$umap_b_1
blood$y <- blood_cells$umap_b_2
blood$state <- dplyr::coalesce(blood$gametocyte_state, blood$stage)
check(all(blood$state %in% names(STATE_COLOURS)), "a blood-stage cell has no state colour")

studies <- read_meta("studies.tsv")
cluster_scores <- read_result("stage_annotation", "cluster_gametocyte_scores.tsv")
marker_panel <- read_result("stage_annotation", "gametocyte_marker_panel.tsv")
author_calls <- read_result("stage_annotation", "author_labels_vs_cluster_call.tsv")
metrics <- read_result("integration", "integration_metrics.tsv")
sharing <- read_result("integration", "per_stage_sharing.tsv")

# the numbers cited in the text (stage 70), so that figures and text agree
values <- read_result("tables_and_values", "values.tsv")
value <- function(key) {
  check(key %in% values$key, "value ", key, " missing from values.tsv")
  values$formatted[values$key == key]
}

# --- figures ----------------------------------------------------------------

save_figure(figure_overview(meta, studies, value), "fig1_overview", 7.5, 4.8)
save_figure(figure_blood_stages(blood, marker_panel, author_calls), "fig2_blood_stages", 7.5, 8.6)
save_figure(figure_integration(meta, umap, metrics, sharing), "fig3_integration", 7.5, 8.4)
save_figure(figure_invasion_genes(meta, umap$integrated, expression), "fig4_invasion_genes", 7.5, 7.4)

save_figure(figure_qc(meta), "s1_fig_qc", 7.5, 3.2)
save_figure(figure_stage_call(meta), "s2_fig_stage_call", 7.5, 3.4)
save_figure(figure_blood_clusters(blood), "s3_fig_blood_clusters", 7.5, 4)
save_figure(figure_gametocyte_calls(cluster_scores, PARAMS$annotation), "s4_fig_gametocyte_calls", 7.5, 3.4)
save_figure(figure_technical(meta, umap, studies), "s5_fig_technical", 7.5, 3)

save_session_info(FIGURES_DIR)
