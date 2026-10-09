# Stage 50: data files for the PlaViSca app.
#
# Writes to results/app_data/:
#   cleaned_dataset.rds   list of
#                           mr_data  per-cell metadata and embeddings (the
#                                    columns of metadata/columns.tsv)
#                           pca_df   per study: standard deviation of each PC
#                           hiv_data per study: variable genes (vst)
#   normalize_df.rds      mr_data + log-normalised expression (cells x genes)
#   raw_df.rds            mr_data + raw UMI counts (cells x genes)
#   data_source.csv       one row per study (metadata/studies.tsv + cell numbers)
#   gene_annotation.rds   gene ID, name and description (stage 00)
#
# Cells are sorted by name. Gene IDs are written as in PlasmoDB ("PVP01_...").

STAGE <- "50_export_app_data"
source("R/lib/setup.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyr)
})

out <- results_dir("app_data")
atlas <- readRDS(file.path(stage_dir("40_annotate_stages"), "atlas.rds"))
columns <- read_meta("columns.tsv")
gene_annotation <- readRDS(file.path(stage_dir("00_reference"), "gene_annotation.rds"))
plasmodb_id <- function(seurat_id) gsub("-", "_", seurat_id)

# --- per-cell table ---------------------------------------------------------

embedding <- function(reduction, prefix) {
  x <- Embeddings(atlas, reduction)[, 1:3]
  colnames(x) <- paste0(prefix, 1:3)
  x
}
mr_data <- atlas@meta.data |>
  mutate(barcode = colnames(atlas), n_count_rna = nCount_RNA, n_feature_rna = nFeature_RNA) |>
  cbind(embedding("umap_unintegrated", "umap_u_"), embedding("umap_integrated", "umap_i_"),
    embedding("pca_unintegrated", "pca_u_"), embedding("pca_integrated", "pca_i_"),
    embedding("tsne_integrated", "t_sne_"))
# blood-stage UMAP (stage 35), only for the cells of the blood-stage analysis
blood_cells <- utils::read.delim(file.path(stage_dir("35_blood_embedding"), "blood_cells.tsv"), stringsAsFactors = FALSE)
mr_data$umap_b_1 <- blood_cells$umap_b_1[match(rownames(mr_data), blood_cells$cell)]
mr_data$umap_b_2 <- blood_cells$umap_b_2[match(rownames(mr_data), blood_cells$cell)]
missing <- setdiff(columns$column, colnames(mr_data))
check(length(missing) == 0, "columns of metadata/columns.tsv missing from the atlas: ", paste(missing, collapse = ", "))
mr_data <- mr_data[sort(colnames(atlas)), columns$column]
check_unique(mr_data$barcode, "cell names")

# --- expression tables ------------------------------------------------------

expression_table <- function(layer) {
  expression <- Matrix::t(GetAssayData(atlas, assay = "RNA", layer = layer))[rownames(mr_data), ]
  colnames(expression) <- plasmodb_id(colnames(expression))
  check_same_order(rownames(expression), rownames(mr_data), "cells of the expression table and of mr_data")
  cbind(mr_data, as.matrix(expression))
}
save_object(expression_table("data"), file.path(out, "normalize_df.rds"))
save_object(expression_table("counts"), file.path(out, "raw_df.rds"))

# --- per-study PCA and variable genes (QC panels of the app) ----------------

per_study <- lapply(SplitObject(atlas, split.by = "study_label"), function(object) {
  object <- NormalizeData(object, assay = "RNA", verbose = FALSE)
  object <- RunPCA(object, assay = "RNA", verbose = FALSE)
  object <- FindVariableFeatures(object, assay = "RNA", selection.method = "vst", nfeatures = 2000)
  study <- object$study_label[1]
  list(
    pca = data.frame(pc = seq_along(object@reductions$pca@stdev), std = object@reductions$pca@stdev, study = study),
    variable = HVFInfo(object, assay = "RNA", method = "vst", status = TRUE) |>
      tibble::rownames_to_column("id") |>
      mutate(id = plasmodb_id(id)) |>
      left_join(gene_annotation, by = "id") |>
      transmute(id, name = coalesce(name, id), mean, variable, rank, variance.standardized, study = study)
  )
})
check_same_set(names(per_study), STUDIES, "studies in the per-study tables")

save_object(list(
  mr_data = mr_data,
  pca_df = bind_rows(lapply(per_study, `[[`, "pca")),
  hiv_data = bind_rows(lapply(per_study, `[[`, "variable"))
), file.path(out, "cleaned_dataset.rds"))

# --- data sources -----------------------------------------------------------

Sys.setlocale("LC_TIME", "C")
studies <- read_meta("studies.tsv")
cells <- table(mr_data$study_label)
data_source <- data.frame(
  Publication_date = sub("^0", "", format(as.Date(studies$publication_date), "%d-%b-%y")),
  Study_label = studies$study_label,
  study_pmid = studies$study_pmid,
  Authors = studies$authors,
  Journal = studies$journal,
  Title = gsub("Plasmodium vivax", "<em>Plasmodium vivax</em>", studies$title, fixed = TRUE),
  DOI = sprintf("<a href='https://doi.org/%s' target='_blank'>https://doi.org/%s</a>", studies$doi, studies$doi),
  Number_of_cells = as.integer(cells[studies$study_label])
)
readr::write_csv(data_source, file.path(out, "data_source.csv"))
record_output(file.path(out, "data_source.csv"))

save_object(gene_annotation, file.path(out, "gene_annotation.rds"))
