# Stage 40: life-cycle stage of each cell.
#
#  1. Blood-stage cells eligible for inference: cells of host blood samples
#     with at least idc_min_umi UMIs (config/params.yaml).
#  2. Asexual stage by similarity to a reference time course (SingleR):
#     the reference is the SMRU1 time course of Zhu et al. 2016 (RNA-seq,
#     7 time points), with genes mapped to PvP01 by the PlasmoDB one-to-one
#     orthologs. Each eligible cell gets the time point it most resembles
#     (idc_timepoint, SingleR label; idc_delta, SingleR delta.next;
#     idc_pruned, TRUE where SingleR's pruning would discard the label),
#     which is converted to a stage with metadata/idc_timepoints.tsv
#     (idc_stage).
#  3. Gametocytes, per cluster of the blood-only clustering of stage 35
#     (blood_cluster): the detection of the marker sets of
#     metadata/gametocyte_markers.tsv in each cluster, against control genes
#     of similar expression (R/lib/gametocytes.R), and the cluster's median
#     UMIs relative to the other clusters (UMIs grow as gametocytes mature)
#     assign each cluster to the female or male branch
#     (cluster_gametocyte_call) and to a state along it (gametocyte_state:
#     developing gametocyte, early/late female, early/late male), by the rule
#     documented in config/params.yaml (annotation). Module scores of the
#     female and male sets (female_gam_score, male_gam_score) are kept per
#     cell for display.
#  4. Final stage. Eligible cells take the cluster gametocyte call, else
#     the stage from the reference time course; all other cells take the
#     stage given by the authors (author_stage). The authors' gametocyte
#     labels of Sa2020 (author_sex) and Hazzard2024's sexual/asexual groups
#     (author_pseudogroup) stay as fields for comparison. stage_source says
#     which source was used, and broad_stage is the broad stage of
#     metadata/stages.tsv.
#
# Writes atlas.rds and, in results/stage_annotation/, the tables behind the
# numbers and figures reported for the annotation (the figures themselves are
# made by stage 80).

STAGE <- "40_annotate_stages"
source("R/lib/setup.R")
source("R/lib/gametocytes.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(SummarizedExperiment)
  library(readxl)
})

out <- stage_dir(STAGE)
res <- results_dir("stage_annotation")
p <- PARAMS$annotation
stages <- read_meta("stages.tsv")

atlas <- readRDS(file.path(stage_dir("30_integrate"), "atlas_integrated.rds"))
meta <- atlas@meta.data

# --- 1. eligible cells ------------------------------------------------------

eligible <- meta$tissue == "Host blood" & meta$nCount_RNA >= p$idc_min_umi
eligibility <- meta |>
  mutate(eligible = eligible) |>
  dplyr::count(study_label, tissue, eligible, name = "cells")
write_table(eligibility, file.path(res, "eligible_cells.tsv"))
message(sprintf("%d of %d cells eligible for blood-stage inference", sum(eligible), ncol(atlas)))

# --- 2. reference time course -----------------------------------------------

zhu <- as.data.frame(read_excel(input_path(p$idc_reference), sheet = p$idc_reference_sheet))
orthologs <- as.data.frame(read_excel(input_path(p$idc_orthologs)))
check_unique(orthologs$InputOrtholog, "PVX genes in the ortholog table")
check_unique(orthologs$GeneID, "PvP01 genes in the ortholog table")
timepoints <- setdiff(colnames(zhu), "Gid")
reference <- inner_join(zhu, orthologs[, c("InputOrtholog", "GeneID")], by = c("Gid" = "InputOrtholog"))
check(!anyNA(reference[, timepoints]), "missing values in the reference time course")
reference_matrix <- as.matrix(reference[, timepoints])
rownames(reference_matrix) <- reference$GeneID
reference_summary <- data.frame(
  reference_genes = nrow(zhu), with_pvp01_ortholog = nrow(reference),
  without_ortholog = nrow(zhu) - nrow(reference), timepoints = paste(timepoints, collapse = ","))
message(sprintf("reference: %d of %d genes have a PvP01 ortholog", nrow(reference), nrow(zhu)))

test <- as.SingleCellExperiment(DietSeurat(subset(atlas, cells = colnames(atlas)[eligible]), layers = "data", assays = "RNA"))
reference_summary$genes_shared_with_atlas <- length(intersect(rownames(test), rownames(reference_matrix)))
write_table(reference_summary, file.path(res, "idc_reference.tsv"))

reference_se <- SummarizedExperiment(assays = list(logcounts = reference_matrix))
prediction <- SingleR(test = test, ref = reference_se, labels = colnames(reference_matrix), assay.type.test = 1)
rm(test)

meta$idc_timepoint <- NA_character_
meta$idc_delta <- NA_real_
meta$idc_pruned <- NA
# every eligible cell keeps its best time point; the cells whose label
# SingleR's pruning would discard are flagged
meta[eligible, "idc_timepoint"] <- prediction$labels
meta[eligible, "idc_delta"] <- prediction$delta.next
meta[eligible, "idc_pruned"] <- is.na(prediction$pruned.labels)
message(sprintf("%d of %d eligible cells have a label that SingleR's pruning would discard", sum(is.na(prediction$pruned.labels)), sum(eligible)))
timepoint_stage <- read_meta("idc_timepoints.tsv")
check_same_set(timepoint_stage$timepoint_h, timepoints, "time points of metadata/idc_timepoints.tsv and of the reference")
meta$idc_stage <- timepoint_stage$stage[match(meta$idc_timepoint, timepoint_stage$timepoint_h)]

# --- 3. gametocytes ---------------------------------------------------------

atlas_cluster <- paste0("cluster_res_", sub(".", "_", formatC(p$atlas_cluster_resolution, format = "f", digits = 2), fixed = TRUE))
meta$cluster <- meta[[atlas_cluster]]

blood_cells <- utils::read.delim(file.path(stage_dir("35_blood_embedding"), "blood_cells.tsv"), stringsAsFactors = FALSE)
check_same_set(blood_cells$cell, rownames(meta)[eligible], "cells of stage 35 and cells eligible for blood-stage inference")
meta$blood_cluster <- NA_character_
meta$blood_cluster[eligible] <- as.character(blood_cells$blood_cluster[match(rownames(meta)[eligible], blood_cells$cell)])

markers <- read_meta("gametocyte_markers.tsv")
check(all(markers$gene_id %in% rownames(atlas)), "gametocyte marker genes missing from the atlas")
module_score <- function(set) {
  scored <- AddModuleScore(atlas, features = list(markers$gene_id[markers$marker_set == set]), name = "score")
  scored$score1
}
meta$female_gam_score <- module_score("female")
meta$male_gam_score <- module_score("male")

marker_sets <- list(female = markers$gene_id[markers$marker_set == "female"],
  male = markers$gene_id[markers$marker_set == "male"],
  male_mature = markers$gene_id[markers$marker_set == "panel_male"])
enrichment <- marker_enrichment(GetAssayData(atlas, layer = "counts")[, eligible], GetAssayData(atlas, layer = "data")[, eligible],
  meta$blood_cluster[eligible], marker_sets, p$marker_expression_bins, p$marker_control_genes, p$marker_control_seed)
cluster_umi <- tapply(meta$nCount_RNA[eligible], meta$blood_cluster[eligible], median)
cluster_scores <- data.frame(blood_cluster = rownames(enrichment), cells = as.integer(table(meta$blood_cluster[eligible])[rownames(enrichment)]),
  enrichment, median_umi = as.numeric(cluster_umi[rownames(enrichment)])) |>
  mutate(relative_umi = median_umi / median(median_umi)) |>
  call_gametocyte_clusters(p) |>
  arrange(call, state, desc(female))
print(as.data.frame(cluster_scores), digits = 2)
write_table(cluster_scores, file.path(res, "cluster_gametocyte_scores.tsv"))
check(any(cluster_scores$call %in% "Female gametocyte") && any(cluster_scores$call %in% "Male gametocyte"),
  "the rule finds no female or no male cluster; see cluster_gametocyte_scores.tsv")

meta$cluster_gametocyte_call <- NA_character_
meta$gametocyte_state <- NA_character_
meta$cluster_gametocyte_call[eligible] <- cluster_scores$call[match(meta$blood_cluster[eligible], cluster_scores$blood_cluster)]
meta$gametocyte_state[eligible] <- cluster_scores$state[match(meta$blood_cluster[eligible], cluster_scores$blood_cluster)]

# --- 4. final stage ---------------------------------------------------------

# our own call for the host-blood cells we can classify, the authors'
# stage for the others
meta$stage <- ifelse(eligible, coalesce(meta$cluster_gametocyte_call, meta$idc_stage), meta$author_stage)
meta$stage_source <- case_when(
  eligible & !is.na(meta$cluster_gametocyte_call) ~ "gametocyte_markers",
  eligible & !is.na(meta$idc_stage) ~ "idc_reference",
  !eligible & !is.na(meta$author_stage) ~ "author",
  TRUE ~ NA_character_
)
check(all(meta$stage %in% c(NA, stages$stage)), "stages not in metadata/stages.tsv")
meta$broad_stage <- stages$broad_stage[match(meta$stage, stages$stage)]

# comparison of the authors' gametocyte labels with the cluster call (all
# cells of the two studies with author sexual-stage labels)
author_vs_call <- meta |>
  filter(study_label %in% c("Sa2020", "Hazzard2024")) |>
  transmute(study_label,
    author_label = coalesce(author_stage, author_pseudogroup),
    cluster_gametocyte_call = coalesce(cluster_gametocyte_call, "none")) |>
  dplyr::count(study_label, author_label, cluster_gametocyte_call, name = "cells")
write_table(author_vs_call, file.path(res, "author_labels_vs_cluster_call.tsv"))

# --- checks -----------------------------------------------------------------

blood_only <- c("idc_timepoint", "cluster_gametocyte_call", "gametocyte_state", "blood_cluster")
for (column in blood_only) {
  check(all(is.na(meta[[column]][!eligible])), column, " is set for a cell not eligible for blood-stage inference")
}
gametocyte <- meta$stage %in% c("Developing gametocyte", "Female gametocyte", "Male gametocyte")
check(all(meta$tissue[gametocyte] == "Host blood"), "gametocyte stage outside host blood")
check(all(meta$broad_stage[meta$tissue == "Host liver"] %in% c(NA, "Liver stage")), "a liver cell has a non-liver stage")

# --- curated marker panel per gametocyte state -----------------------------
# fraction of cells detecting each gene of the curated sexual-stage panel and
# its mean normalised expression, in the asexual cells and in each
# gametocyte state (eligible cells only)

state_levels <- c("Developing gametocyte", "Early female", "Late female", "Early male", "Late male")
panel <- markers[startsWith(markers$marker_set, "panel_"), ]
panel_data <- GetAssayData(atlas, layer = "data")[panel$gene_id, eligible]
group <- factor(ifelse(is.na(meta$gametocyte_state[eligible]), "Asexual", meta$gametocyte_state[eligible]),
  levels = c("Asexual", state_levels))
dots <- do.call(rbind, lapply(levels(group), function(g) {
  x <- panel_data[, group == g, drop = FALSE]
  data.frame(state = g, gene_id = panel$gene_id, gene = panel$gene, marker_set = sub("^panel_", "", panel$marker_set),
    detected = Matrix::rowMeans(x > 0), mean = Matrix::rowMeans(x))
}))
write_table(dots, file.path(res, "gametocyte_marker_panel.tsv"))

atlas@meta.data <- meta
save_object(atlas, file.path(out, "atlas.rds"))
save_session_info(out)
