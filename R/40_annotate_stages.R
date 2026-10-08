# Stage 40: life-cycle stage of each cell.
#
#  1. Blood-stage cells eligible for inference: cells of host blood samples
#     with at least idc_min_umi UMIs (config/params.yaml).
#  2. Asexual stage by similarity to a reference time course (SingleR):
#     the reference is the SMRU1 time course of Zhu et al. 2016 (RNA-seq,
#     7 time points), with genes mapped to PvP01 by the PlasmoDB one-to-one
#     orthologs. Each eligible cell gets the time point it most resembles
#     (idc_timepoint, pruned SingleR label; idc_delta, SingleR delta.next),
#     which is converted to a stage with metadata/idc_timepoints.tsv
#     (idc_stage).
#  3. Gametocytes by cluster: module scores of female and male markers
#     (metadata/gametocyte_markers.tsv) are averaged over the eligible cells
#     of each cluster; a cluster is called female (male) gametocytes if its
#     mean female (male) score is at least gametocyte_min_score and exceeds
#     the other score by at least gametocyte_min_margin
#     (cluster_gametocyte_call, set for eligible cells only).
#  4. Final stage, first available of: the stage given by the authors
#     (author_stage), the cluster gametocyte call, the stage from the
#     reference time course. stage_source says which one was used, and
#     broad_stage is the broad stage of metadata/stages.tsv.
#
# Writes atlas.rds and, in results/stage_annotation/, the tables behind the
# numbers reported for the annotation.

STAGE <- "40_annotate_stages"
source("R/lib/setup.R")
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
meta[eligible, "idc_timepoint"] <- prediction$pruned.labels
meta[eligible, "idc_delta"] <- prediction$delta.next
timepoint_stage <- read_meta("idc_timepoints.tsv")
check_same_set(timepoint_stage$timepoint_h, timepoints, "time points of metadata/idc_timepoints.tsv and of the reference")
meta$idc_stage <- timepoint_stage$stage[match(meta$idc_timepoint, timepoint_stage$timepoint_h)]

pdf(file.path(res, "idc_delta_distribution.pdf"), width = 10, height = 6)
print(plotDeltaDistribution(prediction, ncol = 4))
dev.off()
record_output(file.path(res, "idc_delta_distribution.pdf"))

# --- 3. gametocyte clusters -------------------------------------------------

cluster_column <- paste0("cluster_res_", sub(".", "_", formatC(p$cluster_resolution, format = "f", digits = 2), fixed = TRUE))
meta$cluster <- meta[[cluster_column]]

markers <- read_meta("gametocyte_markers.tsv")
check(all(markers$gene_id %in% rownames(atlas)), "gametocyte marker genes missing from the atlas")
atlas <- AddModuleScore(atlas, features = list(markers$gene_id[markers$sex == "female"]), name = "female_gam_score")
atlas <- AddModuleScore(atlas, features = list(markers$gene_id[markers$sex == "male"]), name = "male_gam_score")
meta$female_gam_score <- atlas$female_gam_score1
meta$male_gam_score <- atlas$male_gam_score1

cluster_scores <- meta[eligible, ] |>
  group_by(cluster) |>
  summarise(eligible_cells = n(), mean_female_score = mean(female_gam_score), mean_male_score = mean(male_gam_score), .groups = "drop") |>
  mutate(call = case_when(
    mean_female_score >= p$gametocyte_min_score & mean_female_score - mean_male_score >= p$gametocyte_min_margin ~ "Female gametocyte",
    mean_male_score >= p$gametocyte_min_score & mean_male_score - mean_female_score >= p$gametocyte_min_margin ~ "Male gametocyte",
    TRUE ~ NA_character_
  ))
print(as.data.frame(cluster_scores))
write_table(cluster_scores, file.path(res, "cluster_gametocyte_scores.tsv"))
check(any(!is.na(cluster_scores$call)), "no cluster passes the gametocyte thresholds; see cluster_gametocyte_scores.tsv")

meta$cluster_gametocyte_call <- NA_character_
meta$cluster_gametocyte_call[eligible] <- cluster_scores$call[match(meta$cluster[eligible], cluster_scores$cluster)]

# --- 4. final stage ---------------------------------------------------------

meta$stage <- coalesce(meta$author_stage, meta$cluster_gametocyte_call, meta$idc_stage)
meta$stage_source <- case_when(
  !is.na(meta$author_stage) ~ "author",
  !is.na(meta$cluster_gametocyte_call) ~ "gametocyte_markers",
  !is.na(meta$idc_stage) ~ "idc_reference",
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

blood_only <- c("idc_timepoint", "cluster_gametocyte_call")
for (column in blood_only) {
  check(all(is.na(meta[[column]][!eligible])), column, " is set for a cell not eligible for blood-stage inference")
}
gametocyte <- meta$stage %in% c("Female gametocyte", "Male gametocyte")
check(all(meta$tissue[gametocyte] == "Host blood"), "gametocyte stage outside host blood")
check(all(meta$broad_stage[meta$tissue == "Host liver"] %in% c(NA, "Liver stage")), "a liver cell has a non-liver stage")
stage_counts <- meta |>
  mutate(stage = coalesce(stage, "(none)")) |>
  dplyr::count(study_label, stage, name = "cells")
check_expected_counts(stage_counts, "stage_cells")

atlas@meta.data <- meta
save_object(atlas, file.path(out, "atlas.rds"))
save_session_info(out)
