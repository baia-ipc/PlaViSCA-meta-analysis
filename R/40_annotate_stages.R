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
#  4. Final stage, first available of: the stage given by the authors
#     (author_stage), the cluster gametocyte call, the stage from the
#     reference time course. stage_source says which one was used, and
#     broad_stage is the broad stage of metadata/stages.tsv.
#
# Writes atlas.rds and, in results/stage_annotation/, the tables behind the
# numbers reported for the annotation.

STAGE <- "40_annotate_stages"
source("R/lib/setup.R")
source("R/lib/gametocytes.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
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

pdf(file.path(res, "idc_delta_distribution.pdf"), width = 10, height = 6)
print(plotDeltaDistribution(prediction, ncol = 4))
dev.off()
record_output(file.path(res, "idc_delta_distribution.pdf"))

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

blood_only <- c("idc_timepoint", "cluster_gametocyte_call", "gametocyte_state", "blood_cluster")
for (column in blood_only) {
  check(all(is.na(meta[[column]][!eligible])), column, " is set for a cell not eligible for blood-stage inference")
}
gametocyte <- meta$stage %in% c("Developing gametocyte", "Female gametocyte", "Male gametocyte")
check(all(meta$tissue[gametocyte] == "Host blood"), "gametocyte stage outside host blood")
check(all(meta$broad_stage[meta$tissue == "Host liver"] %in% c(NA, "Liver stage")), "a liver cell has a non-liver stage")

# --- figures: gametocytes on the blood-only UMAP ---------------------------
# (a) female and (b) male marker scores, (c) UMIs, (d) gametocyte states;
# and the curated marker panel per state.

blood_meta <- meta[eligible, ]
blood_meta$umap_b_1 <- blood_cells$umap_b_1[match(rownames(blood_meta), blood_cells$cell)]
blood_meta$umap_b_2 <- blood_cells$umap_b_2[match(rownames(blood_meta), blood_cells$cell)]
state_levels <- c("Developing gametocyte", "Early female", "Late female", "Early male", "Late male")
blood_levels <- c("Ring/Early trophozoite", "Trophozoite", "Schizont (Blood stage)", "Merozoite", state_levels)
blood_colours <- c("Ring/Early trophozoite" = "#A9D18E", "Trophozoite" = "#F4D35E", "Schizont (Blood stage)" = "#8EC5E8",
  "Merozoite" = "#2E75B6", "Developing gametocyte" = "#C0504D", "Early female" = "#759F8B",
  "Late female" = "#4E6E59", "Early male" = "#7F88B3", "Late male" = "#42517C")
blood_labels <- c("Ring/Early trophozoite" = "Ring/early trophozoite", "Trophozoite" = "Trophozoite", "Schizont (Blood stage)" = "Schizont",
  "Merozoite" = "Late schizont / merozoite", setNames(state_levels, state_levels))
set.seed(1)
blood_meta <- blood_meta[sample(nrow(blood_meta)), ]
umap_panel <- function(colour, title, continuous = TRUE, point_size = 0.2) {
  g <- ggplot(blood_meta, aes(umap_b_1, umap_b_2, colour = .data[[colour]])) +
    ggrastr::rasterise(geom_point(size = point_size, stroke = 0, shape = 16), dpi = 300, dev = "ragg") +
    coord_equal() + labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    theme_classic(base_size = 8, base_family = "Helvetica") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 9, face = "bold"))
  if (continuous) g + scale_colour_viridis_c(limits = quantile(blood_meta[[colour]], c(0.01, 0.99)), oob = scales::squish)
  else g
}
# blood stage: the gametocyte state where there is one, else the stage
blood_meta$blood_stage <- factor(coalesce(blood_meta$gametocyte_state, blood_meta$stage), levels = blood_levels)
blood_meta$log10_umi <- log10(blood_meta$nCount_RNA)
figure <- (umap_panel("female_gam_score", "Female marker score") | umap_panel("male_gam_score", "Male marker score")) /
  (umap_panel("log10_umi", "UMIs (log10)") |
    (umap_panel("blood_stage", "Blood stage", continuous = FALSE) +
      scale_colour_manual(values = blood_colours, labels = blood_labels, na.value = "#D9D9D9", na.translate = FALSE, drop = FALSE) +
      guides(colour = guide_legend(override.aes = list(size = 2.5))))) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(res, "gametocytes_umap.pdf"), figure, width = 7.5, height = 7, device = cairo_pdf)
record_output(file.path(res, "gametocytes_umap.pdf"))
ggsave(file.path(res, "gametocytes_umap.png"), figure, width = 7.5, height = 7, dpi = 300)
record_output(file.path(res, "gametocytes_umap.png"))

# supplementary: (a) the blood-stage clusters of stage 35, labelled at their
# median position, and (b) the blood stage
cluster_labels <- blood_meta |>
  group_by(blood_cluster) |>
  summarise(umap_b_1 = median(umap_b_1), umap_b_2 = median(umap_b_2), .groups = "drop")
label_layer <- geom_label(data = cluster_labels, aes(umap_b_1, umap_b_2, label = blood_cluster), inherit.aes = FALSE,
  size = 2.2, label.padding = unit(0.08, "lines"), alpha = 0.8)
blood_meta$blood_cluster <- factor(blood_meta$blood_cluster, levels = sort(as.integer(unique(blood_meta$blood_cluster))))
cluster_figure <- (umap_panel("blood_cluster", sprintf("Blood-stage clusters (resolution %s)", PARAMS$blood_embedding$cluster_resolution),
    continuous = FALSE, point_size = 0.5) + guides(colour = "none") + label_layer) |
  (umap_panel("blood_stage", "Blood stage", continuous = FALSE, point_size = 0.5) + label_layer +
    scale_colour_manual(values = blood_colours, labels = blood_labels, na.value = "#D9D9D9", na.translate = FALSE, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(res, "blood_clusters_umap.pdf"), cluster_figure, width = 7.5, height = 4, device = cairo_pdf)
record_output(file.path(res, "blood_clusters_umap.pdf"))
ggsave(file.path(res, "blood_clusters_umap.png"), cluster_figure, width = 7.5, height = 4, dpi = 300)
record_output(file.path(res, "blood_clusters_umap.png"))

panel <- markers[startsWith(markers$marker_set, "panel_"), ]
panel_data <- GetAssayData(atlas, layer = "data")[panel$gene_id, rownames(blood_meta)]
group <- factor(ifelse(is.na(blood_meta$gametocyte_state), "Asexual", blood_meta$gametocyte_state), levels = c("Asexual", state_levels))
dots <- do.call(rbind, lapply(levels(group), function(g) {
  x <- panel_data[, group == g, drop = FALSE]
  data.frame(state = g, gene_id = panel$gene_id, gene = panel$gene, marker_set = sub("^panel_", "", panel$marker_set),
    detected = Matrix::rowMeans(x > 0), mean = Matrix::rowMeans(x))
}))
write_table(dots, file.path(res, "gametocyte_marker_panel.tsv"))
dots$label <- factor(paste0(dots$gene, " (", sub("PVP01-", "", dots$gene_id), ")"), levels = rev(unique(paste0(panel$gene, " (", sub("PVP01-", "", panel$gene_id), ")"))))
dots$state <- factor(dots$state, levels = levels(group))
dotplot <- ggplot(dots, aes(state, label, size = detected, colour = mean)) + geom_point() +
  scale_colour_viridis_c() + scale_size_area(max_size = 4) +
  facet_grid(marker_set ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, size = "Fraction detected", colour = "Mean expression") +
  theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
for (ext in c("pdf", "png")) {
  ggsave(file.path(res, paste0("gametocyte_marker_panel.", ext)), dotplot, width = 5.5, height = 1.5 + 0.2 * nrow(panel), dpi = 300)
  record_output(file.path(res, paste0("gametocyte_marker_panel.", ext)))
}

atlas@meta.data <- meta
save_object(atlas, file.path(out, "atlas.rds"))
save_session_info(out)
