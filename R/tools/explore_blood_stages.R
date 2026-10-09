# Exploration (not a pipeline stage): blood-stage cells in their own space.
#
# Usage: Rscript R/tools/explore_blood_stages.R
#        (PLAVISCA_RUN=<run> selects the run whose atlas is used)
#
# Takes the blood cells eligible for blood-stage inference from the stage 40
# atlas, recomputes variable genes, PCA, Harmony (by study) and a 2D UMAP on
# them only, clusters them at several resolutions, and plots the UMAP coloured
# by stage, reference time point, study, gametocyte marker scores, AP2-G and
# the authors' labels. Writes to <run>/explore/blood_stages/.

STAGE <- "explore_blood_stages"
source("R/lib/setup.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

out <- file.path(RUN_DIR, "explore", "blood_stages")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
resolutions <- c(0.2, 0.5, 1, 2)
ap2g <- "PVP01-1440800"

atlas <- readRDS(file.path(stage_dir("40_annotate_stages"), "atlas.rds"))
eligible <- colnames(atlas)[atlas$tissue == "Host blood" & atlas$nCount_RNA >= PARAMS$annotation$idc_min_umi]
blood <- DietSeurat(subset(atlas, cells = eligible), layers = "counts", assays = "RNA")
rm(atlas)
message(sprintf("%d blood cells: %s", ncol(blood), paste(names(table(blood$study_label)), table(blood$study_label), collapse = ", ")))

blood[["RNA"]] <- split(blood[["RNA"]], f = blood$study_label)
blood <- NormalizeData(blood)
blood <- FindVariableFeatures(blood, nfeatures = 2000)
blood <- ScaleData(blood)
blood <- RunPCA(blood, npcs = 50, seed.use = 42)
blood <- IntegrateLayers(blood, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
blood <- JoinLayers(blood)
blood <- RunUMAP(blood, reduction = "harmony", dims = 1:30, n.components = 2, seed.use = 42, reduction.name = "umap_blood")
blood <- FindNeighbors(blood, reduction = "harmony", dims = 1:30)
for (r in resolutions) blood <- FindClusters(blood, resolution = r, random.seed = 0)

meta <- blood@meta.data
meta$ap2g_umi <- GetAssayData(blood, layer = "counts")[ap2g, ]
meta$idc_timepoint <- factor(meta$idc_timepoint, levels = c("6", "12", "24", "32", "36", "42", "48"))
meta$author_label <- dplyr::coalesce(meta$author_sex, meta$author_pseudogroup)
umap <- Embeddings(blood, "umap_blood")
meta$umap_1 <- umap[, 1]
meta$umap_2 <- umap[, 2]
saveRDS(blood, file.path(out, "blood.rds"))

# --- figures ----------------------------------------------------------------

panel <- function(colour, title, continuous = FALSE) {
  d <- meta[sample(nrow(meta)), ]
  p <- ggplot(d, aes(umap_1, umap_2, colour = .data[[colour]])) +
    ggrastr::rasterise(geom_point(size = 0.05, stroke = 0, shape = 16), dpi = 200, dev = "ragg") +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    theme_classic(base_size = 9) + theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (continuous) p + scale_colour_viridis_c(limits = quantile(d[[colour]], c(0.01, 0.99), na.rm = TRUE), oob = scales::squish)
  else p + guides(colour = guide_legend(override.aes = list(size = 2.5)))
}
set.seed(1)
meta$ap2g_detected <- ifelse(meta$ap2g_umi > 0, "AP2-G detected", "not detected")
overview <- (panel("stage", "Current stage") | panel("idc_timepoint", "Reference time point (h)")) /
  (panel("study_label", "Study") | panel("author_label", "Authors' label")) /
  (panel("female_gam_score", "Female marker score", TRUE) | panel("male_gam_score", "Male marker score", TRUE)) /
  (panel("ap2g_detected", "AP2-G (PVP01_1440800)") | panel("nCount_RNA", "UMIs", TRUE))
ggsave(file.path(out, "blood_umap_overview.png"), overview, width = 11, height = 16, dpi = 150)

clusters <- wrap_plots(lapply(resolutions, function(r) {
  col <- paste0("RNA_snn_res.", r)
  panel(col, sprintf("Clusters, resolution %s (%d)", r, length(unique(meta[[col]])))) + theme(legend.position = "none") +
    geom_text(data = aggregate(cbind(umap_1, umap_2) ~ cl, data.frame(cl = meta[[col]], umap_1 = meta$umap_1, umap_2 = meta$umap_2), median),
      aes(umap_1, umap_2, label = cl), inherit.aes = FALSE, size = 2.5)
}), ncol = 2)
ggsave(file.path(out, "blood_umap_clusters.png"), clusters, width = 11, height = 11, dpi = 150)

# --- per-cluster summaries --------------------------------------------------

for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  s <- meta |>
    group_by(cluster = .data[[col]]) |>
    summarise(cells = n(),
      mean_female_score = mean(female_gam_score), mean_male_score = mean(male_gam_score),
      ap2g_detected = mean(ap2g_umi > 0), median_umi = median(nCount_RNA),
      main_timepoint = names(which.max(table(idc_timepoint))),
      current_stage = names(which.max(table(stage))),
      authors_female_or_sexual = sum(author_label %in% c("Female gametocyte", "Sexual")),
      authors_male = sum(author_label %in% "Male gametocyte"),
      Sa2020 = sum(study_label == "Sa2020"), Hazzard2022 = sum(study_label == "Hazzard2022"), Hazzard2024 = sum(study_label == "Hazzard2024"),
      .groups = "drop")
  write.table(s, file.path(out, sprintf("clusters_res_%s.tsv", r)), sep = "\t", quote = FALSE, row.names = FALSE)
}
message("written to ", out)
