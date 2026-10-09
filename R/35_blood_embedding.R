# Stage 35: embedding and clusters of the blood-stage cells alone.
#
# The cells eligible for blood-stage inference (host blood, at least
# annotation$idc_min_umi UMIs) are analysed apart from the liver and
# sporozoite cells, so that their own variation (asexual cycle, sexual
# branch) shapes the embedding:
#   - variable genes, PCA and Harmony (one batch per study) on these cells
#   - 2D UMAP (umap_blood) and a nearest-neighbour graph on Harmony PCs
#   - Louvain clusters (blood_cluster), used by stage 40 to call
#     gametocytes cluster by cluster
# Parameters: config/params.yaml, blood_embedding.
#
# Writes blood.rds (Seurat object of these cells) and blood_cells.tsv (cell,
# cluster, UMAP coordinates).

STAGE <- "35_blood_embedding"
source("R/lib/setup.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

out <- stage_dir(STAGE)
p <- PARAMS$blood_embedding
dims <- seq_len(p$n_pcs_used)

atlas <- readRDS(file.path(stage_dir("30_integrate"), "atlas_integrated.rds"))
eligible <- colnames(atlas)[atlas$tissue == "Host blood" & atlas$nCount_RNA >= PARAMS$annotation$idc_min_umi]
blood <- DietSeurat(subset(atlas, cells = eligible), layers = "counts", assays = "RNA")
rm(atlas)
message(sprintf("%d blood cells: %s", ncol(blood),
  paste(names(table(blood$study_label)), table(blood$study_label), collapse = ", ")))

blood[["RNA"]] <- split(blood[["RNA"]], f = blood$study_label)
blood <- NormalizeData(blood)
blood <- FindVariableFeatures(blood, nfeatures = p$n_variable_genes)
blood <- ScaleData(blood)
blood <- RunPCA(blood, npcs = 50, seed.use = p$seed)
blood <- IntegrateLayers(blood, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
blood <- JoinLayers(blood)
blood <- RunUMAP(blood, reduction = "harmony", dims = dims, n.components = 2, seed.use = p$seed, reduction.name = "umap_blood")
blood <- FindNeighbors(blood, reduction = "harmony", dims = dims)
blood <- FindClusters(blood, resolution = p$cluster_resolution, random.seed = 0)
blood$blood_cluster <- blood$seurat_clusters
message(sprintf("%d clusters at resolution %s", length(levels(blood$blood_cluster)), p$cluster_resolution))

umap <- Embeddings(blood, "umap_blood")
write_table(data.frame(cell = colnames(blood), blood_cluster = blood$blood_cluster,
  umap_b_1 = umap[, 1], umap_b_2 = umap[, 2]), file.path(out, "blood_cells.tsv"))
save_object(blood, file.path(out, "blood.rds"))
save_session_info(out)
