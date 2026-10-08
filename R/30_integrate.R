# Stage 30: merge the studies, normalise, integrate, embed and cluster.
#
#  1. Remove the rRNA genes from every study (only the Mancio-Silva2022
#     matrix still has them; the STARsolo import removes them earlier).
#  2. Define the genes "measured" in each study: the genes it was counted
#     against, whether or not they are detected:
#       - STARsolo studies: the STARsolo gene list (stage 00)
#       - Ruberto2022_1: the genes in the authors' matrix
#       - Mancio-Silva2022: the genes of the authors' matrix
#     Variable genes are chosen among the genes measured in all six studies,
#     so that a gene missing from one study's annotation cannot look variable.
#  3. Merge the studies (one count layer per study), log-normalise, choose
#     the variable genes, scale, PCA, UMAP (unintegrated).
#  4. Harmony through IntegrateLayers(): with one layer per study, the study
#     is the batch variable. Then UMAP, t-SNE, nearest-neighbour graph and
#     Louvain clusters at several resolutions on the Harmony embedding.
#
# Reproducibility: PCA, UMAP and t-SNE use fixed seeds; Harmony uses the
# random number state left by the unintegrated UMAP, so the order of the
# steps below must not change.
#
# Writes atlas_integrated.rds and the tables rrna_genes_removed.tsv and
# measured_genes.tsv.

STAGE <- "30_integrate"
source("R/lib/setup.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

out <- stage_dir(STAGE)
p <- PARAMS$integration
dims <- seq_len(p$n_pcs_used)
check_same_set(p$merge_order, STUDIES, "studies in integration$merge_order and metadata/studies.tsv")

genes <- utils::read.delim(file.path(stage_dir("00_reference"), "genes.tsv"), stringsAsFactors = FALSE)
rrna <- genes$seurat_id[genes$is_rrna]
starsolo_genes <- gsub("_", "-", utils::read.delim(file.path(stage_dir("00_reference"), "starsolo_genes.tsv"))$gene_id)

objects <- lapply(setNames(p$merge_order, p$merge_order), function(study) {
  readRDS(file.path(stage_dir("20_cell_metadata"), paste0(study, ".rds")))
})

# --- 1. rRNA genes ----------------------------------------------------------

rrna_removed <- do.call(rbind, lapply(names(objects), function(study) {
  counts <- LayerData(objects[[study]], assay = "RNA", layer = "counts")
  removed <- intersect(rownames(counts), rrna)
  data.frame(study_label = study, genes_before = nrow(counts), rrna_genes_removed = length(removed),
    rrna_umi_removed = sum(counts[removed, , drop = FALSE]), total_umi_before = sum(counts),
    genes_after = nrow(counts) - length(removed))
}))
for (study in names(objects)) {
  object <- objects[[study]]
  keep <- setdiff(rownames(object), rrna)
  if (length(keep) < nrow(object)) object <- subset(object, features = keep)
  counts <- LayerData(object, assay = "RNA", layer = "counts")
  object$nCount_RNA <- unname(Matrix::colSums(counts))
  object$nFeature_RNA <- unname(Matrix::colSums(counts > 0))
  objects[[study]] <- object
}
print(rrna_removed, row.names = FALSE)
write_table(rrna_removed, file.path(out, "rrna_genes_removed.tsv"))

# --- 2. measured genes ------------------------------------------------------

starsolo_genes <- setdiff(starsolo_genes, rrna)
ruberto_genes <- read_meta("ruberto2022_1_genes.tsv")
measured <- lapply(setNames(names(objects), names(objects)), function(study) {
  present <- rownames(objects[[study]])
  if (identical(PARAMS$import[[study]]$source, "starsolo")) {
    check(all(present %in% starsolo_genes), study, " has genes outside the STARsolo gene list")
    starsolo_genes
  } else if (study == "Ruberto2022_1") {
    intersect(present, ruberto_genes$gene_id[ruberto_genes$in_author_matrix == "TRUE"])
  } else {
    present
  }
})
all_genes <- Reduce(union, measured)
shared_genes <- Reduce(intersect, measured)
message(sprintf("genes: %d in total, %d measured in all studies; measured per study: %s",
  length(all_genes), length(shared_genes), paste(names(measured), lengths(measured), sep = " ", collapse = ", ")))

measured_table <- data.frame(gene_id = all_genes)
for (study in names(measured)) measured_table[[study]] <- all_genes %in% measured[[study]]
measured_table$in_all_studies <- all_genes %in% shared_genes
write_table(measured_table, file.path(out, "measured_genes.tsv"))

# --- 3. merge, normalise, variable genes, PCA, UMAP -------------------------

atlas <- Reduce(function(x, y) merge(x, y), objects)
rm(objects)
check_unique(colnames(atlas), "cell names in the merged studies")

# IntegrateLayers() treats each count layer as one batch: check that there
# is exactly one layer per study
layers <- Layers(atlas, assay = "RNA", search = "counts")
layer_study <- vapply(layers, function(layer) {
  studies <- unique(atlas$study_label[colnames(LayerData(atlas, layer = layer, assay = "RNA"))])
  check(length(studies) == 1, "count layer ", layer, " has cells of ", length(studies), " studies")
  studies
}, "")
check_same_set(layer_study, STUDIES, "studies of the count layers")
check_unique(layer_study, "studies of the count layers")
message("one count layer per study: ", paste(layers, layer_study, sep = " = ", collapse = ", "))

atlas <- NormalizeData(atlas)

shared <- subset(atlas, features = shared_genes)
shared <- FindVariableFeatures(shared, selection.method = p$variable_gene_method, nfeatures = p$n_variable_genes)
variable_genes <- VariableFeatures(shared)
rm(shared)
check(all(variable_genes %in% shared_genes), "a variable gene is not measured in all studies")
message(sprintf("%d variable genes chosen among the %d shared genes", length(variable_genes), length(shared_genes)))

VariableFeatures(atlas) <- variable_genes
atlas <- ScaleData(atlas, features = variable_genes)
atlas <- RunPCA(atlas, features = variable_genes, seed.use = p$pca_seed)
atlas <- RunUMAP(atlas, reduction = "pca", dims = dims, n.components = p$umap_components,
  n.threads = p$threads, seed.use = p$umap_seed)
# copied as they are (keeping their keys PC_ and umap_)
atlas@reductions$pca_unintegrated <- atlas@reductions$pca
atlas@reductions$umap_unintegrated <- atlas@reductions$umap

# --- 4. Harmony, embeddings and clusters on the integrated PCA --------------

atlas <- IntegrateLayers(atlas, method = HarmonyIntegration, orig.reduction = "pca",
  new.reduction = "pca_integrated", verbose = FALSE)
atlas <- JoinLayers(atlas)

atlas <- RunUMAP(atlas, reduction = "pca_integrated", dims = dims, n.components = p$umap_components,
  n.threads = p$threads, seed.use = p$umap_seed)
atlas@reductions$umap_integrated <- atlas@reductions$umap

atlas <- FindNeighbors(atlas, reduction = "pca_integrated", dims = dims)
resolutions <- seq(p$cluster_resolutions$from, p$cluster_resolutions$to, by = p$cluster_resolutions$by)
for (resolution in resolutions) {
  atlas <- FindClusters(atlas, resolution = resolution, graph.name = "RNA_snn")
  seurat_column <- paste0("RNA_snn_res.", resolution)
  column <- paste0("cluster_res_", sub(".", "_", formatC(resolution, format = "f", digits = 2), fixed = TRUE))
  atlas@meta.data[[column]] <- atlas@meta.data[[seurat_column]]
  atlas@meta.data[[seurat_column]] <- NULL
}

atlas <- RunTSNE(atlas, reduction = "pca_integrated", dims = dims, seed.use = p$tsne_seed,
  dim.embed = p$tsne_components, reduction.name = "tsne_integrated", reduction.key = "tSNE_")

atlas@reductions$pca <- NULL
atlas@reductions$umap <- NULL
atlas@meta.data$seurat_clusters <- NULL

check(!any(rownames(atlas) %in% rrna), "rRNA genes in the integrated atlas")
message(sprintf("atlas: %d cells, %d genes", ncol(atlas), nrow(atlas)))
save_object(atlas, file.path(out, "atlas_integrated.rds"))
save_session_info(out)
