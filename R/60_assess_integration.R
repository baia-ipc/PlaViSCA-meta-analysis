# Stage 60: how well the integration mixes studies and keeps stages apart.
#
# Compares the PCA before integration (pca_unintegrated) with the Harmony
# PCA (pca_integrated), on PCs 1..n_pcs, with Euclidean distances and the k
# nearest neighbours of each cell (itself excluded). Parameters:
# config/params.yaml, assessment.
#
# Global metrics, on one random subsample of n_cells cells (the same cells
# for both embeddings):
#   knn_other_study   mean fraction of a cell's neighbours from another study
#                     (higher = studies mix better)
#   knn_other_stage   same for the stage; cells without a stage are left out
#                     (lower = stages stay apart)
#   silhouette_study  mean silhouette width by study (lower = better mixing)
#   silhouette_stage  mean silhouette width by stage (higher = stages apart)
#
# Per-stage sharing, for each stage present in more than one study, on up to
# max_cells_per_stage cells of the stage: fraction of a cell's neighbours
# among the cells of the same stage that come from another study (observed),
# against the value under complete mixing (ratio 1 = complete mixing). Under
# complete mixing the neighbours of a cell of study i are a random draw of
# the other n - 1 cells, so the expected value is sum_i p_i (n - n_i) / (n - 1),
# where n_i is the number of cells of study i among the n cells used and
# p_i = n_i / n.
# The random draws continue the same random number stream as the global
# subsample, stage after stage.
#
# Writes to results/integration/: count tables (study x stage, study x broad
# stage, cluster x study, cluster x stage), integration_metrics.tsv,
# per_stage_sharing.tsv and per_stage_pairwise_sharing.tsv. The figure is
# made by stage 80.

STAGE <- "60_assess_integration"
source("R/lib/setup.R")

out <- results_dir("integration")
p <- PARAMS$assessment
dims <- seq_len(p$n_pcs)

atlas <- readRDS(file.path(stage_dir("40_annotate_stages"), "atlas.rds"))
meta <- atlas@meta.data
embedding <- function(name, dims) {
  x <- atlas@reductions[[name]]@cell.embeddings
  check_same_order(rownames(x), rownames(meta), paste("cells of", name, "and of the metadata"))
  x[, dims, drop = FALSE]
}
pcs <- list(unintegrated = embedding("pca_unintegrated", dims), integrated = embedding("pca_integrated", dims))
rm(atlas)

study <- meta$study_label
stage <- meta$stage
n_cells <- nrow(meta)

# --- count tables -----------------------------------------------------------

write_count_table <- function(rows, cols, rows_name, cols_name, file) {
  tab <- addmargins(table(rows, cols, useNA = "ifany"))
  label <- function(x) ifelse(is.na(x), "NA", x)
  df <- data.frame(label(rownames(tab)), unclass(tab), check.names = FALSE, row.names = NULL)
  colnames(df) <- c(paste0(rows_name, "\\", cols_name), label(colnames(tab)))
  write_table(df, file.path(out, file))
}
write_count_table(study, stage, "study_label", "stage", "study_by_stage.tsv")
write_count_table(study, meta$broad_stage, "study_label", "broad_stage", "study_by_broad_stage.tsv")
write_count_table(meta$cluster, study, "cluster", "study_label", "cluster_by_study.tsv")
write_count_table(meta$cluster, stage, "cluster", "stage", "cluster_by_stage.tsv")

# --- metrics ----------------------------------------------------------------

# k nearest neighbours (row indices, self excluded) and the distances
nearest_neighbours <- function(x, k) {
  d <- dist(x)
  dm <- as.matrix(d)
  diag(dm) <- Inf
  list(nn = t(apply(dm, 1, function(r) order(r)[seq_len(k)])), dist = d)
}

# mean fraction of neighbours with a different label; cells and neighbours
# without a label are left out
fraction_other <- function(nn, labels) {
  neighbour_labels <- matrix(labels[nn], nrow = nrow(nn))
  per_cell <- rowMeans(neighbour_labels != labels, na.rm = TRUE)
  mean(per_cell[!is.na(labels)], na.rm = TRUE)
}

mean_silhouette <- function(d, labels) {
  keep <- !is.na(labels)
  sil <- cluster::silhouette(as.integer(factor(labels[keep])), dmatrix = as.matrix(d)[keep, keep])
  mean(sil[, "sil_width"])
}

set.seed(p$seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
sample_cells <- sample(n_cells, p$n_cells)

global <- lapply(pcs, function(x) {
  knn <- nearest_neighbours(x[sample_cells, , drop = FALSE], p$k_neighbours)
  c(knn_other_study = fraction_other(knn$nn, study[sample_cells]),
    knn_other_stage = fraction_other(knn$nn, stage[sample_cells]),
    silhouette_study = mean_silhouette(knn$dist, study[sample_cells]),
    silhouette_stage = mean_silhouette(knn$dist, stage[sample_cells]))
})
metrics <- data.frame(
  metric = names(global$unintegrated),
  unintegrated = unname(global$unintegrated),
  integrated = unname(global$integrated),
  better_if = c("higher (study mixing)", "lower (stage separation)", "lower (study mixing)", "higher (stage separation)"),
  n_cells = p$n_cells, k = p$k_neighbours, pcs = sprintf("1:%d", p$n_pcs), seed = p$seed
)
print(metrics)
write_table(metrics, file.path(out, "integration_metrics.tsv"))

# fraction of neighbours of study A cells that come from study B, for each
# ordered pair of studies of a stage (no random draws)
pairwise_sharing <- function(stage_name, knns, labels, shares) {
  pairs <- expand.grid(from_study = names(shares), to_study = names(shares), stringsAsFactors = FALSE)
  pairs <- pairs[pairs$from_study != pairs$to_study, ]
  fraction <- function(nn, a, b) {
    from <- labels == a
    if (!any(from)) return(NA_real_)
    mean(rowMeans(matrix(labels[nn[from, , drop = FALSE]], nrow = sum(from)) == b))
  }
  data.frame(stage = stage_name, from_study = pairs$from_study, to_study = pairs$to_study,
    n_from_sampled = vapply(pairs$from_study, function(a) sum(labels == a), 1L),
    observed_unintegrated = mapply(function(a, b) fraction(knns$unintegrated, a, b), pairs$from_study, pairs$to_study),
    observed_integrated = mapply(function(a, b) fraction(knns$integrated, a, b), pairs$from_study, pairs$to_study),
    expected_complete_mixing = vapply(pairs$to_study, function(b) sum(labels == b) / (length(labels) - 1), 1),
    row.names = NULL)
}

shared_stages <- names(which(rowSums(table(stage, study) > 0) > 1))
pairwise <- list()
per_stage <- do.call(rbind, lapply(shared_stages, function(s) {
  cells <- which(stage == s)
  n_stage <- length(cells)
  if (n_stage > p$max_cells_per_stage) cells <- sample(cells, p$max_cells_per_stage)
  knns <- lapply(pcs, function(x) nearest_neighbours(x[cells, , drop = FALSE], p$k_neighbours)$nn)
  observed <- vapply(knns, fraction_other, 1, labels = study[cells])
  n_by_study <- table(study[cells])
  shares <- n_by_study / length(cells)
  expected <- sum(shares * (length(cells) - n_by_study) / (length(cells) - 1))
  pairwise[[s]] <<- pairwise_sharing(s, knns, study[cells], shares)
  data.frame(stage = s, n_cells_stage = n_stage, n_cells_sampled = length(cells), n_studies = length(shares),
    observed_unintegrated = observed[["unintegrated"]], observed_integrated = observed[["integrated"]],
    expected_complete_mixing = expected,
    ratio_unintegrated = observed[["unintegrated"]] / expected, ratio_integrated = observed[["integrated"]] / expected)
}))
print(per_stage)
write_table(per_stage, file.path(out, "per_stage_sharing.tsv"))
pairwise <- do.call(rbind, pairwise)
pairwise$ratio_unintegrated <- pairwise$observed_unintegrated / pairwise$expected_complete_mixing
pairwise$ratio_integrated <- pairwise$observed_integrated / pairwise$expected_complete_mixing
write_table(pairwise, file.path(out, "per_stage_pairwise_sharing.tsv"))
