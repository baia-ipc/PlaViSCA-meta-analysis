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
# against the value under complete mixing, 1 - sum(p_i^2), where p_i is the
# share of study i among all cells of the stage (ratio 1 = complete mixing).
# The random draws continue the same random number stream as the global
# subsample, stage after stage.
#
# Writes to results/integration/: count tables (study x stage, study x broad
# stage, cluster x study, cluster x stage), integration_metrics.tsv,
# per_stage_sharing.tsv, per_stage_pairwise_sharing.tsv and the figure
# integration.pdf/.png (UMAPs before/after integration, per-stage sharing).

STAGE <- "60_assess_integration"
source("R/lib/setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

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
umaps <- list(unintegrated = embedding("umap_unintegrated", 1:2), integrated = embedding("umap_integrated", 1:2))
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
    expected_complete_mixing = as.numeric(shares[pairs$to_study]), row.names = NULL)
}

shared_stages <- names(which(rowSums(table(stage, study) > 0) > 1))
pairwise <- list()
per_stage <- do.call(rbind, lapply(shared_stages, function(s) {
  cells <- which(stage == s)
  n_stage <- length(cells)
  if (n_stage > p$max_cells_per_stage) cells <- sample(cells, p$max_cells_per_stage)
  knns <- lapply(pcs, function(x) nearest_neighbours(x[cells, , drop = FALSE], p$k_neighbours)$nn)
  observed <- vapply(knns, fraction_other, 1, labels = study[cells])
  shares <- prop.table(table(study[stage %in% s]))
  expected <- 1 - sum(shares^2)
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

# --- figure -----------------------------------------------------------------
# PLOS figure format: one image, at most 7.5 in wide, Helvetica, >= 8 pt.
# Colours: Okabe-Ito for studies, Paul Tol "muted" for stages (both
# colour-blind safe, no shared hues).

study_colours <- c("Hazzard2022" = "#E69F00", "Hazzard2024" = "#56B4E9", "Mancio-Silva2022" = "#009E73",
  "Ruberto2022_1" = "#CC79A7", "Ruberto2022_2" = "#D55E00", "Sa2020" = "#0072B2")
stage_order <- c("Sporozoite", "Hypnozoite", "Schizont (Liver stage)", "Ring", "Trophozoite",
  "Schizont (Blood stage)", "Merozoite", "Female gametocyte", "Male gametocyte")
tol_muted <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
no_stage <- "Unassigned (NA)"
embedding_colours <- c("Unintegrated" = "#BBBBBB", "Integrated (Harmony)" = "#004488")

check(all(study %in% names(study_colours)), "no figure colour for a study")
check(all(stage %in% c(NA, stage_order)), "no figure colour for a stage")
present_stages <- intersect(stage_order, unique(stage))
check(length(present_stages) <= length(tol_muted), "more stages than stage colours")
stage_colours <- c(setNames(tol_muted[seq_along(present_stages)], present_stages), setNames("#DDDDDD", no_stage))
stage_shown <- ifelse(is.na(stage), no_stage, stage)

# draw the studies interleaved in proportion to their size, so that no study
# is painted over another
rank_in_study <- ave(seq_len(n_cells), study, FUN = seq_along)
draw_order <- order(rank_in_study / as.numeric(table(study)[study]))

text_theme <- theme(
  text = element_text(family = "Helvetica", size = 8), axis.text = element_text(size = 8),
  axis.title = element_text(size = 8), plot.title = element_text(size = 9, face = "bold", hjust = 0),
  legend.text = element_text(size = 8), legend.title = element_text(size = 8, face = "bold"),
  legend.key.size = unit(10, "pt"), plot.tag = element_text(size = 11, face = "bold"))

umap_panel <- function(which, colour_by, title) {
  df <- data.frame(x = umaps[[which]][draw_order, 1], y = umaps[[which]][draw_order, 2],
    study = factor(study[draw_order], levels = names(study_colours)),
    stage = factor(stage_shown[draw_order], levels = names(stage_colours)))
  colours <- if (colour_by == "study") scale_colour_manual(name = "Study", values = study_colours, drop = TRUE)
    else scale_colour_manual(name = "Stage", values = stage_colours, drop = TRUE)
  ggplot(df, aes(x = x, y = y, colour = .data[[colour_by]])) +
    ggrastr::rasterise(geom_point(size = 0.05, stroke = 0, shape = 16), dpi = 300, dev = "ragg") +
    colours + guides(colour = guide_legend(override.aes = list(size = 2.5))) + coord_equal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 8, base_family = "Helvetica") + text_theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}

sharing <- rbind(
  data.frame(stage = per_stage$stage, embedding = "Unintegrated", ratio = per_stage$ratio_unintegrated),
  data.frame(stage = per_stage$stage, embedding = "Integrated (Harmony)", ratio = per_stage$ratio_integrated))
sharing$stage <- factor(sharing$stage, levels = rev(intersect(stage_order, per_stage$stage)))
sharing$embedding <- factor(sharing$embedding, levels = names(embedding_colours))
arrows <- data.frame(stage = factor(per_stage$stage, levels = levels(sharing$stage)),
  from = per_stage$ratio_unintegrated, to = per_stage$ratio_integrated)

sharing_panel <- ggplot(sharing, aes(x = ratio, y = stage)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_segment(data = arrows, aes(x = from, xend = to, y = stage, yend = stage), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(colour = embedding), size = 2.2) +
  scale_colour_manual(name = "Embedding", values = embedding_colours) +
  scale_x_continuous(limits = c(0, max(1.05, ceiling(max(sharing$ratio) * 20) / 20)), breaks = seq(0, 1.5, 0.25),
    expand = expansion(mult = c(0, 0.01))) +
  labs(title = "Cross-study neighbour sharing within stage",
    x = sprintf("Observed / expected fraction of %d nearest neighbours from another study\n(1 = complete mixing)", p$k_neighbours),
    y = NULL) +
  theme_classic(base_size = 8, base_family = "Helvetica") + text_theme +
  theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))

figure <- (umap_panel("unintegrated", "study", "Unintegrated: study") | umap_panel("unintegrated", "stage", "Unintegrated: stage")) /
  (umap_panel("integrated", "study", "Integrated: study") | umap_panel("integrated", "stage", "Integrated: stage")) /
  free(sharing_panel) +
  plot_layout(heights = c(1, 1, 0.75), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "right")

pdf(file.path(out, "integration.pdf"), width = 7.5, height = 8, family = "Helvetica", useDingbats = FALSE)
print(figure)
invisible(dev.off())
record_output(file.path(out, "integration.pdf"))
ragg::agg_png(file.path(out, "integration.png"), width = 7.5, height = 8, units = "in", res = 300)
print(figure)
invisible(dev.off())
record_output(file.path(out, "integration.png"))
