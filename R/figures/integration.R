# Integration of the studies (stage 30) and its assessment (stage 60).
#   Fig 3  (A-D) UMAP before and after Harmony, by study and by stage,
#          (E) global mixing and separation metrics, (F) cross-study
#          neighbour sharing within each stage
#   S5 Fig technical variables on the embeddings: platform and UMIs

EMBEDDING_COLOURS <- c("Unintegrated" = "#BBBBBB", "Integrated (Harmony)" = "#004488")

# meta: per-cell metadata; xy: 2D embedding (rows = cells of meta);
# cells of each study interleaved, so that no study hides another
atlas_points <- function(meta, xy) {
  df <- meta
  df$x <- xy[, 1]
  df$y <- xy[, 2]
  df[interleaved_order(df$study_label), ]
}

study_panel <- function(meta, xy, title) {
  df <- atlas_points(meta, xy)
  df$study_label <- factor(df$study_label, levels = names(STUDY_COLOURS))
  embedding_plot(df, "study_label", title) + discrete_scale(STUDY_COLOURS, STUDY_NAMES, "Study")
}

stage_panel <- function(meta, xy, title) {
  df <- atlas_points(meta, xy)
  df$stage <- factor(df$stage, levels = names(STAGE_COLOURS))
  embedding_plot(df, "stage", title) + discrete_scale(STAGE_COLOURS, STAGE_NAMES, "Stage")
}

# dumbbell: one row per item, a grey point before and a dark point after
before_after_plot <- function(df, title, x_label) {
  long <- rbind(data.frame(item = df$item, embedding = "Unintegrated", value = df$before),
    data.frame(item = df$item, embedding = "Integrated (Harmony)", value = df$after))
  long$embedding <- factor(long$embedding, levels = names(EMBEDDING_COLOURS))
  long$item <- factor(long$item, levels = rev(df$item))
  df$item <- factor(df$item, levels = rev(df$item))
  ggplot(long, aes(value, item)) +
    geom_segment(data = df, aes(x = before, xend = after, y = item, yend = item), colour = "grey75", linewidth = 0.5,
      arrow = arrow(length = unit(3, "pt"), type = "closed")) +
    geom_point(aes(colour = embedding), size = 1.8) +
    scale_colour_manual(name = "Embedding", values = EMBEDDING_COLOURS) +
    labs(title = title, x = x_label, y = NULL) +
    theme_figure() + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
}

metrics_panel <- function(metrics) {
  names <- c(
    knn_other_study = "Other-study neighbors",
    silhouette_study = "Silhouette, study",
    knn_other_stage = "Other-stage neighbors",
    silhouette_stage = "Silhouette, stage")
  check_same_set(metrics$metric, names(names), "metrics of integration_metrics.tsv and of the figure")
  metrics <- metrics[match(names(names), metrics$metric), ]
  before_after_plot(data.frame(item = unname(names), before = metrics$unintegrated, after = metrics$integrated),
    "Mixing and stage separation", sprintf("Value (%s cells, k = %s)",
      formatC(metrics$n_cells[1], format = "d", big.mark = ","), metrics$k[1])) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    scale_x_continuous(breaks = c(-0.2, 0, 0.2))
}

sharing_panel <- function(sharing) {
  sharing <- sharing[order(match(sharing$stage, names(STAGE_COLOURS))), ]
  before_after_plot(data.frame(item = STAGE_NAMES[sharing$stage], before = sharing$ratio_unintegrated,
    after = sharing$ratio_integrated), "Cross-study sharing within stage",
    "Observed / expected neighbors from\nanother study (1 = complete mixing)") +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    scale_x_continuous(limits = c(0, max(1.05, sharing$ratio_integrated)), breaks = c(0, 0.5, 1))
}

figure_integration <- function(meta, umap, metrics, sharing) {
  (study_panel(meta, umap$unintegrated, "Unintegrated: study") | stage_panel(meta, umap$unintegrated, "Unintegrated: stage")) /
    (study_panel(meta, umap$integrated, "Integrated: study") | stage_panel(meta, umap$integrated, "Integrated: stage")) /
    (free(metrics_panel(metrics)) | sharing_panel(sharing)) +
    plot_layout(heights = c(1, 1, 0.8), guides = "collect") + tag_panels() & theme(legend.position = "right")
}

# --- S5 Fig -----------------------------------------------------------------

PLATFORM_COLOURS <- c("10x (v3)" = "#BBBBBB", "10x (v2)" = "#0072B2", "Seq-Well" = "#D55E00")

figure_technical <- function(meta, umap, studies) {
  meta$platform <- studies$platform[match(meta$study_label, studies$study_label)]
  check(all(meta$platform %in% names(PLATFORM_COLOURS)), "a platform has no figure colour")
  meta$log10_umi <- log10(meta$nCount_RNA)
  platform_panel <- function(xy, title) {
    df <- meta
    df$x <- xy[, 1]
    df$y <- xy[, 2]
    df$platform <- factor(df$platform, levels = names(PLATFORM_COLOURS))
    # the large 10x v3 group below the two small ones
    embedding_plot(df[order(df$platform), ], "platform", title) + discrete_scale(PLATFORM_COLOURS, title = "Platform")
  }
  umi_panel <- function(xy, title) {
    df <- meta
    df$x <- xy[, 1]
    df$y <- xy[, 2]
    set.seed(1)
    embedding_plot(df[sample(nrow(df)), ], "log10_umi", title) + continuous_scale(df$log10_umi, "UMIs\n(log10)")
  }
  (platform_panel(umap$unintegrated, "Unintegrated: platform") | platform_panel(umap$integrated, "Integrated: platform") |
    umi_panel(umap$integrated, "Integrated: UMIs per cell")) +
    plot_layout(guides = "collect") + tag_panels()
}
