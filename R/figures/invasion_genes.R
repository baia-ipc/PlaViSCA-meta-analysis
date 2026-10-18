# Fig 4: invasion- and vaccine-related genes across the life cycle, on the
# integrated UMAP. (A) broad stage, (B-D) expression of DBP, CSP and AMA1,
# (E) fraction of cells detecting each gene per stage, one point per study.

# minimum cells of a study in a stage for a point in panel E
INVASION_MIN_CELLS <- 20

broad_stage_panel <- function(meta, xy) {
  df <- meta
  df$x <- xy[, 1]
  df$y <- xy[, 2]
  df$broad_stage <- factor(df$broad_stage, levels = names(BROAD_STAGE_COLOURS))
  embedding_plot(df[interleaved_order(df$broad_stage), ], "broad_stage", "Broad stage") +
    discrete_scale(BROAD_STAGE_COLOURS, c("Blood stage" = "Blood stages", "Liver stage" = "Liver stages",
      "Sporozoite stage" = "Sporozoites"), "Stage")
}

# cells with expression drawn last (on top), the others in grey below
gene_panel <- function(xy, values, gene) {
  df <- data.frame(x = xy[, 1], y = xy[, 2], value = ifelse(values > 0, values, NA))
  df <- df[order(!is.na(df$value), df$value), ]
  embedding_plot(df, "value", gene) +
    continuous_scale(df$value, "Expression", limits = c(0, stats::quantile(df$value, 0.99, na.rm = TRUE))) +
    theme(plot.title = element_text(face = "bold.italic"))
}

detection_panel <- function(meta, expression) {
  df <- do.call(rbind, lapply(colnames(expression), function(gene) {
    meta |>
      dplyr::mutate(detected = expression[, gene] > 0) |>
      dplyr::filter(!is.na(stage)) |>
      dplyr::group_by(stage, study_label) |>
      dplyr::summarise(cells = dplyr::n(), detected = mean(detected), .groups = "drop") |>
      dplyr::mutate(gene = gene)
  }))
  df <- df[df$cells >= INVASION_MIN_CELLS, ]
  df$gene <- factor(df$gene, levels = colnames(expression))
  df$stage <- factor(df$stage, levels = names(STAGE_COLOURS))
  df$study_label <- factor(df$study_label, levels = names(STUDY_COLOURS))
  ggplot(df, aes(stage, detected, colour = study_label)) +
    geom_point(position = position_dodge(width = 0.6), size = 1.3) +
    facet_grid(gene ~ .) +
    scale_colour_manual(name = "Study", values = STUDY_COLOURS, labels = STUDY_NAMES) +
    scale_x_discrete(labels = STAGE_NAMES) +
    scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(title = "Cells detecting the gene, per stage and study", x = NULL, y = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.2))) +
    theme_figure() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0, face = "bold.italic"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(fill = NA, colour = "grey80", linewidth = 0.3),
      plot.margin = margin(5.5, 5.5, 5.5, 30))
}

figure_invasion_genes <- function(meta, xy, expression) {
  genes <- colnames(expression)
  (broad_stage_panel(meta, xy) | gene_panel(xy, expression[, genes[1]], genes[1])) /
    (gene_panel(xy, expression[, genes[2]], genes[2]) | gene_panel(xy, expression[, genes[3]], genes[3])) /
    detection_panel(meta, expression) +
    plot_layout(heights = c(1, 1, 1.1)) + tag_panels()
}
