# Fig 1: the resource. (A) Processing workflow, (B) cells per study,
# (C) harmonized stage composition of each study.

# studies grouped by compartment, as in the stage-count table
OVERVIEW_STUDY_ORDER <- c("Sa2020", "Hazzard2022", "Hazzard2024", "Ruberto2022_1", "Mancio-Silva2022", "Ruberto2022_2")

workflow_panel <- function(studies, value) {
  n_reads <- sum(grepl("STARsolo", studies$count_source))
  n_matrices <- nrow(studies) - n_reads
  steps <- data.frame(
    title = c("Public data", "Counts and cells", "Integration", "Stage annotation", "PlaViSCA"),
    body = c(
      sprintf("%s studies\n%d from raw reads\n%d from authors'\ncount matrices", value("Studies"), n_reads, n_matrices),
      "STARsolo\nPvP01, PlasmoDB-68\ncells as in the\nsource studies",
      sprintf("log-normalization\n%s variable genes\nPCA, Harmony\nUMAP, t-SNE", value("VariableGenes")),
      "liver, sporozoites:\nsource labels\nblood: SingleR\nand marker genes",
      sprintf("%s cells\n%s genes\nweb application\ndownloadable data", value("CellsTotal"), value("GenesTotal"))))
  width <- 1.75
  gap <- 0.3
  steps$xmin <- (seq_len(nrow(steps)) - 1) * (width + gap)
  steps$xmax <- steps$xmin + width
  steps$x <- (steps$xmin + steps$xmax) / 2
  arrows <- data.frame(x = head(steps$xmax, -1) + 0.04, xend = tail(steps$xmin, -1) - 0.04)
  fill <- c(rep("#F2F2F2", nrow(steps) - 1), "#DCE6F2")
  ggplot(steps) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1), fill = fill, colour = "grey55", linewidth = 0.3) +
    geom_text(aes(x = x, y = 0.85, label = title), size = 7.5 / .pt, fontface = "bold", family = "Helvetica") +
    geom_text(aes(x = x, y = 0.4, label = body), size = 6.5 / .pt, lineheight = 0.95, family = "Helvetica", colour = "grey15") +
    geom_segment(data = arrows, aes(x = x, xend = xend, y = 0.5, yend = 0.5), linewidth = 0.4, colour = "grey35",
      arrow = arrow(length = unit(4, "pt"), type = "closed")) +
    scale_x_continuous(expand = expansion(add = 0.02)) + scale_y_continuous(expand = expansion(add = 0.02)) +
    theme_void(base_family = "Helvetica") + theme(plot.tag = element_text(size = 11, face = "bold", hjust = 0))
}

cells_per_study_panel <- function(meta) {
  n <- as.data.frame(table(study = factor(meta$study_label, levels = rev(OVERVIEW_STUDY_ORDER))), responseName = "cells")
  ggplot(n, aes(cells, study, fill = study)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = formatC(cells, format = "d", big.mark = ",")), hjust = -0.1, size = 7 / .pt, family = "Helvetica") +
    scale_fill_manual(values = STUDY_COLOURS, guide = "none") +
    scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()), breaks = c(1e3, 1e4, 1e5), expand = expansion(mult = c(0, 0.25))) +
    scale_y_discrete(labels = STUDY_NAMES) +
    labs(title = "Cells per study", x = "Cells (log scale)", y = NULL) +
    theme_figure() + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
}

stage_composition_panel <- function(meta) {
  n <- as.data.frame(table(study = factor(meta$study_label, levels = rev(OVERVIEW_STUDY_ORDER)),
    stage = factor(meta$stage, levels = names(STAGE_COLOURS))), responseName = "cells")
  n <- n[n$cells > 0, ]
  ggplot(n, aes(cells, study, fill = stage)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7, colour = "white", linewidth = 0.2) +
    scale_fill_manual(name = "Harmonized stage", values = STAGE_COLOURS, labels = STAGE_NAMES, drop = FALSE) +
    scale_x_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.01))) +
    scale_y_discrete(labels = NULL) +
    labs(title = "Harmonized stage", x = "Share of the study's cells", y = NULL) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_figure() + theme(axis.ticks.y = element_blank())
}

figure_overview <- function(meta, studies, value) {
  free(workflow_panel(studies, value)) /
    (cells_per_study_panel(meta) | stage_composition_panel(meta)) +
    plot_layout(heights = c(1, 2.2)) + tag_panels()
}
