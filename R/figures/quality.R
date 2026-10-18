# Cell quality and confidence of the stage call.
#   S1 Fig genes and UMIs per cell, per study
#   S2 Fig reference-based blood-stage call (SingleR): margin of the best
#          time point over the next (delta) and share of labels that
#          SingleR's pruning would discard

QC_STUDY_ORDER <- c("Sa2020", "Hazzard2022", "Hazzard2024", "Ruberto2022_1", "Mancio-Silva2022", "Ruberto2022_2")

figure_qc <- function(meta) {
  df <- rbind(data.frame(study = meta$study_label, measure = "Genes per cell", value = meta$nFeature_RNA),
    data.frame(study = meta$study_label, measure = "UMIs per cell", value = meta$nCount_RNA))
  df$study <- factor(df$study, levels = QC_STUDY_ORDER)
  ggplot(df, aes(study, value, fill = study)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.85) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_fill_manual(values = STUDY_COLOURS, guide = "none") +
    scale_x_discrete(labels = STUDY_NAMES) +
    facet_wrap(~measure, scales = "free_y") +
    labs(x = NULL, y = "Log scale") +
    theme_figure() + theme(axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
}

figure_stage_call <- function(meta) {
  called <- meta[!is.na(meta$idc_timepoint), ]
  timepoints <- sort(unique(as.numeric(called$idc_timepoint)))
  called$timepoint <- factor(as.numeric(called$idc_timepoint), levels = timepoints, labels = paste(timepoints, "h"))
  timepoint_stage <- read_meta("idc_timepoints.tsv")
  stage_of <- setNames(timepoint_stage$stage, paste(timepoint_stage$timepoint_h, "h"))
  check(all(levels(called$timepoint) %in% names(stage_of)), "time point without a stage in idc_timepoints.tsv")
  fill <- STAGE_COLOURS[stage_of[levels(called$timepoint)]]
  names(fill) <- levels(called$timepoint)

  pruned <- called |>
    dplyr::group_by(timepoint) |>
    dplyr::summarise(cells = dplyr::n(), pruned = mean(idc_pruned), .groups = "drop")
  delta <- ggplot(called, aes(timepoint, idc_delta, fill = timepoint)) +
    geom_violin(scale = "width", colour = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    geom_text(data = pruned, aes(timepoint, Inf, label = formatC(cells, format = "d", big.mark = ",")), inherit.aes = FALSE,
      vjust = 1.2, size = 6.5 / .pt, family = "Helvetica", colour = "grey30") +
    scale_fill_manual(values = fill, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = "Margin over the next time point", x = "Most similar reference time point", y = "SingleR delta") +
    theme_figure()

  pruned_plot <- ggplot(pruned, aes(timepoint, pruned, fill = timepoint)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = fill, guide = "none") +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Labels flagged by SingleR pruning", x = "Most similar reference time point", y = "Share of cells") +
    theme_figure() + theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

  (delta | pruned_plot) + plot_layout(widths = c(1.6, 1)) + tag_panels()
}
