# Blood-stage cells in their own embedding (stage 35) and the gametocyte
# calls of stage 40.
#   Fig 2  (A) blood stage and gametocyte state, (B) UMIs, (C, D) female and
#          male marker scores, (E) curated marker panel per state,
#          (F) our calls in the cells labelled by the source studies
#   S3 Fig blood-stage clusters and their stage
#   S4 Fig how clusters are called: marker enrichment and relative UMIs

# blood: one row per blood-stage cell with x, y (blood-only UMAP) and the
# per-cell metadata; drawn in random order
blood_points <- function(blood) {
  set.seed(1)
  blood[sample(nrow(blood)), ]
}

blood_state_panel <- function(blood, title = "Blood stage", point_size = 0.15) {
  blood$state <- factor(blood$state, levels = names(STATE_COLOURS))
  embedding_plot(blood_points(blood), "state", title, point_size) + discrete_scale(STATE_COLOURS, STATE_NAMES)
}

# the colour bar sits right of the plot, clear of the cells
blood_value_panel <- function(blood, column, title, legend) {
  embedding_plot(blood_points(blood), column, title, point_size = 0.15) + continuous_scale(blood[[column]], legend) +
    theme(legend.position = "right")
}

marker_panel_plot <- function(marker_panel) {
  states <- c("Asexual", "Developing gametocyte", "Early female", "Late female", "Early male", "Late male")
  marker_panel$state <- factor(marker_panel$state, levels = states,
    labels = c("Asexual", "Developing", "Early female", "Late female", "Early male", "Late male"))
  gene_label <- paste0(marker_panel$gene, " (", sub("^PVP01-", "", marker_panel$gene_id), ")")
  marker_panel$label <- factor(gene_label, levels = rev(unique(gene_label)))
  marker_panel$marker_set <- factor(marker_panel$marker_set, levels = c("female", "male"), labels = c("Female", "Male"))
  ggplot(marker_panel, aes(state, label, size = detected, colour = mean)) +
    geom_point() +
    scale_colour_gradientn(name = "Mean\nexpression", colours = EXPRESSION_RAMP,
      guide = guide_colourbar(barwidth = unit(5, "pt"), barheight = unit(40, "pt"))) +
    scale_size_area(name = "Fraction\ndetected", max_size = 3.5, breaks = c(0.25, 0.5, 0.75)) +
    facet_grid(marker_set ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Curated sexual-stage genes", x = NULL, y = NULL) +
    theme_figure(7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = -90),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.25), legend.box = "vertical",
      # the title starts at the gene labels, a little right of the panel tag:
      # from the narrow panel it would run into the next panel
      plot.title.position = "plot", plot.title = element_text(margin = margin(l = 14, b = 5.5)))
}

# author_calls: study_label, author_label, cluster_gametocyte_call, cells
author_concordance_plot <- function(author_calls) {
  group <- with(author_calls, dplyr::case_when(
    study_label == "Sa2020" & author_label %in% "Female gametocyte" ~ "Sá 2020: female",
    study_label == "Sa2020" & author_label %in% "Male gametocyte" ~ "Sá 2020: male",
    study_label == "Sa2020" ~ "Sá 2020: other",
    study_label == "Hazzard2024" & author_label == "Sexual" ~ "Hazzard 2024: sexual",
    study_label == "Hazzard2024" ~ "Hazzard 2024: asexual"))
  check(!anyNA(group), "author label without a group in the concordance figure")
  call_colours <- c("Asexual" = NO_VALUE_COLOUR, STAGE_COLOURS[c("Developing gametocyte", "Female gametocyte", "Male gametocyte")])
  call_names <- c("Asexual" = "Asexual stage", "Developing gametocyte" = "Developing gametocyte",
    "Female gametocyte" = "Female gametocyte", "Male gametocyte" = "Male gametocyte")
  df <- data.frame(group, call = ifelse(author_calls$cluster_gametocyte_call == "none", "Asexual", author_calls$cluster_gametocyte_call),
    cells = author_calls$cells) |>
    dplyr::group_by(group, call) |> dplyr::summarise(cells = sum(cells), .groups = "drop")
  totals <- df |> dplyr::group_by(group) |> dplyr::summarise(cells = sum(cells))
  order <- c("Sá 2020: female", "Sá 2020: male", "Sá 2020: other", "Hazzard 2024: sexual", "Hazzard 2024: asexual")
  df$group <- factor(df$group, levels = rev(order))
  df$call <- factor(df$call, levels = names(call_colours))
  totals$group <- factor(totals$group, levels = rev(order))
  ggplot(df, aes(cells, group, fill = call)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7, colour = "white", linewidth = 0.2) +
    geom_text(data = totals, aes(x = 1.02, y = group, label = paste0("n = ", formatC(cells, format = "d", big.mark = ","))),
      inherit.aes = FALSE, hjust = 0, size = 6.5 / .pt, family = "Helvetica") +
    scale_fill_manual(name = "Call in this resource", values = call_colours, labels = call_names) +
    scale_x_continuous(labels = scales::label_percent(), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.45))) +
    coord_cartesian(clip = "off") +
    labs(title = "Source-study labels", x = "Share of cells", y = NULL) +
    guides(fill = guide_legend(ncol = 2, title.position = "top")) +
    theme_figure(7) + theme(legend.position = "bottom")
}

figure_blood_stages <- function(blood, marker_panel, author_calls) {
  blood$log10_umi <- log10(blood$nCount_RNA)
  ((blood_state_panel(blood) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.2)))) |
    blood_value_panel(blood, "log10_umi", "UMIs per cell", "log10")) /
    (blood_value_panel(blood, "female_gam_score", "Female marker score", "Score") |
      blood_value_panel(blood, "male_gam_score", "Male marker score", "Score")) /
    (free(marker_panel_plot(marker_panel)) + author_concordance_plot(author_calls) + plot_layout(widths = c(1.4, 1))) +
    plot_layout(heights = c(1.2, 1.2, 1.1)) + tag_panels()
}

# --- S3 Fig -----------------------------------------------------------------

figure_blood_clusters <- function(blood) {
  blood$blood_cluster <- factor(blood$blood_cluster, levels = sort(as.integer(unique(blood$blood_cluster))))
  labels <- blood |>
    dplyr::group_by(blood_cluster) |>
    dplyr::summarise(x = stats::median(x), y = stats::median(y), .groups = "drop")
  label_layer <- geom_label(data = labels, aes(x, y, label = blood_cluster), inherit.aes = FALSE, size = 6 / .pt,
    label.padding = unit(0.06, "lines"), label.size = 0.15, alpha = 0.85, family = "Helvetica")
  n_clusters <- nlevels(blood$blood_cluster)
  set.seed(1)
  cluster_colours <- setNames(grDevices::hcl.colors(n_clusters, "Dark 3")[sample.int(n_clusters, n_clusters)], levels(blood$blood_cluster))
  clusters <- embedding_plot(blood_points(blood), "blood_cluster", sprintf("Blood-stage clusters (resolution %s)",
    PARAMS$blood_embedding$cluster_resolution), point_size = 0.15) +
    scale_colour_manual(values = cluster_colours, guide = "none") + label_layer
  (clusters | blood_state_panel(blood, point_size = 0.15) + label_layer) + tag_panels()
}

# --- S4 Fig -----------------------------------------------------------------

# cluster_scores: one row per blood-stage cluster with female, male,
# male_mature (enrichment, percentage points), relative_umi, call, state.
# p: annotation parameters (thresholds of the rule).
figure_gametocyte_calls <- function(cluster_scores, p) {
  cluster_scores$state[is.na(cluster_scores$state) | cluster_scores$state == ""] <- "Asexual"
  cluster_scores$state <- factor(cluster_scores$state,
    levels = c("Asexual", "Developing gametocyte", "Early female", "Late female", "Early male", "Late male"))
  state_colours <- c("Asexual" = "grey65", STATE_COLOURS[c("Developing gametocyte", "Early female", "Late female", "Early male", "Late male")])
  state_names <- c("Asexual" = "Asexual", STATE_NAMES[c("Developing gametocyte", "Early female", "Late female", "Early male", "Late male")])
  # only the gametocyte clusters are labelled (all clusters: S3 Fig)
  points <- list(geom_point(size = 1.6),
    geom_text(data = function(d) d[d$state != "Asexual", ], aes(label = blood_cluster), size = 6 / .pt, nudge_y = 3,
      colour = "grey25", family = "Helvetica"),
    scale_colour_manual(name = "Cluster call", values = state_colours, labels = state_names, drop = FALSE))
  threshold <- function(...) list(linetype = "dashed", colour = "grey45", linewidth = 0.3, ...)

  branch <- ggplot(cluster_scores, aes(female, male, colour = state)) +
    do.call(geom_vline, threshold(xintercept = p$gametocyte_min_enrichment)) +
    do.call(geom_hline, threshold(yintercept = p$gametocyte_min_enrichment)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey80", linewidth = 0.3) +
    points + labs(title = "Branch", x = "Female markers (pp)", y = "Male markers (pp)") + theme_figure()

  female <- cluster_scores[cluster_scores$state %in% c("Developing gametocyte", "Early female", "Late female"), ]
  female_state <- ggplot(female, aes(relative_umi, female, colour = state)) +
    do.call(geom_vline, threshold(xintercept = c(p$developing_max_rel_umi, p$late_female_min_rel_umi))) +
    points + scale_x_log10() + guides(colour = "none") +
    labs(title = "Female branch", x = "Relative median UMIs (log scale)", y = "Female markers (pp)") +
    theme_figure()

  male <- cluster_scores[cluster_scores$state %in% c("Early male", "Late male"), ]
  male_state <- ggplot(male, aes(male, male_mature, colour = state)) +
    do.call(geom_hline, threshold(yintercept = p$late_male_min_enrichment)) +
    points + expand_limits(y = 0, x = 0) + guides(colour = "none") +
    labs(title = "Male branch", x = "Male markers (pp)", y = "Known male markers (pp)") +
    theme_figure()

  (branch | female_state | male_state) + plot_layout(guides = "collect", widths = c(1.2, 1, 0.8)) + tag_panels()
}
