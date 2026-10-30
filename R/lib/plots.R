# Shared style of the manuscript figures (stage 80).
#
# Provides:
#   STUDY_COLOURS, STUDY_NAMES      colour and display name of each study
#   STAGE_COLOURS, STAGE_NAMES      colour and display name of each stage
#   STATE_COLOURS                   blood stages with the gametocyte states
#   BROAD_STAGE_COLOURS             blood, liver and sporozoite stages
#   EXPRESSION_RAMP                 sequential scale for expression, UMIs, scores
#   theme_figure(), theme_embedding()
#   embedding_plot(df, colour, ...) points of a 2D embedding
#   tag_panels()                    A, B, C... tags of a patchwork figure
#   save_figure(plot, name, width, height)  PDF and PNG in results/figures/
#
# PLOS figure format: at most 7.5 x 8.75 in, Helvetica (Arial), text >= 8 pt
# in print. Colours: the stage palette groups the stages by compartment
# (asexual blood cycle from yellow to dark red, gametocytes in pale blue,
# green and blue, liver stages in lilac and purple, sporozoites in dark teal);
# all pairs of STAGE_COLOURS, of STATE_COLOURS and of BROAD_STAGE_COLOURS are
# separable with protanopia and deuteranopia (OKLab delta E >= 8, simulated
# as in Machado et al. 2009). Studies use the Okabe-Ito colours.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

STUDY_COLOURS <- c(
  "Sa2020" = "#0072B2", "Hazzard2022" = "#E69F00", "Hazzard2024" = "#56B4E9",
  "Ruberto2022_1" = "#CC79A7", "Mancio-Silva2022" = "#009E73", "Ruberto2022_2" = "#D55E00")
STUDY_NAMES <- c(
  "Sa2020" = "S\u00e1 2020", "Hazzard2022" = "Hazzard 2022", "Hazzard2024" = "Hazzard 2024",
  "Ruberto2022_1" = "Ruberto 2022 (1)", "Mancio-Silva2022" = "Mancio-Silva 2022", "Ruberto2022_2" = "Ruberto 2022 (2)")

# in the order of metadata/stages.tsv
STAGE_COLOURS <- c(
  "Ring/Early trophozoite" = "#DDDD55", "Trophozoite" = "#EE8844", "Schizont (Blood stage)" = "#BB4400",
  "Merozoite" = "#771100", "Developing gametocyte" = "#AACCDD", "Female gametocyte" = "#118866",
  "Male gametocyte" = "#3377CC", "Hypnozoite" = "#CC88DD", "Schizont (Liver stage)" = "#773388",
  "Sporozoite" = "#114444")
STAGE_NAMES <- c(
  "Ring/Early trophozoite" = "Ring/early trophozoite", "Trophozoite" = "Trophozoite",
  "Schizont (Blood stage)" = "Blood schizont", "Merozoite" = "Late schizont/merozoite",
  "Developing gametocyte" = "Developing gametocyte", "Female gametocyte" = "Female gametocyte",
  "Male gametocyte" = "Male gametocyte", "Hypnozoite" = "Hypnozoite",
  "Schizont (Liver stage)" = "Liver schizont", "Sporozoite" = "Sporozoite")

# blood stages with the early and late gametocyte states: a lighter and a
# darker shade of the colour of each sex
STATE_COLOURS <- c(STAGE_COLOURS[c("Ring/Early trophozoite", "Trophozoite", "Schizont (Blood stage)", "Merozoite",
  "Developing gametocyte")],
  "Early female" = "#00BB99", "Late female" = "#226655", "Early male" = "#3399EE", "Late male" = "#114499")
STATE_NAMES <- c(STAGE_NAMES[1:5], setNames(names(STATE_COLOURS)[6:9], names(STATE_COLOURS)[6:9]))

BROAD_STAGE_COLOURS <- c("Blood stage" = "#BB4400", "Liver stage" = "#773388", "Sporozoite stage" = "#114444")

EXPRESSION_RAMP <- c("#E3E3E3", "#B9C9DE", "#6F92C2", "#2F5DA0", "#0B2A63")
NO_VALUE_COLOUR <- "#D9D9D9"

# --- themes -----------------------------------------------------------------

theme_figure <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(size = base_size, colour = "grey20"),
      axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size), legend.key.size = unit(9, "pt"),
      strip.background = element_blank(), strip.text = element_text(size = base_size, face = "bold"),
      plot.tag = element_text(size = 11, face = "bold"))
}

theme_embedding <- function() {
  theme_figure() + theme(axis.text = element_blank(), axis.ticks = element_blank())
}

# --- building blocks --------------------------------------------------------

# df: columns x, y and the colour column, already in drawing order (last on
# top). A discrete colour gets a legend with large keys.
embedding_plot <- function(df, colour, title = NULL, point_size = 0.05, axis_name = "UMAP") {
  ggplot(df, aes(x, y, colour = .data[[colour]])) +
    ggrastr::rasterise(geom_point(size = point_size, stroke = 0, shape = 16), dpi = 400, dev = "ragg") +
    coord_equal() + labs(title = title, x = paste(axis_name, 1), y = paste(axis_name, 2)) +
    theme_embedding()
}

discrete_scale <- function(colours, labels = setNames(names(colours), names(colours)), title = NULL, ...) {
  list(scale_colour_manual(name = title, values = colours, labels = labels, drop = TRUE, na.value = NO_VALUE_COLOUR, ...),
    guides(colour = guide_legend(override.aes = list(size = 2.2))))
}

# continuous values squished to their 1-99% quantiles, so that a few extreme
# cells do not flatten the scale
continuous_scale <- function(values, title = NULL, limits = stats::quantile(values, c(0.01, 0.99), na.rm = TRUE)) {
  list(scale_colour_gradientn(name = title, colours = EXPRESSION_RAMP, limits = limits, oob = scales::squish,
    na.value = NO_VALUE_COLOUR),
    guides(colour = guide_colourbar(barwidth = unit(5, "pt"), barheight = unit(45, "pt"))))
}

# order so that the cells of each group are interleaved in proportion to the
# group sizes: no group is painted over another
interleaved_order <- function(group) {
  rank_in_group <- ave(seq_along(group), group, FUN = seq_along)
  order(rank_in_group / as.numeric(table(group)[group]))
}

# the tag style is set by theme_figure()
tag_panels <- function() {
  plot_annotation(tag_levels = "A")
}

FIGURES_DIR <- NULL

# PDF and PNG named `name` in FIGURES_DIR; the TIFF for the PLOS submission
# (RGB, LZW, 600 dpi) named by its label in the manuscript (`label`, e.g.
# "Fig1", "S1_Fig") in FIGURES_DIR/plos/
save_figure <- function(plot, name, label, width, height) {
  check(width <= 7.5 && height <= 8.75, "figure ", name, " larger than the PLOS maximum of 7.5 x 8.75 in")
  pdf_path <- file.path(FIGURES_DIR, paste0(name, ".pdf"))
  png_path <- file.path(FIGURES_DIR, paste0(name, ".png"))
  tiff_path <- file.path(FIGURES_DIR, "plos", paste0(label, ".tif"))
  dir.create(dirname(tiff_path), showWarnings = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(tiff_path, plot, width = width, height = height, dpi = 600, device = ragg::agg_tiff,
         bg = "white", compression = "lzw")
  check(file.size(tiff_path) <= 10e6, "figure ", label, " larger than the PLOS maximum of 10 MB")
  for (path in c(pdf_path, png_path, tiff_path)) record_output(path)
  invisible(pdf_path)
}
