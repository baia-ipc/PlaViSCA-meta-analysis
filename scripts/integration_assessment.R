# integration_assessment.R
#
# Read-only post-integration assessment of the PlaViSca v1 release object
# (plavisca-v1-dec01-dec02-2026-09-21/pv_all_studies.rds). Reproduces the
# metrics reported in the manuscript (Methods "Assessment of integration",
# Results "Stage-resolved integration across studies", Table 3) and Fig 2.
#
# Never modifies the release object or any pipeline artifact: it only reads
# the object (readRDS + slot(), no Seurat attach needed) and writes TSV/PDF/
# PNG outputs to results/integration_assessment/, recording each output in
# build_manifest.tsv via record_build_manifest().
#
# Run (from scripts/):  pixi run integration-assessment
#
# Input override: PLAVISCA_RELEASE_RDS=<path>. The input's SHA256 must match
# EXPECTED_INPUT_SHA256 and the computed metrics must match the published
# values (PUBLISHED_*), otherwise the script fails loudly - the published
# numbers belong to exactly this object and this sampling sequence.

.plavisca_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PREPROCESS_ROOT", unset = NA_character_),
    getwd(),
    dirname(getwd()),
    "/home/sopheap/pvsca_b/pre_process_data"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) {
    stop(
      "Could not locate the pre_process_data project root (looked for scripts/",
      "pipeline_lib.R under $PLAVISCA_PREPROCESS_ROOT, the current working ",
      "directory, its parent, and the historical hard-coded path). Set the ",
      "PLAVISCA_PREPROCESS_ROOT environment variable to the pre_process_data ",
      "directory's absolute path, or run this script from that directory.",
      call. = FALSE
    )
  }
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Every package this script relies on must come from the pixi environment's
# own (locked) library, never from a personal library on R_LIBS_USER.
locked_library <- normalizePath(file.path(R.home(), "library"))
for (p in c("SeuratObject", "Matrix", "cluster", "ggplot2", "patchwork", "ggrastr", "ragg", "openssl")) {
  p_path <- tryCatch(normalizePath(find.package(p)), error = function(e) NA_character_)
  if (is.na(p_path) || !startsWith(p_path, locked_library)) {
    pipeline_fail(sprintf(
      "package '%s' resolved to %s, outside the pixi environment library %s - run via `pixi run integration-assessment` (which clears R_LIBS_USER)",
      p, p_path, locked_library
    ))
  }
}

SCRIPT_PATH <- "scripts/integration_assessment.R"
RELEASE_RDS <- Sys.getenv(
  "PLAVISCA_RELEASE_RDS",
  unset = "/home/baia/prj/plavisca/deployments/data/releases/plavisca-v1-dec01-dec02-2026-09-21/pv_all_studies.rds"
)
PUBLISHED_INPUT_SHA256 <- "32068878fc92952e9be539befcc8ff2f42a6818716a912a55f5c45930cc830b9"
# Any release other than the published 2026-09-21 object must be named
# explicitly by its SHA256 and written to its own output directory; the
# published values below are then reported for comparison, not enforced.
EXPECTED_INPUT_SHA256 <- Sys.getenv("PLAVISCA_RELEASE_SHA256", unset = PUBLISHED_INPUT_SHA256)
IS_PUBLISHED_RELEASE <- identical(EXPECTED_INPUT_SHA256, PUBLISHED_INPUT_SHA256)
OUT_DIR <- Sys.getenv("PLAVISCA_ASSESSMENT_OUT_DIR", unset = "results/integration_assessment")
if (!IS_PUBLISHED_RELEASE && OUT_DIR == "results/integration_assessment") {
  pipeline_fail("a non-published release needs its own PLAVISCA_ASSESSMENT_OUT_DIR so the published outputs are never overwritten")
}

SEED <- 1L
N_SUBSAMPLE <- 8000L
N_PER_STAGE_MAX <- 3000L
K_NEIGHBOURS <- 30L
DIMS <- 1:30
CLUSTER_COLUMN <- "RNA_snn_res.0.05"
EMBEDDINGS <- c(unintegrated = "pca_unintegrated", integrated = "pca_integrated")
UMAPS <- c(unintegrated = "umap_unintegrated", integrated = "umap_integrated")

# Published values (manuscript Table 3 / Results), used as a regression gate.
PUBLISHED_TOLERANCE <- 0.005
PUBLISHED_GLOBAL <- data.frame(
  metric = c("knn_other_study", "knn_other_stage", "silhouette_study", "silhouette_stage"),
  unintegrated = c(0.0441, 0.2996, -0.0737, -0.0546),
  integrated = c(0.222, 0.336, -0.229, -0.141)
)
PUBLISHED_PER_STAGE <- data.frame(
  stage = c(
    "Hypnozoite", "Merozoite", "Ring", "Schizont (Blood stage)",
    "Schizont (Liver stage)", "Sporozoite", "Trophozoite"
  ),
  unintegrated = c(0.038, 0.025, 0.011, 0.021, 0.002, 0.170, 0.008),
  integrated = c(0.165, 0.081, 0.047, 0.160, 0.065, 0.356, 0.112)
)
PUBLISHED_STAGE_TOTALS <- c(
  "Ring" = 8521L, "Trophozoite" = 17646L, "Schizont (Blood stage)" = 57712L,
  "Merozoite" = 6340L, "Hypnozoite" = 1419L, "Schizont (Liver stage)" = 1499L,
  "Sporozoite" = 12805L, "NA" = 7L
)
PUBLISHED_TOTAL_CELLS <- 105949L

# ---------------------------------------------------------------------------
# Input: verify identity, then read slots without attaching Seurat
# ---------------------------------------------------------------------------
if (!file.exists(RELEASE_RDS)) {
  pipeline_fail(sprintf("release object not found: %s", RELEASE_RDS))
}
input_sha256 <- paste(openssl::sha256(file(RELEASE_RDS, raw = TRUE))) # paste() drops the "hash" class
if (!identical(input_sha256, EXPECTED_INPUT_SHA256)) {
  pipeline_fail(sprintf(
    "SHA256 of %s is %s, expected %s - refusing to report metrics for a different object",
    RELEASE_RDS, input_sha256, EXPECTED_INPUT_SHA256
  ))
}

obj <- readRDS(RELEASE_RDS)
md <- slot(obj, "meta.data")
reductions <- slot(obj, "reductions")
n_cells <- nrow(md)
assert_cell_count(n_cells, PUBLISHED_TOTAL_CELLS, label = "release pv_all_studies.rds")

required_columns <- c("study_label", "harmonized_life_cycle_stage", "parasite_broad_stage", CLUSTER_COLUMN)
missing_columns <- setdiff(required_columns, colnames(md))
if (length(missing_columns) > 0) {
  pipeline_fail(sprintf("release object meta.data lacks column(s): %s", paste(missing_columns, collapse = ", ")))
}
missing_reductions <- setdiff(c(EMBEDDINGS, UMAPS), names(reductions))
if (length(missing_reductions) > 0) {
  pipeline_fail(sprintf("release object lacks reduction(s): %s", paste(missing_reductions, collapse = ", ")))
}

get_embedding <- function(name, dims) {
  emb <- slot(reductions[[name]], "cell.embeddings")
  assert_rowname_order_equal(rownames(emb), rownames(md), label_a = name, label_b = "meta.data")
  if (max(dims) > ncol(emb)) {
    pipeline_fail(sprintf("reduction '%s' has %d dims, %d requested", name, ncol(emb), max(dims)))
  }
  emb[, dims, drop = FALSE]
}
pcs <- lapply(EMBEDDINGS, get_embedding, dims = DIMS)
umaps <- lapply(UMAPS, get_embedding, dims = 1:2)

study <- md$study_label
stage <- md$harmonized_life_cycle_stage
broad_stage <- md$parasite_broad_stage
cluster <- md[[CLUSTER_COLUMN]]
assert_no_missing_mandatory(study, label = "study_label")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- function(name) file.path(OUT_DIR, name)

# ---------------------------------------------------------------------------
# Composition tables (with margins)
# ---------------------------------------------------------------------------
write_count_table <- function(row, col, row_name, col_name, path) {
  tab <- addmargins(table(row, col, useNA = "ifany"))
  labels <- function(x) ifelse(is.na(x), "NA", x)
  df <- data.frame(labels(rownames(tab)), unclass(tab)[, , drop = FALSE], check.names = FALSE, row.names = NULL)
  colnames(df) <- c(paste0(row_name, "\\", col_name), labels(colnames(tab)))
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  tab
}

study_by_stage <- write_count_table(study, stage, "study_label", "harmonized_life_cycle_stage", out_path("study_by_stage.tsv"))
invisible(write_count_table(study, broad_stage, "study_label", "parasite_broad_stage", out_path("study_by_broad_stage.tsv")))
invisible(write_count_table(cluster, study, CLUSTER_COLUMN, "study_label", out_path("cluster_by_study.tsv")))
invisible(write_count_table(cluster, stage, CLUSTER_COLUMN, "harmonized_life_cycle_stage", out_path("cluster_by_stage.tsv")))

stage_totals <- study_by_stage["Sum", ]
names(stage_totals)[is.na(names(stage_totals))] <- "NA"
if (IS_PUBLISHED_RELEASE && (
  !identical(as.integer(stage_totals[names(PUBLISHED_STAGE_TOTALS)]), unname(PUBLISHED_STAGE_TOTALS)) ||
    stage_totals[["Sum"]] != PUBLISHED_TOTAL_CELLS ||
    length(stage_totals) != length(PUBLISHED_STAGE_TOTALS) + 1L)) {
  pipeline_fail("study x stage totals differ from the published Table 3 totals")
}

# ---------------------------------------------------------------------------
# Neighbourhood / silhouette metrics
# ---------------------------------------------------------------------------

#' Indices of the k nearest neighbours (Euclidean, self excluded) of every
#' row of `emb`, as an n x k matrix of row indices into `emb`. Also returns
#' the dist object for silhouette reuse.
knn_indices <- function(emb, k) {
  d <- dist(emb)
  dm <- as.matrix(d)
  diag(dm) <- Inf
  nn <- t(apply(dm, 1, function(r) order(r)[seq_len(k)]))
  list(nn = nn, dist = d)
}

#' Mean over cells of the fraction of a cell's neighbours whose label differs
#' from its own. Cells with an NA label are dropped, as are NA-labelled
#' neighbours from each cell's fraction.
mean_frac_other <- function(nn, labels) {
  nb <- matrix(labels[nn], nrow = nrow(nn))
  per_cell <- rowMeans(nb != labels, na.rm = TRUE)
  mean(per_cell[!is.na(labels)], na.rm = TRUE)
}

mean_silhouette <- function(d, labels) {
  keep <- !is.na(labels)
  dm <- as.matrix(d)[keep, keep]
  sil <- cluster::silhouette(as.integer(factor(labels[keep])), dmatrix = dm)
  mean(sil[, "sil_width"])
}

set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
idx <- sample(n_cells, N_SUBSAMPLE)

global_metrics <- lapply(names(EMBEDDINGS), function(e) {
  knn <- knn_indices(pcs[[e]][idx, , drop = FALSE], K_NEIGHBOURS)
  c(
    knn_other_study = mean_frac_other(knn$nn, study[idx]),
    knn_other_stage = mean_frac_other(knn$nn, stage[idx]),
    silhouette_study = mean_silhouette(knn$dist, study[idx]),
    silhouette_stage = mean_silhouette(knn$dist, stage[idx])
  )
})
names(global_metrics) <- names(EMBEDDINGS)

integration_metrics <- data.frame(
  metric = names(global_metrics$unintegrated),
  unintegrated = unname(global_metrics$unintegrated),
  integrated = unname(global_metrics$integrated),
  better_if = c(
    "higher (study mixing)", "lower (stage preservation)",
    "lower (study mixing)", "higher (stage preservation)"
  ),
  n_cells = N_SUBSAMPLE,
  k = K_NEIGHBOURS,
  dims = sprintf("%d:%d", min(DIMS), max(DIMS)),
  seed = SEED
)
utils::write.table(integration_metrics, out_path("integration_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Per-stage cross-study sharing. Continues the same seeded RNG stream as the
# global subsample above; stage order and sampling calls are fixed so the
# published values are reproduced.
#' For every ordered study pair (A, B) within a stage: mean fraction of an
#' A cell's k neighbours that come from B, against the B share of the stage
#' expected under complete mixing. Reuses the per-stage kNN graphs, so it
#' draws no random numbers.
pairwise_sharing <- function(s, knns, labels, p) {
  studies <- names(p)
  rows <- expand.grid(from_study = studies, to_study = studies, stringsAsFactors = FALSE)
  rows <- rows[rows$from_study != rows$to_study, , drop = FALSE]
  frac <- function(nn, a, b) {
    from <- labels == a
    if (!any(from)) return(NA_real_)
    nb <- matrix(labels[nn[from, , drop = FALSE]], nrow = sum(from))
    mean(rowMeans(nb == b))
  }
  data.frame(
    stage = s,
    from_study = rows$from_study,
    to_study = rows$to_study,
    n_from_sampled = vapply(rows$from_study, function(a) sum(labels == a), integer(1)),
    observed_unintegrated = mapply(function(a, b) frac(knns[["unintegrated"]], a, b), rows$from_study, rows$to_study),
    observed_integrated = mapply(function(a, b) frac(knns[["integrated"]], a, b), rows$from_study, rows$to_study),
    expected_complete_mixing = as.numeric(p[rows$to_study]),
    row.names = NULL
  )
}
pairwise <- list()

stage_study_table <- table(stage, study)
shared_stages <- names(which(rowSums(stage_study_table > 0) > 1))

per_stage <- do.call(rbind, lapply(shared_stages, function(s) {
  sel <- which(stage == s)
  n_stage <- length(sel)
  if (n_stage > N_PER_STAGE_MAX) sel <- sample(sel, N_PER_STAGE_MAX)
  knns <- lapply(setNames(names(EMBEDDINGS), names(EMBEDDINGS)), function(e) knn_indices(pcs[[e]][sel, , drop = FALSE], K_NEIGHBOURS)$nn)
  observed <- vapply(names(EMBEDDINGS), function(e) mean_frac_other(knns[[e]], study[sel]), numeric(1))
  p <- prop.table(table(study[which(stage == s)]))
  expected <- 1 - sum(p^2)
  pairwise[[s]] <<- pairwise_sharing(s, knns, study[sel], p)
  data.frame(
    stage = s,
    n_cells_stage = n_stage,
    n_cells_sampled = length(sel),
    n_studies = length(p),
    observed_unintegrated = observed[["unintegrated"]],
    observed_integrated = observed[["integrated"]],
    expected_complete_mixing = expected,
    ratio_unintegrated = observed[["unintegrated"]] / expected,
    ratio_integrated = observed[["integrated"]] / expected,
    check.names = FALSE
  )
}))
utils::write.table(per_stage, out_path("per_stage_cross_study_sharing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pairwise <- do.call(rbind, pairwise)
pairwise$ratio_unintegrated <- pairwise$observed_unintegrated / pairwise$expected_complete_mixing
pairwise$ratio_integrated <- pairwise$observed_integrated / pairwise$expected_complete_mixing
utils::write.table(pairwise, out_path("per_stage_pairwise_study_sharing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# Regression gate against the published values
# ---------------------------------------------------------------------------
check_published <- function(observed, published, keys, label) {
  m <- merge(published, observed, by = keys, suffixes = c(".published", ".observed"), all.x = TRUE)
  diffs <- abs(c(
    m$unintegrated.published - m$unintegrated.observed,
    m$integrated.published - m$integrated.observed
  ))
  if (anyNA(diffs) || any(diffs > PUBLISHED_TOLERANCE)) {
    print(m)
    pipeline_fail(sprintf(
      "%s deviate from the published values by more than %.3f (max |diff| = %s) - the sampling or input differs",
      label, PUBLISHED_TOLERANCE, format(max(diffs), digits = 3)
    ))
  }
  invisible(m)
}
observed_per_stage <- data.frame(stage = per_stage$stage, unintegrated = per_stage$observed_unintegrated, integrated = per_stage$observed_integrated)
if (IS_PUBLISHED_RELEASE) {
  check_published(integration_metrics[, c("metric", "unintegrated", "integrated")], PUBLISHED_GLOBAL, "metric", "global integration metrics")
  check_published(observed_per_stage, PUBLISHED_PER_STAGE, "stage", "per-stage cross-study sharing values")
  if (!setequal(per_stage$stage, PUBLISHED_PER_STAGE$stage)) {
    pipeline_fail("set of stages shared by more than one study differs from the published list")
  }
} else {
  # Report-only comparison with the published 2026-09-21 values.
  cmp <- function(label, observed, published, key) {
    m <- merge(published, observed, by = key, suffixes = c("_published", "_observed"), all = TRUE)
    names(m)[1] <- "key"
    m$diff_unintegrated <- m$unintegrated_observed - m$unintegrated_published
    m$diff_integrated <- m$integrated_observed - m$integrated_published
    cbind(table = label, m)
  }
  utils::write.table(
    rbind(
      cmp("global", integration_metrics[, c("metric", "unintegrated", "integrated")], PUBLISHED_GLOBAL, "metric"),
      cmp("per_stage", observed_per_stage, PUBLISHED_PER_STAGE, "stage")
    ),
    out_path("comparison_to_published.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

# ---------------------------------------------------------------------------
# Figure 2 (PLOS: single combined image, <= 7.5 in wide, Helvetica, >= 8 pt)
# ---------------------------------------------------------------------------
if (!requireNamespace("ggrastr", quietly = TRUE) || !requireNamespace("ragg", quietly = TRUE)) {
  pipeline_fail("ggrastr and ragg are required to rasterise the UMAP point layers (both are in scripts/pixi.lock)")
}

FONT_FAMILY <- "Helvetica"
FIG_WIDTH_IN <- 7.5
FIG_HEIGHT_IN <- 8.0

# Okabe-Ito for studies, Paul Tol "muted" for stages: both colour-blind safe,
# and the two schemes do not share hues so study and stage panels cannot be
# confused. Study colours are fixed by name so every panel agrees.
STUDY_COLOURS <- c(
  "Hazzard2022" = "#E69F00", "Hazzard2024" = "#56B4E9", "Mancio-Silva2022" = "#009E73",
  "Ruberto2022_1" = "#CC79A7", "Ruberto2022_2" = "#D55E00", "Sa2020" = "#0072B2"
)
STAGE_ORDER <- c(
  "Sporozoite", "Hypnozoite", "Schizont (Liver stage)", "EEF",
  "Ring", "Trophozoite", "Schizont (Blood stage)", "Merozoite",
  "Asexual (unresolved stage)", "Female gametocyte", "Male gametocyte"
)
TOL_MUTED <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
STAGE_NA_LABEL <- "Unassigned (NA)"
EMBEDDING_COLOURS <- c("Unintegrated" = "#BBBBBB", "Integrated (Harmony)" = "#004488")

unknown_studies <- setdiff(unique(study), names(STUDY_COLOURS))
if (length(unknown_studies) > 0) pipeline_fail(sprintf("no figure colour for study label(s): %s", paste(unknown_studies, collapse = ", ")))
unknown_stages <- setdiff(unique(stage[!is.na(stage)]), STAGE_ORDER)
if (length(unknown_stages) > 0) pipeline_fail(sprintf("no figure colour for stage(s): %s", paste(unknown_stages, collapse = ", ")))
present_stages <- intersect(STAGE_ORDER, unique(stage))
if (length(present_stages) > length(TOL_MUTED)) pipeline_fail("more stages present than colours in the stage palette")
stage_colours <- c(setNames(TOL_MUTED[seq_along(present_stages)], present_stages), setNames("#DDDDDD", STAGE_NA_LABEL))

stage_display <- ifelse(is.na(stage), STAGE_NA_LABEL, stage)

# Deterministic draw order that interleaves studies in proportion to their
# size, so no study is systematically painted over another (no RNG used).
within_study_rank <- ave(seq_len(n_cells), study, FUN = seq_along)
draw_order <- order(within_study_rank / as.numeric(table(study)[study]))

text_theme <- theme(
  text = element_text(family = FONT_FAMILY, size = 8),
  axis.text = element_text(size = 8),
  axis.title = element_text(size = 8),
  plot.title = element_text(size = 9, face = "bold", hjust = 0),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8, face = "bold"),
  legend.key.size = unit(10, "pt"),
  plot.tag = element_text(size = 11, face = "bold")
)

umap_panel <- function(which_umap, colour_by, title) {
  df <- data.frame(
    x = umaps[[which_umap]][draw_order, 1],
    y = umaps[[which_umap]][draw_order, 2],
    study = factor(study[draw_order], levels = names(STUDY_COLOURS)),
    stage = factor(stage_display[draw_order], levels = names(stage_colours))
  )
  if (colour_by == "study") {
    scale <- scale_colour_manual(name = "Study", values = STUDY_COLOURS, drop = TRUE)
  } else {
    scale <- scale_colour_manual(name = "Harmonized stage", values = stage_colours, drop = TRUE)
  }
  ggplot(df, aes(x = x, y = y, colour = .data[[colour_by]])) +
    ggrastr::rasterise(geom_point(size = 0.05, stroke = 0, shape = 16), dpi = 300, dev = "ragg") +
    scale +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    coord_equal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 8, base_family = FONT_FAMILY) +
    text_theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}

sharing_long <- rbind(
  data.frame(stage = per_stage$stage, embedding = "Unintegrated", ratio = per_stage$ratio_unintegrated),
  data.frame(stage = per_stage$stage, embedding = "Integrated (Harmony)", ratio = per_stage$ratio_integrated)
)
sharing_long$stage <- factor(sharing_long$stage, levels = rev(intersect(STAGE_ORDER, per_stage$stage)))
sharing_long$embedding <- factor(sharing_long$embedding, levels = names(EMBEDDING_COLOURS))
sharing_wide <- data.frame(
  stage = factor(per_stage$stage, levels = levels(sharing_long$stage)),
  from = per_stage$ratio_unintegrated,
  to = per_stage$ratio_integrated
)

panel_e <- ggplot(sharing_long, aes(x = ratio, y = stage)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_segment(data = sharing_wide, aes(x = from, xend = to, y = stage, yend = stage), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(colour = embedding), size = 2.2) +
  scale_colour_manual(name = "Embedding", values = EMBEDDING_COLOURS) +
  scale_x_continuous(limits = c(0, max(1.05, ceiling(max(sharing_long$ratio) * 20) / 20)), breaks = seq(0, 1.5, 0.25), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "Cross-study neighbour sharing within stage",
    x = "Observed / expected fraction of 30 nearest neighbours from another study\n(1 = complete mixing)",
    y = NULL
  ) +
  theme_classic(base_size = 8, base_family = FONT_FAMILY) +
  text_theme +
  theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))

fig2 <- (umap_panel("unintegrated", "study", "Unintegrated: study") |
  umap_panel("unintegrated", "stage", "Unintegrated: stage")) /
  (umap_panel("integrated", "study", "Integrated: study") |
    umap_panel("integrated", "stage", "Integrated: stage")) /
  free(panel_e) +
  plot_layout(heights = c(1, 1, 0.75), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

fig_pdf <- out_path("fig2_integration.pdf")
fig_png <- out_path("fig2_integration.png")
grDevices::pdf(fig_pdf, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, family = FONT_FAMILY, useDingbats = FALSE)
print(fig2)
invisible(grDevices::dev.off())
ragg::agg_png(fig_png, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in", res = 300)
print(fig2)
invisible(grDevices::dev.off())

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
pkg_version <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "absent")
provenance_notes <- sprintf(
  "Read-only post-integration assessment (manuscript Table 3/Fig 2); input=%s input_sha256=%s; subsample=%d k=%d dims=%d:%d; Matrix=%s cluster=%s ggplot2=%s patchwork=%s ggrastr=%s ragg=%s",
  RELEASE_RDS, input_sha256, N_SUBSAMPLE, K_NEIGHBOURS, min(DIMS), max(DIMS),
  pkg_version("Matrix"), pkg_version("cluster"), pkg_version("ggplot2"),
  pkg_version("patchwork"), pkg_version("ggrastr"), pkg_version("ragg")
)
outputs <- out_path(c(
  "study_by_stage.tsv", "study_by_broad_stage.tsv", "cluster_by_study.tsv", "cluster_by_stage.tsv",
  "integration_metrics.tsv", "per_stage_cross_study_sharing.tsv", "per_stage_pairwise_study_sharing.tsv",
  "fig2_integration.pdf", "fig2_integration.png",
  if (!IS_PUBLISHED_RELEASE) "comparison_to_published.tsv"
))
for (f in outputs) {
  record_build_manifest(f, SCRIPT_PATH, cell_count = n_cells, seed = SEED, notes = provenance_notes)
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
options(width = 160)
cat(sprintf("[integration_assessment.R] input: %s (sha256 %s, %d cells)\n", RELEASE_RDS, input_sha256, n_cells))
cat(sprintf("[integration_assessment.R] global metrics (%d-cell subsample, PCs %d-%d, k=%d):\n", N_SUBSAMPLE, min(DIMS), max(DIMS), K_NEIGHBOURS))
print(format(integration_metrics[, c("metric", "unintegrated", "integrated", "better_if")], digits = 3), row.names = FALSE)
cat("[integration_assessment.R] per-stage cross-study sharing (observed / expected under complete mixing):\n")
print(format(per_stage[, c(
  "stage", "n_cells_sampled", "observed_unintegrated", "observed_integrated",
  "expected_complete_mixing", "ratio_unintegrated", "ratio_integrated"
)], digits = 3), row.names = FALSE)
if (IS_PUBLISHED_RELEASE) {
  cat(sprintf("[integration_assessment.R] all values within %.3f of the published values\n", PUBLISHED_TOLERANCE))
} else {
  cat("[integration_assessment.R] non-published release: published values reported in comparison_to_published.tsv, not enforced\n")
}
cat(sprintf("[integration_assessment.R] wrote %d files to %s/ and recorded them in build_manifest.tsv\n", length(outputs), OUT_DIR))
