# Stage 70: tables, values and figures for reporting the results.
#
# Collects the numbers produced by the earlier stages in a form ready to be
# cited (e.g. in the manuscript, which is edited separately using these
# files). Writes to results/tables_and_values/:
#   values.tsv / values.tex       the numbers cited in the text; values.tex
#                                 defines one LaTeX command per value
#                                 (\pvCellsTotal, ...)
#   table_datasets.tsv, table_datasets_rows.tex          studies, accessions, cells
#   table_stage_counts.tsv, table_stage_counts_rows.tex  cells per study and stage
#   starsolo_parameters.tsv       STAR version and options, read from the
#                                 STARsolo logs
#   parameters.tsv                config/params.yaml as a table
#   software_versions.tsv         versions of R and of the R packages used
#   columns.tsv                   description of the per-cell columns
#   figure_integration.pdf        copy of results/integration/integration.pdf
#   supplementary_qc.pdf          UMIs and genes per cell, per study
#   supplementary_gametocyte_scores.pdf  module scores per cluster
#   supplementary_idc_delta.pdf   copy of the SingleR delta distribution

STAGE <- "70_tables_and_values"
source("R/lib/setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})

out <- results_dir("tables_and_values")
annotation_dir <- file.path(RESULTS, "stage_annotation")
integration_dir <- file.path(RESULTS, "integration")
read_result <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

studies <- read_meta("studies.tsv")
runs <- read_meta("runs.tsv")
stages <- read_meta("stages.tsv")
mr_data <- readRDS(file.path(RESULTS, "app_data", "cleaned_dataset.rds"))$mr_data

# --- values -----------------------------------------------------------------

values <- list()
add <- function(key, value, description, source, digits = NA, percent = FALSE) {
  check(grepl("^[A-Za-z]+$", key), "value key must be letters only (LaTeX command name): ", key)
  check(is.null(values[[key]]), "value key defined twice: ", key)
  formatted <- if (percent) sprintf("%.0f\\%%", 100 * value)
    else if (!is.na(digits)) paste0(if (value < 0) "$-$" else "", formatC(abs(value), format = "f", digits = digits))
    else formatC(value, format = "d", big.mark = ",")
  values[[key]] <<- data.frame(key = key, value = value, formatted = formatted, description = description, source = source)
}

# cells and genes
add("Studies", length(unique(mr_data$study_label)), "studies", "results/app_data")
add("CellsTotal", nrow(mr_data), "cells in the atlas", "results/app_data")
atlas_genes <- ncol(readRDS(file.path(RESULTS, "app_data", "raw_df.rds"))) - nrow(read_meta("columns.tsv"))
add("GenesTotal", atlas_genes, "genes in the atlas", "results/app_data/raw_df.rds")

measured <- read_result(file.path(stage_dir("30_integrate"), "measured_genes.tsv"))
add("GenesShared", sum(measured$in_all_studies), "genes measured in all studies (variable genes chosen among these)", "outs/30_integrate/measured_genes.tsv")
add("GenesRubertoOne", sum(measured$Ruberto2022_1), "genes measured in Ruberto2022_1 (authors' matrix)", "outs/30_integrate/measured_genes.tsv")
add("GenesMancioSilva", sum(measured[["Mancio-Silva2022"]]), "genes measured in Mancio-Silva2022 (authors' matrix, rRNA removed)", "outs/30_integrate/measured_genes.tsv")
add("VariableGenes", PARAMS$integration$n_variable_genes, "variable genes", "config/params.yaml")

rrna <- read_result(file.path(stage_dir("30_integrate"), "rrna_genes_removed.tsv"))
ms_rrna <- rrna[rrna$study_label == "Mancio-Silva2022", ]
add("RrnaGenesRemoved", sum(rrna$rrna_genes_removed), "rRNA genes removed from the authors' matrices", "outs/30_integrate/rrna_genes_removed.tsv")
add("RrnaShareMancioSilva", ms_rrna$rrna_umi_removed / ms_rrna$total_umi_before, "share of Mancio-Silva2022 UMIs in rRNA genes", "outs/30_integrate/rrna_genes_removed.tsv", percent = TRUE)
add("RrnaGenesGff", sum(utils::read.delim(file.path(stage_dir("00_reference"), "genes.tsv"))$is_rrna), "rRNA genes in the annotation", "outs/00_reference/genes.tsv")
add("DuplicatesMancioSilva", nrow(read_meta("mancio_silva2022_duplicate_cells.tsv")), "Mancio-Silva2022 cells present twice (one copy removed)", "metadata/mancio_silva2022_duplicate_cells.tsv")

# stages
broad <- table(mr_data$broad_stage)
broad_studies <- tapply(mr_data$study_label, mr_data$broad_stage, function(x) length(unique(x)))
for (b in names(broad)) {
  key <- gsub("[^A-Za-z]", "", tools::toTitleCase(b))
  add(paste0("Cells", key), as.integer(broad[[b]]), paste("cells,", b), "results/integration/study_by_broad_stage.tsv")
  add(paste0("Studies", key), as.integer(broad_studies[[b]]), paste("studies with cells of", b), "results/integration/study_by_broad_stage.tsv")
}
stage_n <- table(mr_data$stage)
for (s in names(stage_n)) {
  key <- paste0("Cells", gsub("[^A-Za-z]", "", tools::toTitleCase(s)))
  add(key, as.integer(stage_n[[s]]), paste("cells,", s), "results/integration/study_by_stage.tsv")
}
add("CellsNoStage", sum(is.na(mr_data$stage)), "cells without a stage", "results/integration/study_by_stage.tsv")

# annotation
eligible <- read_result(file.path(annotation_dir, "eligible_cells.tsv"))
add("CellsEligible", sum(eligible$cells[eligible$eligible]), "blood cells eligible for stage inference", "results/stage_annotation/eligible_cells.tsv")
add("MinUmi", PARAMS$annotation$idc_min_umi, "minimum UMIs for stage inference", "config/params.yaml")
reference <- read_result(file.path(annotation_dir, "idc_reference.tsv"))
add("ReferenceGenes", reference$reference_genes, "genes of the reference time course", "results/stage_annotation/idc_reference.tsv")
add("ReferenceGenesMapped", reference$with_pvp01_ortholog, "reference genes with a PvP01 ortholog", "results/stage_annotation/idc_reference.tsv")
add("ReferenceGenesUnmapped", reference$without_ortholog, "reference genes without a PvP01 ortholog", "results/stage_annotation/idc_reference.tsv")
markers <- read_meta("gametocyte_markers.tsv")
add("MarkersFemale", sum(markers$marker_set == "female"), "female gametocyte marker genes", "metadata/gametocyte_markers.tsv")
add("MarkersMale", sum(markers$marker_set == "male"), "male gametocyte marker genes", "metadata/gametocyte_markers.tsv")
add("GametocyteMinScore", PARAMS$annotation$gametocyte_female_min_score, "minimum mean module score of a gametocyte cluster", "config/params.yaml", digits = 2)
add("GametocyteMinMargin", PARAMS$annotation$gametocyte_min_margin, "minimum margin over the other sex", "config/params.yaml", digits = 2)
add("ClusterResolution", PARAMS$blood_embedding$cluster_resolution, "resolution of the blood-stage clusters used for gametocytes", "config/params.yaml", digits = 2)
add("AtlasClusterResolution", PARAMS$annotation$atlas_cluster_resolution, "resolution of the atlas clusters", "config/params.yaml", digits = 2)
add("GametocyteMaleMinScore", PARAMS$annotation$gametocyte_male_min_score, "minimum mean male score of a male gametocyte cluster", "config/params.yaml", digits = 2)
add("LateMaleMinScore", PARAMS$annotation$late_male_min_score, "minimum mean mature-male panel score of a late male cluster", "config/params.yaml", digits = 2)
add("UndecidedMinMaleScore", PARAMS$annotation$undecided_min_male_score, "minimum mean male score for a female-branch cluster to be left without sex", "config/params.yaml", digits = 2)
add("DevelopingMaxUmi", PARAMS$annotation$developing_max_umi, "median UMIs below which a female-branch cluster is developing gametocytes", "config/params.yaml")
add("LateFemaleMinUmi", PARAMS$annotation$late_female_min_umi, "median UMIs from which a female-branch cluster is late female", "config/params.yaml")
add("MarkersPanelMale", sum(markers$marker_set == "panel_male"), "genes of the mature-male panel", "metadata/gametocyte_markers.tsv")
gstate <- table(mr_data$gametocyte_state)
for (g in names(gstate)) add(paste0("CellsState", gsub("[^A-Za-z]", "", tools::toTitleCase(g))), as.integer(gstate[[g]]), paste("cells with gametocyte_state", g), "results/app_data")
add("FemaleMaleRatio", sum(mr_data$stage == "Female gametocyte", na.rm = TRUE) / sum(mr_data$stage == "Male gametocyte", na.rm = TRUE),
  "female to male gametocytes", "results/app_data", digits = 1)

calls <- read_result(file.path(annotation_dir, "author_labels_vs_cluster_call.tsv"))
n_calls <- function(study, author, call = NULL) {
  sel <- calls$study_label == study & calls$author_label %in% author
  if (!is.null(call)) sel <- sel & calls$cluster_gametocyte_call %in% call
  sum(calls$cells[sel])
}
add("SaFemaleAuthor", n_calls("Sa2020", "Female gametocyte"), "Sa2020 cells labelled female gametocytes by the authors", "results/stage_annotation/author_labels_vs_cluster_call.tsv")
add("SaFemaleAgree", n_calls("Sa2020", "Female gametocyte", "Female gametocyte"), "... of which in a female gametocyte cluster", "results/stage_annotation/author_labels_vs_cluster_call.tsv")
add("SaFemaleDeveloping", n_calls("Sa2020", "Female gametocyte", "Developing gametocyte"), "... of which in a developing gametocyte cluster", "results/stage_annotation/author_labels_vs_cluster_call.tsv")
add("SaFemaleInMale", n_calls("Sa2020", "Female gametocyte", "Male gametocyte"), "... of which in a male gametocyte cluster", "results/stage_annotation/author_labels_vs_cluster_call.tsv")
add("SaMaleAuthor", n_calls("Sa2020", "Male gametocyte"), "Sa2020 cells labelled male gametocytes by the authors", "results/stage_annotation/author_labels_vs_cluster_call.tsv")
add("SaMaleAgree", n_calls("Sa2020", "Male gametocyte", "Male gametocyte"), "... of which in a male gametocyte cluster", "results/stage_annotation/author_labels_vs_cluster_call.tsv")
sexual <- n_calls("Hazzard2024", "Sexual")
asexual_groups <- setdiff(unique(calls$author_label[calls$study_label == "Hazzard2024"]), "Sexual")
asexual <- n_calls("Hazzard2024", asexual_groups)
add("HazSexual", sexual, "Hazzard2024 cells in the authors' sexual group", "results/stage_annotation/author_labels_vs_cluster_call.tsv")
add("HazSexualCalled", n_calls("Hazzard2024", "Sexual", c("Female gametocyte", "Male gametocyte")) / sexual,
  "share of them in a gametocyte cluster", "results/stage_annotation/author_labels_vs_cluster_call.tsv", percent = TRUE)
add("HazAsexualCalled", n_calls("Hazzard2024", asexual_groups, c("Female gametocyte", "Male gametocyte")) / asexual,
  "share of the authors' asexual groups in a gametocyte cluster", "results/stage_annotation/author_labels_vs_cluster_call.tsv", percent = TRUE)

# integration assessment
a <- PARAMS$assessment
add("AssessCells", a$n_cells, "cells subsampled for the global metrics", "config/params.yaml")
add("AssessNeighbours", a$k_neighbours, "nearest neighbours", "config/params.yaml")
add("AssessPcs", a$n_pcs, "principal components", "config/params.yaml")
add("AssessMaxStage", a$max_cells_per_stage, "maximum cells per stage for per-stage sharing", "config/params.yaml")
metrics <- read_result(file.path(integration_dir, "integration_metrics.tsv"))
metric_keys <- c(knn_other_study = "KnnStudy", knn_other_stage = "KnnStage", silhouette_study = "SilStudy", silhouette_stage = "SilStage")
for (i in seq_len(nrow(metrics))) {
  key <- metric_keys[[metrics$metric[i]]]
  add(paste0(key, "Before"), metrics$unintegrated[i], paste(metrics$metric[i], "before integration"), "results/integration/integration_metrics.tsv", digits = 2)
  add(paste0(key, "After"), metrics$integrated[i], paste(metrics$metric[i], "after integration"), "results/integration/integration_metrics.tsv", digits = 2)
}
sharing <- read_result(file.path(integration_dir, "per_stage_sharing.tsv"))
for (i in seq_len(nrow(sharing))) {
  key <- paste0("Share", gsub("[^A-Za-z]", "", tools::toTitleCase(sharing$stage[i])))
  add(paste0(key, "Before"), sharing$ratio_unintegrated[i], paste("sharing ratio before integration,", sharing$stage[i]), "results/integration/per_stage_sharing.tsv", digits = 2)
  add(paste0(key, "After"), sharing$ratio_integrated[i], paste("sharing ratio after integration,", sharing$stage[i]), "results/integration/per_stage_sharing.tsv", digits = 2)
}

# --- software and STARsolo --------------------------------------------------

packages <- c("Seurat", "SeuratObject", "harmony", "SingleR", "DropletUtils", "uwot", "Rtsne", "cluster", "Matrix", "ggplot2")
software <- data.frame(software = c("R", packages),
  version = c(paste(R.version$major, R.version$minor, sep = "."), vapply(packages, function(x) as.character(packageVersion(x)), "")))

star_options <- c("readFilesIn", "soloType", "soloCBwhitelist", "soloCBlen", "soloUMIlen", "soloBarcodeReadLength",
  "soloStrand", "soloUMIfiltering", "soloFeatures", "alignIntronMin", "alignIntronMax", "sjdbGTFfile")
starsolo_studies <- studies$study_label[vapply(studies$study_label, function(s) identical(PARAMS$import[[s]]$source, "starsolo"), TRUE)]
starsolo <- do.call(rbind, lapply(starsolo_studies, function(study) {
  pmid <- studies$study_pmid[studies$study_label == study]
  logs <- file.path(INPUTS, PARAMS$import$starsolo_counts, pmid,
    paste0(runs$run_id[runs$study_label == study], "_solo_out"), "Log.out")
  logs <- logs[file.exists(logs)]
  if (length(logs) == 0) return(NULL)
  per_run <- lapply(logs, function(f) {
    lines <- readLines(f, n = 200)
    command <- lines[which(lines == "##### Command Line:") + 1]
    tokens <- strsplit(trimws(command), " +")[[1]]
    option_values <- vapply(star_options, function(o) {
      i <- match(paste0("--", o), tokens)
      if (is.na(i)) return(NA_character_)
      j <- i + 1
      v <- character()
      while (j <= length(tokens) && !startsWith(tokens[j], "--")) { v <- c(v, basename(tokens[j])); j <- j + 1 }
      # read files: keep only the read number (R1, R2), not the run name
      if (o == "readFilesIn") v <- paste0("R", sub("^.*_([12])\\.fastq(\\.gz)?$", "\\1", v))
      paste(v, collapse = " ")
    }, "")
    c(star_version = sub("^STAR version=", "", grep("^STAR version=", lines, value = TRUE)[1]), option_values)
  })
  per_run <- unique(do.call(rbind, per_run))
  check(nrow(per_run) == 1, "STARsolo options differ between the runs of ", study)
  data.frame(study_label = study, runs = length(logs), per_run, check.names = FALSE)
}))
write_table(starsolo, file.path(out, "starsolo_parameters.tsv"))
star_versions <- unique(starsolo$star_version)
check(length(star_versions) == 1, "more than one STAR version: ", paste(star_versions, collapse = ", "))
software <- rbind(data.frame(software = "STAR", version = star_versions), software)
write_table(software, file.path(out, "software_versions.tsv"))
for (i in seq_len(nrow(software))) {
  values[[paste0("Version", software$software[i])]] <- data.frame(key = paste0("Version", gsub("[^A-Za-z]", "", software$software[i])),
    value = NA, formatted = software$version[i], description = paste(software$software[i], "version"), source = "software_versions.tsv")
}

params_table <- function(x, prefix = "") {
  do.call(rbind, lapply(names(x), function(n) {
    v <- x[[n]]
    if (is.list(v)) params_table(v, paste0(prefix, n, "."))
    else data.frame(parameter = paste0(prefix, n), value = paste(v, collapse = ", "))
  }))
}
write_table(params_table(PARAMS), file.path(out, "parameters.tsv"))

values <- do.call(rbind, values)
check_unique(values$key, "value keys")
write_table(values, file.path(out, "values.tsv"))
writeLines(c("% Generated by R/70_tables_and_values.R", sprintf("\\newcommand{\\pv%s}{%s}", values$key, values$formatted)),
  file.path(out, "values.tex"))
record_output(file.path(out, "values.tex"))

# --- tables -----------------------------------------------------------------
# The .tex files hold the body rows of the tables (header and caption are
# written where the table is used).

comma <- function(x) formatC(x, format = "d", big.mark = ",")
tex_escape <- function(x) gsub("([&%_#])", "\\\\\\1", x)
tex_rows <- function(cells, path) {
  writeLines(c("% Generated by R/70_tables_and_values.R", paste0(apply(cells, 1, paste, collapse = " & "), " \\\\")), path)
  record_output(path)
}
study_name <- setNames(sprintf("%s \\textit{et al.}, %s%s", studies$first_author, substr(studies$publication_date, 1, 4),
  ifelse(grepl("_[0-9]$", studies$study_label), paste0(" (", sub("^.*_", "", studies$study_label), ")"), "")),
  studies$study_label)

chemistry <- c("737K-august-2016.txt" = "10x (v2)", "3M-february-2018.txt" = "10x (v3)")
platform <- setNames(chemistry[starsolo$soloCBwhitelist], starsolo$study_label)
platform[["Mancio-Silva2022"]] <- "Seq-Well"
datasets <- studies |>
  transmute(study_label, sample_types, platform = unname(platform[study_label]), project_accession,
    cells = as.integer(table(mr_data$study_label)[study_label]))
check_no_na(datasets$platform, "platform of each study")
write_table(datasets, file.path(out, "table_datasets.tsv"))
tex_rows(rbind(
  cbind(sprintf("%s~\\cite{%s}", study_name[datasets$study_label], datasets$study_label), tex_escape(datasets$sample_types),
    tex_escape(datasets$platform), datasets$project_accession, comma(datasets$cells)),
  c("\\midrule\nTotal", "", "", "", comma(sum(datasets$cells)))), file.path(out, "table_datasets_rows.tex"))

# cells per study and stage; studies grouped by compartment (blood, liver,
# sporozoite), in the order of metadata/studies.tsv within a group
stage_order <- stages$stage[stages$stage %in% mr_data$stage]
counts <- as.data.frame.matrix(table(factor(mr_data$study_label, levels = studies$study_label), factor(mr_data$stage, levels = stage_order)))
counts$Total <- as.integer(table(mr_data$study_label)[rownames(counts)])
first_stage <- apply(counts[, stage_order] > 0, 1, which.max)
counts <- counts[order(first_stage, match(rownames(counts), studies$study_label)), ]
write_table(cbind(study_label = rownames(counts), counts), file.path(out, "table_stage_counts.tsv"))
unlabelled <- counts$Total > rowSums(counts[, stage_order])
tex_rows(rbind(
  cbind(paste0(study_name[rownames(counts)], ifelse(unlabelled, "\\textsuperscript{a}", "")), as.matrix(as.data.frame(lapply(counts, comma)))),
  c("\\midrule\nTotal", comma(colSums(counts)))), file.path(out, "table_stage_counts_rows.tex"))

columns <- read_meta("columns.tsv")
write_table(columns, file.path(out, "columns.tsv"))

# --- figures ----------------------------------------------------------------

invisible(file.copy(file.path(integration_dir, "integration.pdf"), file.path(out, "figure_integration.pdf"), overwrite = TRUE))
record_output(file.path(out, "figure_integration.pdf"))
invisible(file.copy(file.path(annotation_dir, "idc_delta_distribution.pdf"), file.path(out, "supplementary_idc_delta.pdf"), overwrite = TRUE))
record_output(file.path(out, "supplementary_idc_delta.pdf"))

qc <- mr_data |>
  select(study_label, UMIs = n_count_rna, Genes = n_feature_rna) |>
  pivot_longer(c(UMIs, Genes), names_to = "measure")
ggsave(file.path(out, "supplementary_qc.pdf"), width = 7.5, height = 4,
  ggplot(qc, aes(study_label, value)) + geom_violin(scale = "width", fill = "grey85") +
    geom_boxplot(width = 0.15, outlier.shape = NA) + scale_y_log10() + facet_wrap(~measure, scales = "free_y") +
    labs(x = NULL, y = "per cell (log scale)") + theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
record_output(file.path(out, "supplementary_qc.pdf"))

scores <- mr_data |>
  filter(tissue == "Host blood", n_count_rna >= PARAMS$annotation$idc_min_umi) |>
  select(cluster, Female = female_gam_score, Male = male_gam_score) |>
  pivot_longer(c(Female, Male), names_to = "markers", values_to = "score")
ggsave(file.path(out, "supplementary_gametocyte_scores.pdf"), width = 7.5, height = 3.5,
  ggplot(scores, aes(cluster, score, fill = markers)) + geom_boxplot(outlier.size = 0.2) +
    geom_hline(yintercept = PARAMS$annotation$gametocyte_min_score, linetype = "dashed") +
    scale_fill_manual(values = c(Female = "#DAA520", Male = "#1F78B4")) +
    labs(x = "Cluster", y = "Module score", fill = "Marker set") + theme_classic(base_size = 8))
record_output(file.path(out, "supplementary_gametocyte_scores.pdf"))
