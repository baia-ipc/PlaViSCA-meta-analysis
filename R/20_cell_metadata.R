# Stage 20: per-cell metadata.
#
# For each study, adds to the cells of stage 10:
#   - study and run information (metadata/studies.tsv, metadata/runs.tsv)
#   - the annotations given by the authors of the study:
#       author_stage        stage assigned by the authors, when they assign one:
#                           Sa2020 gametocytes (S5 Data), Ruberto2022_1 liver
#                           forms, Mancio-Silva2022 liver forms (Table S2), and
#                           all cells of the sporozoite samples
#                           (metadata/runs.tsv, sample_stage)
#       author_sex          Sa2020 gametocyte sex (S5 Data)
#       author_state, author_sex_state, author_ap2g   Mancio-Silva2022 Table S2
#       author_pseudogroup  Hazzard2024 sexual/asexual groups
# (author_cell_id, the cell's name in the authors' data, is set in stage 10.)

STAGE <- "20_cell_metadata"
source("R/lib/setup.R")
source("R/lib/cells.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
})

out <- stage_dir(STAGE)
in_dir <- stage_dir("10_import_counts")
runs <- read_meta("runs.tsv")
studies <- read_meta("studies.tsv")
stage_names <- read_meta("stages.tsv")$stage

# run-level columns copied to every cell of the run
RUN_COLUMNS <- c("host_species", "tissue", "strain", "parasite_origin_location", "infection_design",
  "day_post_infection", "treatment", "replicate_id", "replicate_type", "sc_technology", "sequencer")

# --- author annotations, one function per study -----------------------------

annotate_sa2020 <- function(meta, study_runs) {
  # gametocyte sex from the authors' S5 Data (Sample = "<author sample>_<barcode>")
  s5 <- read_excel(input_path(PARAMS$annotation$sa2020_gametocytes))
  s5 <- s5[s5$Type %in% c("Female", "Male"), ]
  s5_cell <- sample_barcode_cells(s5$Sample, study_runs, "Sa2020 S5 Data")$cell_id
  check_unique(s5_cell, "Sa2020 S5 Data cells")
  found <- s5_cell %in% rownames(meta)
  message(sprintf("Sa2020: %d of %d gametocytes of the S5 Data are among our cells", sum(found), length(found)))
  meta$author_stage <- paste(s5$Type, "gametocyte")[match(rownames(meta), s5_cell)]
  meta$author_sex <- meta$author_stage
  meta
}

annotate_ruberto2022_1 <- function(meta) {
  liver_forms <- c(Hypnozoites = "Hypnozoite", Schizonts = "Schizont (Liver stage)")
  check(all(meta$author_liver_form %in% names(liver_forms)), "unexpected Ruberto2022_1 liver forms: ",
    paste(setdiff(unique(meta$author_liver_form), names(liver_forms)), collapse = ", "))
  meta$author_stage <- unname(liver_forms[meta$author_liver_form])
  meta$author_liver_form <- NULL
  meta
}

annotate_mancio_silva2022 <- function(meta) {
  table_s2 <- read_excel(input_path(PARAMS$import[["Mancio-Silva2022"]]$cell_table), skip = 1)
  i <- match(rownames(meta), table_s2$Updated_names)
  check_no_na(i, "Table S2 row of each Mancio-Silva2022 cell")
  meta$author_state <- table_s2$State[i]
  meta$author_sex_state <- table_s2$sex.state[i]
  meta$author_ap2g <- table_s2$AP2G[i]
  # Table S2 State: "Replicative" are liver schizonts
  meta$author_stage <- ifelse(meta$author_state == "Replicative", "Schizont (Liver stage)", meta$author_state)
  meta
}

annotate_hazzard2024 <- function(meta) {
  processed <- utils::read.delim(input_path(PARAMS$import$Hazzard2024$author_cells))
  meta$author_pseudogroup <- processed$PseudoGroup[match(meta$author_cell_id, rownames(processed))]
  check_no_na(meta$author_pseudogroup, "Hazzard2024 PseudoGroup")
  meta
}

# --- run --------------------------------------------------------------------

for (study in STUDIES) {
  object <- readRDS(file.path(in_dir, paste0(study, ".rds")))
  study_runs <- runs[runs$study_label == study, ]
  meta <- object@meta.data

  meta$study_label <- study
  meta$study_pmid <- studies$study_pmid[studies$study_label == study]
  meta$count_source <- studies$count_source[studies$study_label == study]
  run_info <- study_runs[match(meta$run_id, study_runs$run_id), RUN_COLUMNS]
  run_info$day_post_infection <- as.integer(run_info$day_post_infection)
  rownames(run_info) <- NULL
  meta <- cbind(meta, run_info)

  meta$author_stage <- study_runs$sample_stage[match(meta$run_id, study_runs$run_id)]
  meta <- switch(study,
    Sa2020 = annotate_sa2020(meta, study_runs),
    Ruberto2022_1 = annotate_ruberto2022_1(meta),
    "Mancio-Silva2022" = annotate_mancio_silva2022(meta),
    Hazzard2024 = annotate_hazzard2024(meta),
    meta
  )
  check(all(meta$author_stage %in% c(NA, stage_names)), study, ": author stage not in metadata/stages.tsv: ",
    paste(setdiff(unique(meta$author_stage), c(NA, stage_names)), collapse = ", "))
  check_no_na(meta$tissue, paste(study, "tissue"))

  object@meta.data <- meta
  save_object(object, file.path(out, paste0(study, ".rds")))
  message(sprintf("%s: %d cells; author stage for %d cells", study, ncol(object), sum(!is.na(meta$author_stage))))
}
