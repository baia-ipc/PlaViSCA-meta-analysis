# Alternative entry point: start from per-study count objects instead of the
# STARsolo matrices, skipping the slow part of stage 10.
#
# Usage: Rscript R/tools/import_count_objects.R [study ...]
#        (default: every study with a file in the input folder)
#
# Input: <workspace>/inputs/study_count_objects/<study>.rds, a Seurat object
#        per study whose RNA "counts" layer holds the raw UMI counts of the
#        study's cells, named "<last 3 digits of run>_<barcode>", with rRNA
#        genes already removed. Checksums: SHA256SUMS in the same folder.
#
# Writes the same <study>.rds files as stage 10 (counts, run_id,
# author_cell_id). Where the authors define the cells of a study
# (Hazzard2024, Ruberto2022_2), the cells must be exactly the authors' cells.

STAGE <- "10_import_counts"
source("R/lib/setup.R")
source("R/lib/cells.R")
suppressPackageStartupMessages(library(Seurat))

out <- stage_dir(STAGE)
runs <- read_meta("runs.tsv")
genes <- utils::read.delim(file.path(stage_dir("00_reference"), "genes.tsv"), stringsAsFactors = FALSE)
in_dir <- file.path(INPUTS, "study_count_objects")

studies <- commandArgs(TRUE)
if (length(studies) == 0) studies <- intersect(STUDIES, sub("\\.rds$", "", list.files(in_dir, pattern = "\\.rds$")))
check(all(studies %in% STUDIES), "unknown study: ", paste(setdiff(studies, STUDIES), collapse = ", "))

summaries <- list()
for (study in studies) {
  study_runs <- runs[runs$study_label == study, ]
  input <- readRDS(input_path("study_count_objects", paste0(study, ".rds")))
  object <- CreateSeuratObject(GetAssayData(input, assay = "RNA", layer = "counts"), min.cells = 0, min.features = 0)
  rm(input)
  object$run_id <- run_of_cell(colnames(object), study_runs)
  check(!any(rownames(object) %in% genes$seurat_id[genes$is_rrna]), study, " contains rRNA genes")

  cfg <- PARAMS$import[[study]]
  object$author_cell_id <- NA_character_
  if (!is.null(cfg$author_cells)) {
    author <- switch(study,
      Hazzard2024 = hazzard2024_author_cells(input_path(cfg$author_cells), study_runs),
      Ruberto2022_2 = ruberto2022_2_author_cells(input_path(cfg$author_cells), study_runs)
    )
    check_same_set(colnames(object), author$cell_id, paste(study, "cells and author cells"))
    object$author_cell_id <- author$author_cell_id[match(colnames(object), author$cell_id)]
  }

  save_object(object, file.path(out, paste0(study, ".rds")))
  kept <- table(factor(object$run_id, levels = study_runs$run_id))
  summaries[[study]] <- data.frame(study_label = study, run_id = names(kept), matrix = "count object",
    barcodes_read = NA, cells_called = NA, cells_kept = as.integer(kept))
  message(sprintf("%s: %d cells, %d genes", study, ncol(object), nrow(object)))
}

summary_path <- file.path(out, "cells_per_run.tsv")
all_summaries <- do.call(rbind, summaries)
if (file.exists(summary_path)) {
  previous <- utils::read.delim(summary_path, stringsAsFactors = FALSE)
  all_summaries <- rbind(previous[!previous$study_label %in% studies, ], all_summaries)
}
write_table(all_summaries, summary_path)
