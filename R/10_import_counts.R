# Stage 10: import the count matrices and decide which cells each study has.
#
# Usage: Rscript R/10_import_counts.R [study ...]   (default: all studies)
#
# For each study writes <study>.rds: a Seurat object with the raw UMI counts
# and, per cell, run_id and author_cell_id (the cell's name in the authors'
# data, where the authors provide one).
#
# Studies counted with STARsolo by us (settings in config/params.yaml, import):
#   - read the GeneFull matrix of each run ("filtered" = STARsolo cell calls,
#     "raw" = all barcodes) and remove the rRNA genes
#   - plot the barcode rank (knee) curve of each run
#   - runs are read and made into objects in parallel (import: threads)
#   - Hazzard2022: keep the barcodes with at least a number of unique reads
#     over the whole genome, by tissue (the authors' rule; the reads per
#     barcode come from starsolo/run_starsolo.sh --genome-reads)
#   - Hazzard2024, Ruberto2022_2, Sa2020: keep the cells retained by the
#     authors (Sa2020: the cells of their main analysis, a table in metadata/)
# Studies whose counts come from the authors' matrix:
#   - Ruberto2022_1: raw UMI counts of the authors' object (kallisto|bustools,
#     PlasmoDB-51), on our gene list: genes the authors did not count are
#     added as zero rows (metadata/ruberto2022_1_genes.tsv says which and why)
#   - Mancio-Silva2022: raw UMI counts of the authors' Zenodo object, minus
#     14 cells present twice (metadata/mancio_silva2022_duplicate_cells.tsv)
#
# Alternative entry point: R/tools/import_count_objects.R writes the same
# outputs from existing per-study count objects, skipping the STARsolo import.

STAGE <- "10_import_counts"
source("R/lib/setup.R")
source("R/lib/cells.R")
suppressPackageStartupMessages(library(Seurat))

out <- stage_dir(STAGE)
runs <- read_meta("runs.tsv")
study_table <- read_meta("studies.tsv")
genes <- utils::read.delim(file.path(stage_dir("00_reference"), "genes.tsv"), stringsAsFactors = FALSE)
rrna_genes <- genes$gene_id[genes$is_rrna]

studies <- commandArgs(TRUE)
if (length(studies) == 0) studies <- STUDIES
check(all(studies %in% STUDIES), "unknown study: ", paste(setdiff(studies, STUDIES), collapse = ", "))

# Apply f to each element of x, on PARAMS$import$threads cores; stop if any
# element fails.
parallel_map <- function(x, f) {
  result <- parallel::mclapply(x, f, mc.cores = PARAMS$import$threads, mc.preschedule = FALSE)
  failed <- vapply(result, function(r) is.null(r) || inherits(r, "try-error"), TRUE)
  check(!any(failed), "failed for ", paste(names(x)[failed], collapse = ", "), ": ",
    paste(unique(unlist(lapply(result[failed], as.character))), collapse = "; "))
  result
}

# --- studies counted with STARsolo ------------------------------------------

import_starsolo <- function(study, cfg) {
  suppressPackageStartupMessages(library(DropletUtils))
  study_runs <- runs[runs$study_label == study, ]
  pmid <- study_table$study_pmid[study_table$study_label == study]
  # the import keeps every run of the study listed in metadata/runs.tsv
  matrices <- parallel_map(setNames(study_runs$run_id, study_runs$run_id), function(run) {
    path <- input_path(PARAMS$import$starsolo_counts, pmid, paste0(run, "_solo_out"), "Solo.out", "GeneFull", cfg$matrix)
    sce <- read10xCounts(path, col.names = TRUE)
    sce[!rownames(sce) %in% rrna_genes, ]
  })

  pdf(file.path(out, paste0(study, "_knee_plots.pdf")), width = 12, height = 6)
  par(mfrow = c(1, 2))
  for (run in names(matrices)) {
    ranks <- barcodeRanks(matrices[[run]])
    keep <- !duplicated(ranks$rank)
    plot(ranks$rank[keep], ranks$total[keep], log = "xy", xlab = "Rank", ylab = "Total UMI count", main = run)
    abline(h = S4Vectors::metadata(ranks)$inflection, col = "darkgreen", lty = 2)
    abline(h = S4Vectors::metadata(ranks)$knee, col = "dodgerblue", lty = 2)
    legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "dodgerblue"), lty = 2)
  }
  dev.off()
  record_output(file.path(out, paste0(study, "_knee_plots.pdf")))

  barcodes_read <- vapply(matrices, ncol, 1L)
  if (!is.null(cfg$min_genome_reads)) {
    matrices <- setNames(lapply(names(matrices), function(run) {
      tissue <- study_runs$tissue[study_runs$run_id == run]
      threshold <- cfg$min_genome_reads[[tissue]]
      check(!is.null(threshold), study, ": no min_genome_reads for tissue ", tissue)
      reads <- utils::read.delim(input_path(PARAMS$import$starsolo_counts, pmid, paste0(run, "_solo_out"),
        "genome_reads_per_barcode.tsv"), stringsAsFactors = FALSE)
      cells <- reads$barcode[reads$unique_reads >= threshold]
      check(all(cells %in% colnames(matrices[[run]])), run, ": barcodes of genome_reads_per_barcode.tsv not in the matrix")
      matrices[[run]][, colnames(matrices[[run]]) %in% cells]
    }), names(matrices))
  }

  objects <- parallel_map(setNames(names(matrices), names(matrices)), function(run) {
    counts <- SingleCellExperiment::counts(matrices[[run]])
    colnames(counts) <- cell_id(run, colnames(counts))
    object <- CreateSeuratObject(counts, min.cells = cfg$min_cells, min.features = cfg$min_features)
    object$run_id <- run
    object
  })
  run_names <- names(matrices)
  rm(matrices)
  cells_called <- vapply(objects, function(o) as.integer(ncol(o)), 1L)

  author <- author_cells(study, cfg, study_runs)
  if (!is.null(author)) {
    check_unique(author$cell_id, paste(study, "author cells"))
    missing <- setdiff(author$cell_id, unlist(lapply(objects, colnames), use.names = FALSE))
    check(length(missing) == 0, sprintf("%d %s author cells are not in the STARsolo matrices", length(missing), study))
    # keep only the author cells of each run before merging: merging every
    # barcode of the raw matrices is slow. The genes of each run object stay.
    objects <- lapply(names(objects), function(run) {
      keep <- colnames(objects[[run]])[colnames(objects[[run]]) %in% author$cell_id]
      check(length(keep) > 0, run, " has no ", study, " author cells")
      subset(objects[[run]], cells = keep)
    })
  }

  object <- JoinLayers(Reduce(function(x, y) merge(x, y), unname(objects)))
  object$author_cell_id <- NA_character_
  if (!is.null(author)) {
    object <- subset(object, cells = author$cell_id)
    object$author_cell_id <- author$author_cell_id[match(colnames(object), author$cell_id)]
  }

  list(object = object, summary = data.frame(
    study_label = study, run_id = run_names, matrix = cfg$matrix,
    barcodes_read = as.integer(barcodes_read),
    cells_called = as.integer(cells_called),
    cells_kept = as.integer(table(factor(object$run_id, levels = run_names)))
  ))
}

# --- Ruberto2022_1: authors' matrix -----------------------------------------

import_ruberto2022_1 <- function(cfg) {
  author_object <- readRDS(input_path(cfg$object))
  author_counts <- GetAssayData(author_object, assay = "RNA", layer = "counts")
  check(all(author_counts@x >= 0 & author_counts@x == round(author_counts@x)),
    "the authors' counts layer contains non-integer or negative values")

  # author cell names are "<author sample>_<barcode>", e.g. "51inf_AAAC..."
  author_ids <- colnames(author_counts)
  run <- author_sample_runs(runs[runs$study_label == "Ruberto2022_1", ])[sub("_.*$", "", author_ids)]
  check_no_na(run, "run of each Ruberto2022_1 author cell")
  cell_ids <- cell_id(run, sub("^.*_", "", author_ids))
  check_unique(cell_ids, "Ruberto2022_1 cell names")

  # put the counts on our gene list; genes the authors did not count are zero
  gene_table <- read_meta("ruberto2022_1_genes.tsv")
  in_author <- gene_table$gene_id[gene_table$in_author_matrix == "TRUE"]
  check_same_set(rownames(author_counts), in_author, "genes of the authors' matrix and metadata/ruberto2022_1_genes.tsv")
  counts <- Matrix::Matrix(0, nrow = nrow(gene_table), ncol = length(cell_ids),
    dimnames = list(gene_table$gene_id, cell_ids), sparse = TRUE)
  counts[rownames(author_counts), ] <- author_counts
  message(sprintf("Ruberto2022_1: %d of %d genes counted by the authors, %d set to zero",
    length(in_author), nrow(gene_table), nrow(gene_table) - length(in_author)))

  object <- CreateSeuratObject(counts)
  object$run_id <- unname(run)
  object$author_cell_id <- author_ids
  object$author_liver_form <- author_object@meta.data[author_ids, "LiverForm"]
  list(object = object, summary = data.frame(study_label = "Ruberto2022_1", run_id = names(table(run)),
    matrix = "authors", barcodes_read = NA, cells_called = as.integer(table(run)), cells_kept = as.integer(table(run))))
}

# --- Mancio-Silva2022: authors' matrix --------------------------------------

import_mancio_silva2022 <- function(cfg) {
  author_object <- readRDS(input_path(cfg$object))
  counts <- GetAssayData(author_object, assay = "RNA", layer = "counts")
  # gene names are "PVP01-0000000" followed by optional text; keep the ID
  rownames(counts) <- sub("^([^ -]*-[^ -]*).*", "\\1", rownames(counts))
  object <- CreateSeuratObject(counts)
  n_author <- ncol(object)

  # 14 cells are present twice: the copy with a ".1" suffix is removed
  duplicates <- read_meta("mancio_silva2022_duplicate_cells.tsv")
  check_unique(duplicates$removed_cell_id, "removed Mancio-Silva2022 cells")
  check(all(duplicates$removed_cell_id == paste0(duplicates$kept_cell_id, ".1")), "each removed cell must be the '.1' copy of a kept cell")
  check(all(c(duplicates$removed_cell_id, duplicates$kept_cell_id) %in% colnames(object)),
    "removed and kept Mancio-Silva2022 cells must all be in the authors' object")
  removed_counts <- GetAssayData(object, layer = "counts")[, duplicates$removed_cell_id]
  kept_counts <- GetAssayData(object, layer = "counts")[, duplicates$kept_cell_id]
  check(identical(unname(removed_counts), unname(kept_counts)), "a removed Mancio-Silva2022 cell does not have the same counts as its kept copy")
  object <- subset(object, cells = setdiff(colnames(object), duplicates$removed_cell_id))
  message(sprintf("Mancio-Silva2022: %d cells in the authors' object, %d duplicates removed", n_author, nrow(duplicates)))

  # cell names are "<author sample>_<barcode>"; the sample gives the SRA runs
  sample <- sub("_[ACGTN]+$", "", colnames(object))
  run <- author_sample_runs(runs[runs$study_label == "Mancio-Silva2022", ])[sample]
  check_no_na(run, "run of each Mancio-Silva2022 cell")
  object$run_id <- unname(run)
  object$author_cell_id <- colnames(object)
  list(object = object, summary = data.frame(study_label = "Mancio-Silva2022", run_id = names(table(run)),
    matrix = "authors", barcodes_read = NA, cells_called = NA, cells_kept = as.integer(table(run))))
}

# --- run --------------------------------------------------------------------

summaries <- list()
for (study in studies) {
  cfg <- PARAMS$import[[study]]
  message("importing ", study)
  result <- switch(cfg$source,
    starsolo = import_starsolo(study, cfg),
    author_matrix = switch(study,
      Ruberto2022_1 = import_ruberto2022_1(cfg),
      "Mancio-Silva2022" = import_mancio_silva2022(cfg))
  )
  object <- result$object
  check_unique(colnames(object), paste(study, "cell names"))
  check(all(object$run_id %in% runs$run_id[runs$study_label == study]), study, " cells assigned to a run not listed for it in metadata/runs.tsv")
  save_object(object, file.path(out, paste0(study, ".rds")))
  summaries[[study]] <- result$summary
  message(sprintf("%s: %d cells, %d genes", study, ncol(object), nrow(object)))
}

summary_path <- file.path(out, "cells_per_run.tsv")
all_summaries <- do.call(rbind, summaries)
if (file.exists(summary_path)) {
  previous <- utils::read.delim(summary_path, stringsAsFactors = FALSE)
  all_summaries <- rbind(previous[!previous$study_label %in% studies, ], all_summaries)
}
write_table(all_summaries, summary_path)
