#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Seurat))

repo <- normalizePath(Sys.getenv("PLAVISCA_CANDIDATE_ROOT", getwd()))
previous_root <- Sys.getenv("PLAVISCA_PREVIOUS_CANDIDATE_ROOT", "/home/baia/prj/plavisca/pre_process_data")
out_dir <- file.path(repo, "audit/phase3_preflight")
write_tsv <- function(x, name) {
    write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
                row.names = FALSE, na = "NA")
}
sha256 <- function(path) {
    strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][[1]]
}

so <- readRDS(file.path(repo, "pv_all_studies.rds"))
md <- so@meta.data
expected <- c(
    "Mancio-Silva2022" = 1480L,
    "Sa2020" = 9766L,
    "Hazzard2022" = 3294L,
    "Hazzard2024" = 80024L,
    "Ruberto2022_1" = 1438L,
    "Ruberto2022_2" = 9947L
)
observed <- table(md$study_label)
study_counts <- data.frame(
    study_label = names(expected), expected_cells = unname(expected),
    observed_cells = as.integer(observed[names(expected)])
)
study_counts$status <- ifelse(study_counts$expected_cells == study_counts$observed_cells,
                              "PASS", "FAIL")
write_tsv(study_counts, "final_study_cell_counts.tsv")

dec01 <- read.delim(file.path(repo, "audit/mancio_silva2022/dec01_removed_duplicate_cells.tsv"))
d5 <- c("D5Seq1_CCCCGATTGACG", "D5Seq2_CCCCGATTGACG")
d5_md <- md[d5, , drop = FALSE]
dec02 <- data.frame(
    cell_id = d5,
    present = d5 %in% rownames(md),
    study_label = d5_md$study_label,
    source_sample_id = d5_md$source_orig_ident,
    run_id = d5_md$run_id,
    pair_specific_duplicate_or_suspicion_flag_present = vapply(d5, function(id) {
        any(grepl("d5.*duplicate|duplicate.*d5|suspicion", names(md), ignore.case = TRUE))
    }, logical(1)),
    decision_id = "DEC02"
)
dec02$status <- ifelse(dec02$present &
    !dec02$pair_specific_duplicate_or_suspicion_flag_present, "PASS", "FAIL")
write_tsv(dec02, "dec02_retained_distinct_cells.tsv")

embedding_names <- c("pca_unintegrated", "pca_integrated", "umap_unintegrated",
                     "umap_integrated", "tsne_integrated")
embedding_hashes <- do.call(rbind, lapply(embedding_names, function(nm) {
    path <- tempfile(pattern = paste0(nm, "_"), fileext = ".rds")
    saveRDS(Embeddings(so, nm), path, version = 3)
    row <- data.frame(reduction = nm, dimensions = paste(dim(Embeddings(so, nm)), collapse = "x"),
                      sha256 = sha256(path))
    unlink(path)
    row
}))
write_tsv(embedding_hashes, "final_embedding_hashes.tsv")

# The previous (phase-2/closure, 105,963-cell) atlas is archived in
# archive/phase2-atlas-pv_all_studies-2026-09-20.tar; extract it and point
# PLAVISCA_PREVIOUS_CANDIDATE_ROOT at the directory holding pv_all_studies.rds.
previous_atlas <- file.path(previous_root, "pv_all_studies.rds")
if (!file.exists(previous_atlas)) {
    stop("previous atlas not found: ", previous_atlas,
         " - extract archive/phase2-atlas-pv_all_studies-2026-09-20.tar and set PLAVISCA_PREVIOUS_CANDIDATE_ROOT")
}
old_so <- readRDS(previous_atlas)
old_md <- old_so@meta.data
comparison <- data.frame(
    metric = c("total_cells", paste0("study_cells:", names(expected)),
               "DEC01_removed_ids", "DEC02_retained_cells",
               "paired_missing_stage_cells"),
    previous = c(ncol(old_so), as.integer(table(old_md$study_label)[names(expected)]),
                 0L, sum(d5 %in% rownames(old_md)),
                 sum(is.na(old_md$parasite_broad_stage) & is.na(old_md$harmonized_life_cycle_stage))),
    final = c(ncol(so), as.integer(observed[names(expected)]),
              sum(!dec01$removed_cell_id %in% rownames(md)), sum(d5 %in% rownames(md)),
              sum(is.na(md$parasite_broad_stage) & is.na(md$harmonized_life_cycle_stage)))
)
comparison$classification <- c(
    "EXPECTED_DEC01_MEMBERSHIP_CHANGE",
    ifelse(names(expected) == "Mancio-Silva2022", "EXPECTED_DEC01_MEMBERSHIP_CHANGE", "EXPECTED_UNCHANGED"),
    "EXPECTED_DEC01_REMOVAL", "EXPECTED_DEC02_RETAINED", "EXPECTED_DEC01_NULLABILITY_CHANGE"
)
write_tsv(comparison, "previous_vs_new_candidate_comparison.tsv")

stopifnot(
    ncol(so) == 105949L,
    all(study_counts$status == "PASS"),
    nrow(dec01) == 14L,
    !anyDuplicated(dec01$removed_cell_id),
    !any(dec01$removed_cell_id %in% rownames(md)),
    all(dec01$retained_counterpart %in% rownames(md)),
    all(dec02$status == "PASS"),
    length(unique(dec02$source_sample_id)) == 2L,
    length(unique(dec02$run_id)) == 2L
)
message("DEC01/DEC02 final membership validation PASS")
