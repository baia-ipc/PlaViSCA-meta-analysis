#!/usr/bin/env Rscript
# AT01-AT31 acceptance matrix. Every status is derived from its evidence
# file (never stamped): run after flatten_data.R and every validation script
# listed in audit/phase3_preflight/REPORT.md, plus the app suite, whose
# results are passed in via PLAVISCA_APP_TEST_RESULTS (run_tests.R TSV).

options(stringsAsFactors = FALSE)
pp <- "audit/phase3_preflight"
out <- file.path(pp, "acceptance_matrix_AT01_AT31.tsv")
ev <- function(f) if (grepl("/", f)) f else file.path(pp, f)
rd <- function(f) read.delim(ev(f), check.names = FALSE, stringsAsFactors = FALSE)
all_pass <- function(f, col = "status") { x <- rd(f); nrow(x) > 0 && all(x[[col]] == "PASS") }
unit_log <- readLines(ev("unit_study_tests.log"))
unit_pass <- function(tag) {
  hit <- grep(sprintf("^\\[%s\\]", tag), unit_log, value = TRUE)
  length(hit) > 0 && all(grepl("^\\[[^]]+\\] PASS", hit))
}
annot <- rd("annotation_validation.tsv")
annot_pass <- function(id) any(annot$check_id == id) && all(annot$status[annot$check_id == id] == "PASS")
summary_tbl <- rd("candidate_observed_summary.tsv")
summary_val <- function(m) as.numeric(summary_tbl$value[summary_tbl$metric == m])
mixing <- rd("batch_mixing_within_comparable_populations.tsv")
multi <- mixing[mixing$n_studies_present > 1, ]
per_study <- rd("per_study_neighborhood_composition.tsv")
repro <- rd("reproducibility_check.tsv")
at29 <- repro[repro$test == "AT29_same_seed_determinism", ]
at30 <- repro[repro$test == "AT30_gametocyte_call_stability_under_perturbation", ]
AT30_MIN_STABILITY <- 0.95
app_results_path <- Sys.getenv("PLAVISCA_APP_TEST_RESULTS", unset = NA_character_)
app <- if (!is.na(app_results_path) && file.exists(app_results_path)) read.delim(app_results_path) else NULL
app_file_pass <- function(pattern) {
  !is.null(app) && any(grepl(pattern, app$file)) && all(app$status[grepl(pattern, app$file)] == "PASS")
}
sa_lineage <- rd("audit/mancio_silva2022/sa2020_authoritative_cell_lineage.tsv")
impossible <- rd("impossible_state_check.tsv")
schema <- read.csv(ev("app_candidate_schema.tsv"), sep = "\t")

checks <- list(
  AT01 = list(all_pass("cell_key_validation.tsv"), "cell_key_validation.tsv", "All cell-indexed exports have unique IDs in identical canonical order."),
  AT02 = list(all_pass("final_study_cell_counts.tsv") && nrow(sa_lineage) == 9766L, "final_study_cell_counts.tsv; audit/mancio_silva2022/sa2020_authoritative_cell_lineage.tsv", "Exact per-study counts; Sa2020 lineage 9,766/9,766."),
  AT03 = list(all_pass("final_study_cell_counts.tsv") && file.exists(ev("membership_vs_previous_release.tsv")) && all_pass("membership_vs_previous_release.tsv"), "final_study_cell_counts.tsv; membership_vs_previous_release.tsv", "Study counts exact and cell-ID set identical to the 2026-09-21 release."),
  AT04 = list(unit_pass("UNIT-01") && !any(grepl("^\\[[^]]+\\] FAIL", unit_log)), "unit_study_tests.log", "Cardinality assertions pass."),
  AT05 = list(unit_pass("UNIT-02"), "unit_study_tests.log", "Keyed metadata, no unsafe recycling."),
  AT06 = list(any(!is.na(schema$literal_NA_count)) && all(schema$literal_NA_count[!is.na(schema$literal_NA_count)] == 0), "app_candidate_schema.tsv", "Zero literal NA strings in exported fields (NA count = column absent, e.g. obsolete names)."),
  AT07 = list(unit_pass("AT07"), "unit_study_tests.log", "Source, inferred and harmonized annotations remain separate fields."),
  AT08 = list(annot_pass("AT08"), "annotation_validation.tsv", "No non-blood IDC label."),
  AT09 = list(annot_pass("AT09"), "annotation_validation.tsv", "Source-defined sporozoites broad-labeled Sporozoite stage."),
  AT10 = list(annot_pass("AT10") && summary_val("impossible_liver_gametocyte") == 0, "annotation_validation.tsv", "No liver gametocyte call or label."),
  AT11 = list(annot_pass("AT11"), "annotation_validation.tsv", "Cells below 10 UMI carry no expression-derived annotation."),
  AT12 = list(all_pass("AT12_hazzard2024_validation.tsv"), "AT12_hazzard2024_validation.tsv", "Hazzard2024 keyed fields 100% concordant."),
  AT13 = list(unit_pass("AT13"), "unit_study_tests.log", "Sa2020 source sex import matches 1,532 Sa2020 cells only."),
  AT14 = list(all_pass("dec02_retained_distinct_cells.tsv") && all_pass("final_study_cell_counts.tsv"), "dec02_retained_distinct_cells.tsv; audit/mancio_silva2022/dec01_removed_duplicate_cells.tsv", "DEC01 removals absent, DEC02 cells retained."),
  AT15 = list(summary_val("ruberto_cells") == 1438 && summary_val("ruberto_near_empty") == 0 && summary_val("ruberto_provenance_complete") == 1438, "candidate_observed_summary.tsv", "Ruberto2022_1 1,438 author-matrix cells with complete provenance."),
  AT16 = list(all_pass("cell_key_validation.tsv"), "cell_key_validation.tsv", "Expression and metadata keys equal before export."),
  AT17 = list(all_pass("AT17_cell_order_validation.tsv"), "AT17_cell_order_validation.tsv", "Exports canonically ordered and aligned."),
  AT18 = list(summary_val("top_genes_studies") == 6, "candidate_observed_summary.tsv", "Top-gene export covers all six studies."),
  AT19 = list(app_file_pass("test_at19_at20"), "app test suite (PLAVISCA_APP_TEST_RESULTS)", "App AT19/AT20 tests pass against this release."),
  AT20 = list(app_file_pass("test_at19_at20"), "app test suite (PLAVISCA_APP_TEST_RESULTS)", "No obsolete stage vocabulary in app code."),
  AT21 = list(annot_pass("AT21"), "annotation_validation.tsv", "Broad stage equals the mapping of harmonized stage."),
  AT22 = list(unit_pass("AT22") && nrow(impossible) == 0, "unit_study_tests.log; impossible_state_check.tsv", "Provenance retained; zero impossible states."),
  AT23 = list(all_pass("AT23_manifest_validation.tsv"), "AT23_manifest_validation.tsv", "Manifest complete; hashes and sizes verified."),
  AT24 = list(nrow(rd("integration_validation.tsv")) > 0 && nrow(per_study) == 6, "integration_validation.tsv", "Per-study pre/post-Harmony neighbour composition regenerated."),
  AT25 = list(all_pass("AT25_marker_preservation_validation.tsv"), "AT25_marker_preservation_validation.tsv", "Marker preservation for both gametocyte sexes, unchanged thresholds."),
  AT26 = list(nrow(multi) > 0 && all(multi$mean_frac_same_study_integrated < multi$mean_frac_same_study_unintegrated), "batch_mixing_within_comparable_populations.tsv", "Same-study neighbour fraction decreases after Harmony in every multi-study broad stage."),
  AT27 = list(all(per_study$retention_ratio_integrated > 1), "per_study_neighborhood_composition.tsv", "Every study keeps neighbourhood enrichment above random."),
  AT28 = list(all_pass("AT28_singleR_uncertainty.tsv"), "AT28_singleR_uncertainty.tsv", "SingleR raw label, delta and pruning flag retained."),
  AT29 = list(nrow(at29) == 1 && isTRUE(as.logical(at29$identical_cluster_labels)), "reproducibility_check.tsv", "Same-graph/resolution clustering rerun identical."),
  AT30 = list(nrow(at30) == 1 && at30$stability_fraction >= AT30_MIN_STABILITY, "reproducibility_check.tsv", sprintf("Gametocyte calls stable under perturbed resolution (>= %.2f).", AT30_MIN_STABILITY)),
  AT31 = list(all_pass("AT31_stage_population_floors.tsv"), "AT31_stage_population_floors.tsv", "Both gametocyte labels and every source-defined population present at or above recorded floors.")
)

tbl <- data.frame(
  test_id = names(checks),
  status = ifelse(vapply(checks, function(x) isTRUE(x[[1]]), logical(1)), "PASS", "FAIL"),
  evidence_file = vapply(checks, function(x) x[[2]], character(1)),
  evidence_detail = vapply(checks, function(x) x[[3]], character(1)),
  row.names = NULL
)
tbl$blocking <- ifelse(tbl$status == "PASS", "NO", "YES")
tbl$evidence_detail[tbl$test_id == "AT30"] <- paste(tbl$evidence_detail[tbl$test_id == "AT30"],
  sprintf("Observed %d/%d.", at30$n_stable_under_perturbation, at30$n_original_gametocyte_cells))
write.table(tbl[, c("test_id", "status", "blocking", "evidence_file", "evidence_detail")], out,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", out, ": ", sum(tbl$status == "PASS"), " PASS / ", sum(tbl$status == "FAIL"), " FAIL")
if (any(tbl$status != "PASS")) quit(status = 1L)
