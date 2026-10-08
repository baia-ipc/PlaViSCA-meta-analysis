#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo <- normalizePath(getwd())
# Study objects and the Ruberto author object are recorded relative to the
# repository root and the released outputs relative to the release directory,
# so the manifest stays valid in any checkout of the standard project layout
# (no worktree-specific absolute paths).
candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", "")
input_path <- function(p) if (nzchar(candidate_root)) file.path(candidate_root, p) else p
release_dir <- Sys.getenv("PLAVISCA_RELEASE_DIR", "../deployments/data/releases/plavisca-v1-gametocytes-2026-10-08")
code_commit <- Sys.getenv("PLAVISCA_FINAL_CODE_COMMIT", unset = NA_character_)
if (is.na(code_commit) || !nzchar(code_commit)) stop("Set PLAVISCA_FINAL_CODE_COMMIT to the committed scientific rebuild SHA")
ancestor_check <- system2("git", c("merge-base", "--is-ancestor", code_commit, "HEAD"))
if (ancestor_check != 0L) stop("Frozen technical code commit is not an ancestor of HEAD")

sha256 <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1]
artifact <- function(type, key, path, cells = NA_integer_, study = ".", note = ".") {
  full <- if (grepl("^/", path)) path else file.path(repo, path)
  data.frame(record_type = type, key = key, value = ".", verification_status = "AVAILABLE_AND_VERIFIED",
             path = path, byte_size = file.info(full)$size, sha256 = sha256(full),
             cell_count = cells, study = study, notes = note, stringsAsFactors = FALSE)
}
fact <- function(type, key, value, status = "AVAILABLE_AND_VERIFIED", note = ".") {
  data.frame(record_type = type, key = key, value = as.character(value), verification_status = status,
             path = ".", byte_size = NA_real_, sha256 = ".", cell_count = NA_integer_, study = ".", notes = note,
             stringsAsFactors = FALSE)
}

rows <- list(
  fact("build_identity", "authoritative_manifest", "audit/phase3_preflight/final_release_manifest.tsv"),
  fact("build_identity", "release_branch", "gametocyte-harmonized-stage"),
  fact("build_identity", "release_id", "plavisca-v1-gametocytes-2026-10-08"),
  fact("build_identity", "previous_release", "plavisca-v1-dec01-dec02-2026-09-21", note = "Previous atlas sha256 32068878fc92952e9be539befcc8ff2f42a6818716a912a55f5c45930cc830b9; identical cell membership"),
  fact("build_identity", "previous_release_merge_commit", "117d0b57dea267aaae11a8b148de858f56a6c7f8"),
  fact("build_identity", "parent_release_closure_commit", "dea656faa2b06381cc971b462de8f911df1028f7"),
  fact("build_identity", "final_preprocessing_code_commit", code_commit, note = "Committed scientific rebuild and validation snapshot; the final manifest/report commit is its descendant."),
  fact("build_identity", "phase2_scientific_build_commit", "54c6b9581ff5268422ddffda42ce457bf09c5b55"),
  fact("build_identity", "independent_QA_commit", "a57be450a84d42e47a67d1120320e8a5536797dd"),
  fact("build_identity", "application_commit", Sys.getenv("PLAVISCA_APP_COMMIT", unset = "NOT_RECORDED"), note = "PlaViSca app commit validated against this release (branch gametocyte-harmonized-stage)"),
  fact("build_identity", "release_date", "2026-10-08"),
  fact("build_identity", "annotation_precedence_change", "harmonized_life_cycle_stage = source-defined label > cluster-level gametocyte call > IDC heuristic bin (DEC07 resolved)"),
  fact("environment", "R", "4.3.3"),
  fact("environment", "Seurat", "5.3.0"),
  fact("environment", "SeuratObject", "5.2.0"),
  fact("environment", "harmony", "1.2.3"),
  fact("environment", "SingleR", "2.4.0"),
  artifact("environment", "pixi_manifest", "scripts/pixi.toml", note = "Closure validation environment definition"),
  artifact("environment", "pixi_lock", "scripts/pixi.lock", note = "Resolved closure validation environment"),
  fact("parameter", "integration_seed", "123", note = "Recorded by Phase-2 integration manifest row"),
  fact("parameter", "integration_feature_count", "4709", note = "Intersection of measured panels: STARsolo panel, Mancio-Silva2022 genes, Ruberto2022_1 author_present; after rRNA removal"),
  fact("parameter", "atlas_feature_count", "6811", note = "After removal of 49 GFF rRNA genes (0 rRNA features)"),
  fact("parameter", "rrna_genes_removed", "49", note = "audit/phase2_rebuild/rrna_removal_by_study.tsv"),
  fact("parameter", "HVG_PCA_feature_count", "2000"),
  fact("parameter", "production_cluster_resolution", "0.05"),
  fact("parameter", "AT30_perturbed_resolution", "0.08"),
  fact("parameter", "total_cell_count", "105949"),
  fact("study_count", "Mancio-Silva2022", "1480"),
  fact("study_count", "Sa2020", "9766"),
  fact("study_count", "Hazzard2022", "3294"),
  fact("study_count", "Hazzard2024", "80024"),
  fact("study_count", "Ruberto2022_1", "1438"),
  fact("study_count", "Ruberto2022_2", "9947"),
  fact("provenance", "Ruberto2022_1_source_of_counts", "author_processed_object_raw_umi"),
  fact("provenance", "Ruberto2022_1_counting_pipeline", "kallisto_bustools"),
  fact("provenance", "Ruberto2022_1_reference", "PlasmoDB-51_PvivaxP01_AnnotatedTranscripts"),
  artifact("scientific_input", "Mancio-Silva2022_study_object", input_path("silva2022.rds"), 1480, "Mancio-Silva2022"),
  artifact("scientific_input", "Sa2020_study_object", input_path("sa2020.rds"), 9766, "Sa2020"),
  artifact("scientific_input", "Hazzard2022_study_object", input_path("hazzard2022.rds"), 3294, "Hazzard2022"),
  artifact("scientific_input", "Hazzard2024_study_object", input_path("hazzard2024.rds"), 80024, "Hazzard2024"),
  artifact("scientific_input", "Ruberto2022_1_study_object", input_path("ruberto2022_1.rds"), 1438, "Ruberto2022_1"),
  artifact("scientific_input", "Ruberto2022_2_study_object", input_path("ruberto2022_2.rds"), 9947, "Ruberto2022_2"),
  artifact("scientific_input", "Ruberto2022_1_author_source", input_path("data/Hep59.1.2.seu_20aug2025.rds"), 1438, "Ruberto2022_1"),
  artifact("reference_input", "PlasmoDB68_GFF", "scripts/ref/PlasmoDB-68_PvivaxP01.gff"),
  artifact("reference_input", "STARsolo_feature_panel", "scripts/ref/starsolo_feature_panel.tsv"),
  artifact("reference_input", "Ruberto2022_1_feature_coverage", "audit/ruberto2022_1/phase2_feature_coverage.tsv"),
  artifact("reference_input", "Hazzard2024_authoritative_run_table", "audit/hazzard2024/authoritative_run_table.tsv"),
  artifact("reference_input", "Sa2020_authoritative_sample_crosswalk", "audit/mancio_silva2022/sa2020_authoritative_sample_crosswalk.tsv"),
  artifact("reference_input", "Sa2020_authoritative_cell_lineage", "audit/mancio_silva2022/sa2020_authoritative_cell_lineage.tsv", 9766, "Sa2020"),
  artifact("decision_evidence", "annotation_precedence", "audit/shared_pipeline/annotation_precedence.tsv"),
  artifact("decision_evidence", "AT31_stage_population_floors", "audit/shared_pipeline/stage_population_floors.tsv"),
  artifact("decision_evidence", "DEC01_removal_manifest", "audit/mancio_silva2022/dec01_removed_duplicate_cells.tsv", 14, "Mancio-Silva2022"),
  artifact("reused_input_checksums", "reused_validated_inputs", "audit/phase3_preflight/dec01_dec02_reused_input_checksums.sha256"),
  artifact("scientific_output", "pv_all_studies", file.path(release_dir, "pv_all_studies.rds"), 105949),
  artifact("application_output", "normalize_df", file.path(release_dir, "normalize_df.rds"), 105949),
  artifact("application_output", "raw_df", file.path(release_dir, "raw_df.rds"), 105949),
  artifact("application_output", "scale_df", file.path(release_dir, "scale_df.rds"), 105949),
  artifact("application_output", "cleaned_dataset", file.path(release_dir, "cleaned_dataset.rds"), 105949,
           note = "AT17 canonical lexicographic cell order; non-cell tables unchanged"),
  artifact("application_output", "data_source", file.path(release_dir, "data_source.csv"), 6,
           note = "Six bibliography rows with verified study_pmid")
)

manifest <- do.call(rbind, rows)
write.table(manifest, file.path(repo, "audit/phase3_preflight/final_release_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("Wrote authoritative final release manifest with ", nrow(manifest), " records")
