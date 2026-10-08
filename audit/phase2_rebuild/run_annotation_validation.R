#!/usr/bin/env Rscript
# Phase 2 Part IV/VII: annotation validation + impossible-state check.
# Reads only pv_all_studies.rds (post-singleR.R). Writes only into
# audit/phase2_rebuild/. No production artifact modified.
# Run with: env -u R_LIBS_USER pixi run --manifest-path scripts/pixi.toml \
#   Rscript --vanilla audit/phase2_rebuild/run_annotation_validation.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

.plavisca_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PREPROCESS_ROOT", unset = NA_character_),
    getwd(),
    "/home/sopheap/pvsca_b/pre_process_data"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) stop("Could not locate pre_process_data root.", call. = FALSE)
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

so <- readRDS("pv_all_studies.rds")
meta <- so@meta.data
out_dir <- "audit/phase3_preflight"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------
## annotation_validation.tsv: one row per check, PASS/FAIL + evidence.
## ------------------------------------------------------------------
checks <- list()
add_check <- function(id, description, passed, evidence) {
  checks[[length(checks) + 1]] <<- data.frame(
    check_id = id, description = description,
    status = ifelse(passed, "PASS", "FAIL"), evidence = evidence,
    stringsAsFactors = FALSE
  )
}

non_blood <- meta$tissue_or_sample_type != "Host blood" & !is.na(meta$tissue_or_sample_type)
add_check(
  "AT08", "Non-blood-stage cells receive no blood-IDC/HPI inference",
  !any(!is.na(meta$idc_reference_similarity_label[non_blood])),
  sprintf("%d non-blood cells, %d with non-NA idc_reference_similarity_label", sum(non_blood), sum(!is.na(meta$idc_reference_similarity_label[non_blood])))
)

spz <- !is.na(meta$source_stage_provenance) & meta$source_stage_provenance == "source_selection_defined" &
  !is.na(meta$source_life_cycle_stage) & meta$source_life_cycle_stage == "Sporozoite"
add_check(
  "AT09", "Source-defined sporozoites never broad-labeled Blood stage",
  sum(spz) > 0 && all(meta$parasite_broad_stage[spz] == "Sporozoite stage"),
  sprintf("%d source-defined sporozoite cells, all broad-labeled %s", sum(spz), paste(unique(meta$parasite_broad_stage[spz]), collapse = ","))
)

liver <- meta$parasite_broad_stage == "Liver stage" & !is.na(meta$parasite_broad_stage)
# Gametocyte labels also reach harmonized_life_cycle_stage (source label or
# cluster call), so liver tissue is checked on that field as well.
harmonized_gametocyte <- meta$harmonized_life_cycle_stage %in% c("Male gametocyte", "Female gametocyte")
liver_tissue_gametocyte <- meta$tissue_or_sample_type %in% "Host liver" & harmonized_gametocyte
add_check(
  "AT10", "Liver-stage cells never exposed as gametocytes",
  sum(liver) > 0 && !any(!is.na(meta$pred_gametocyte_sex[liver])) && !any(liver_tissue_gametocyte),
  sprintf("%d liver-stage cells, %d with non-NA pred_gametocyte_sex, %d liver-tissue cells with a gametocyte harmonized stage",
          sum(liver), sum(!is.na(meta$pred_gametocyte_sex[liver])), sum(liver_tissue_gametocyte))
)

near_zero <- meta$total_umi_count < 10
add_check(
  "AT11", "Zero/near-zero expression cells receive no expression-derived annotation",
  !any(!is.na(meta$idc_reference_similarity_label[near_zero])) && !any(!is.na(meta$pred_gametocyte_sex[near_zero])),
  sprintf("%d near-zero (<10 UMI) cells, %d with idc label, %d with gametocyte call", sum(near_zero), sum(!is.na(meta$idc_reference_similarity_label[near_zero])), sum(!is.na(meta$pred_gametocyte_sex[near_zero])))
)

expected_broad <- map_broad_stage(meta$harmonized_life_cycle_stage)
add_check(
  "AT21", "Broad/detailed lifecycle fields always consistent",
  identical(expected_broad, meta$parasite_broad_stage),
  "parasite_broad_stage matches deterministic map_broad_stage(harmonized_life_cycle_stage) for every cell"
)

n_gametocyte <- sum(!is.na(meta$pred_gametocyte_sex))
add_check(
  "gametocyte-presence", "At least one cluster cleared the marker-score threshold for gametocyte identity",
  n_gametocyte > 0,
  sprintf("%d cells called Male/Female gametocyte", n_gametocyte)
)

idc_eligible_count <- sum(!is.na(meta$idc_reference_similarity_label))
add_check(
  "idc-eligibility-sane", "IDC reference-similarity classification restricted to a plausible blood-stage subset",
  idc_eligible_count > 0 && idc_eligible_count < ncol(so),
  sprintf("%d/%d cells received an idc_reference_similarity_label", idc_eligible_count, ncol(so))
)

annotation_validation_tbl <- bind_rows(checks)
write.table(annotation_validation_tbl, file.path(out_dir, "annotation_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------------
## impossible_state_check.tsv: per-cell enumeration of any impossible
## lifecycle-annotation combination found (empty if none - the
## required, expected outcome).
## ------------------------------------------------------------------
liver_or_spz <- meta$parasite_broad_stage %in% c("Liver stage", "Sporozoite stage")
bad_idc <- liver_or_spz & !is.na(meta$idc_reference_similarity_label)
bad_gam <- (liver_or_spz & !is.na(meta$pred_gametocyte_sex)) |
  (harmonized_gametocyte & !(meta$tissue_or_sample_type %in% "Host blood"))
bad_source_provenance <- !is.na(meta$source_stage_provenance) &
  meta$source_stage_provenance == "source_selection_defined" &
  (!is.na(meta$pred_gametocyte_sex) | harmonized_gametocyte)

impossible_cells <- unique(c(
  rownames(meta)[bad_idc], rownames(meta)[bad_gam], rownames(meta)[bad_source_provenance]
))

if (length(impossible_cells) > 0) {
  impossible_state_tbl <- meta[impossible_cells, , drop = FALSE] %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(
      cell_id, study_label, parasite_broad_stage, harmonized_life_cycle_stage,
      source_stage_provenance, idc_reference_similarity_label, pred_gametocyte_sex,
      violation = case_when(
        cell_id %in% rownames(meta)[bad_idc] ~ "Liver/Sporozoite stage with non-NA idc_reference_similarity_label",
        cell_id %in% rownames(meta)[bad_gam] ~ "Liver/Sporozoite stage or non-blood cell with a gametocyte call or label",
        TRUE ~ "source-selection-defined cell with a gametocyte call or label"
      )
    )
} else {
  impossible_state_tbl <- data.frame(
    cell_id = character(0), study_label = character(0), parasite_broad_stage = character(0),
    harmonized_life_cycle_stage = character(0), source_stage_provenance = character(0),
    idc_reference_similarity_label = character(0), pred_gametocyte_sex = character(0),
    violation = character(0)
  )
}
write.table(impossible_state_tbl, file.path(out_dir, "impossible_state_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("=== annotation_validation.tsv ===\n")
print(annotation_validation_tbl)
cat(sprintf("\n=== impossible_state_check.tsv: %d violation(s) found (0 required for release) ===\n", nrow(impossible_state_tbl)))
if (any(annotation_validation_tbl$status == "FAIL") || nrow(impossible_state_tbl) > 0) {
  cat("BLOCKING: at least one annotation validation check failed or an impossible state was found.\n")
} else {
  cat("All annotation validation checks passed; zero impossible states found.\n")
}
