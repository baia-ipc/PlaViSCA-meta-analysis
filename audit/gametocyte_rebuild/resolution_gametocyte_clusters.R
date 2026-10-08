#!/usr/bin/env Rscript
# Read-only: recompute the singleR.R Part 7 marker-score cluster rule at every
# stored clustering resolution of a post-singleR atlas, and report which
# resolutions yield gametocyte clusters, their size, composition and purity.
# Usage: Rscript resolution_gametocyte_clusters.R <pv_all_studies.rds> <out_dir>
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tidyr) })
args <- commandArgs(TRUE); obj <- readRDS(args[1]); out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
wt <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Production definitions (singleR.R Parts 1 and 7), as in run_AT25_marker_preservation.R.
femaleGams <- (c("PVP01-1207200", "PVP01-0616100", "PVP01-1119300", "PVP01-1465500",
  "PVP01-1259400", "PVP01-1441000", "PVP01-0702600", "PVP01-1306800", "PVP01-0946800", "PVP01-0517400",
  "PVP01-1027600", "PVP01-1003000", "PVP01-1143200", "PVP01-1017500", "PVP01-1024300", "PVP01-0806000",
  "PVP01-1240300", "PVP01-0603600", "PVP01-0712800"))
maleGams <- (c("PVP01-1262200", "PVP01-0530800", "PVP01-1412100", "PVP01-1025600",
  "PVP01-1266500", "PVP01-1229400"))
THR_F <- 0.10; THR_M <- 0.10; MARGIN <- 0.05

obj <- AddModuleScore(obj, features = list(femaleGams), name = "femaleGams_module_score", seed = 1)
obj <- AddModuleScore(obj, features = list(maleGams), name = "maleGams_module_score", seed = 1)
md <- obj@meta.data
md$cell_id <- rownames(md)
md$eligible <- md$tissue_or_sample_type %in% "Host blood" &
  md$source_stage_provenance != "source_selection_defined" & md$total_umi_count >= 10
md$f <- md$femaleGams_module_score1; md$m <- md$maleGams_module_score1

res_cols <- grep("^RNA_snn_res\\.", names(md), value = TRUE)
res_cols <- res_cols[order(as.numeric(sub("RNA_snn_res\\.", "", res_cols)))]
summary_rows <- list(); comp_study <- list(); comp_idc <- list(); calls <- list()
for (rc in res_cols) {
  cl <- as.character(md[[rc]])
  cms <- md %>% mutate(cl = cl) %>% filter(eligible) %>% group_by(cl) %>%
    summarise(n_eligible = n(), mean_female_score = mean(f), mean_male_score = mean(m), .groups = "drop") %>%
    mutate(call = case_when(
      mean_female_score >= THR_F & (mean_female_score - mean_male_score) >= MARGIN ~ "Female gametocyte",
      mean_male_score >= THR_M & (mean_male_score - mean_female_score) >= MARGIN ~ "Male gametocyte",
      TRUE ~ "Asexual"))
  n_total <- table(cl)
  cms$n_total_cells <- as.integer(n_total[cms$cl])
  call_vec <- setNames(cms$call, cms$cl)[cl]; call_vec[!md$eligible | call_vec == "Asexual"] <- NA
  calls[[rc]] <- call_vec
  gam <- cms %>% filter(call != "Asexual")
  for (i in seq_len(nrow(gam))) {
    in_cl <- md$eligible & cl == gam$cl[i]
    own <- if (gam$call[i] == "Female gametocyte") md$f[in_cl] >= THR_F else md$m[in_cl] >= THR_M
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      resolution = sub("RNA_snn_res\\.", "", rc), n_clusters = length(n_total), cluster = gam$cl[i], call = gam$call[i],
      n_total_cells = gam$n_total_cells[i], n_eligible_called = gam$n_eligible[i],
      mean_female_score = round(gam$mean_female_score[i], 4), mean_male_score = round(gam$mean_male_score[i], 4),
      purity_own_sex_score_above_threshold = round(mean(own), 4),
      purity_either_sex_score_above_threshold = round(mean(md$f[in_cl] >= THR_F | md$m[in_cl] >= THR_M), 4))
    comp_study[[length(comp_study) + 1]] <- md[in_cl, ] %>% count(study_label) %>%
      mutate(resolution = sub("RNA_snn_res\\.", "", rc), cluster = gam$cl[i], call = gam$call[i], fraction = round(n / sum(n), 4))
    comp_idc[[length(comp_idc) + 1]] <- md[in_cl, ] %>% count(idc_heuristic_stage_bin) %>%
      mutate(resolution = sub("RNA_snn_res\\.", "", rc), cluster = gam$cl[i], call = gam$call[i], fraction = round(n / sum(n), 4))
  }
  if (nrow(gam) == 0) summary_rows[[length(summary_rows) + 1]] <- data.frame(
    resolution = sub("RNA_snn_res\\.", "", rc), n_clusters = length(n_total), cluster = NA, call = "no gametocyte cluster",
    n_total_cells = NA, n_eligible_called = 0, mean_female_score = NA, mean_male_score = NA,
    purity_own_sex_score_above_threshold = NA, purity_either_sex_score_above_threshold = NA)
}
summ <- bind_rows(summary_rows)
totals <- summ %>% group_by(resolution) %>% summarise(
  female_called = sum(n_eligible_called[call == "Female gametocyte"]),
  male_called = sum(n_eligible_called[call == "Male gametocyte"]), .groups = "drop")
# agreement of each resolution's call with the production resolution (0.05)
ref <- calls[["RNA_snn_res.0.05"]]
totals$agreement_with_res_0.05 <- sapply(paste0("RNA_snn_res.", totals$resolution), function(rc) {
  x <- calls[[rc]]; round(mean((is.na(x) & is.na(ref)) | (!is.na(x) & !is.na(ref) & x == ref)), 4) })
# Sa2020 source agreement at each resolution
sa <- md$study_label == "Sa2020" & !is.na(md$source_sex_annotation)
totals$sa2020_source_agreement <- sapply(paste0("RNA_snn_res.", totals$resolution), function(rc) {
  x <- calls[[rc]][sa]; round(mean(!is.na(x) & x == md$source_sex_annotation[sa]), 4) })
wt(summ, "resolution_gametocyte_clusters.tsv")
wt(totals, "resolution_gametocyte_totals.tsv")
wt(bind_rows(comp_study) %>% select(resolution, cluster, call, study_label, n, fraction), "resolution_gametocyte_cluster_by_study.tsv")
wt(bind_rows(comp_idc) %>% select(resolution, cluster, call, idc_heuristic_stage_bin, n, fraction), "resolution_gametocyte_cluster_by_idc_bin.tsv")
# sanity: the recomputed res-0.05 call must reproduce the released pred_gametocyte_sex
same <- (is.na(ref) & is.na(md$pred_gametocyte_sex)) | (!is.na(ref) & !is.na(md$pred_gametocyte_sex) & ref == md$pred_gametocyte_sex)
cat("res 0.05 recomputed call reproduces released pred_gametocyte_sex:", all(same), sprintf("(%d mismatches)\n", sum(!same)))
options(width = 220); print(as.data.frame(summ)); print(as.data.frame(totals))
