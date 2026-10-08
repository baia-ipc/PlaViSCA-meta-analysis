#!/usr/bin/env Rscript
# AT31: every source-defined population (Sa2020 gametocytes, Ruberto2022_1
# hypnozoites/schizonts, Mancio-Silva2022 forms, sporozoites from
# Ruberto2022_2 and Hazzard2022) and each study's harmonized Female/Male
# gametocyte total keeps at least the cell count recorded in
# audit/shared_pipeline/stage_population_floors.tsv. Fails (exit 1) otherwise.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(SeuratObject))
repo <- normalizePath(getwd())
source(file.path(repo, "scripts/pipeline_lib.R"))
candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", repo)
md <- readRDS(file.path(candidate_root, "pv_all_studies.rds"))@meta.data

result <- check_stage_population_floors(md, file.path(repo, STAGE_POPULATION_FLOORS_PATH))
gametocyte_columns <- all(c("Female gametocyte", "Male gametocyte") %in%
  colnames(table(md$study_label, md$harmonized_life_cycle_stage)))
result <- rbind(result[, c("population", "study_label", "harmonized_life_cycle_stage",
                           "source_stage_provenance", "min_cells", "observed_cells", "status")],
  data.frame(population = "table(study_label, harmonized_life_cycle_stage) has both gametocyte columns",
             study_label = "all", harmonized_life_cycle_stage = "Female gametocyte; Male gametocyte",
             source_stage_provenance = "any", min_cells = NA, observed_cells = NA,
             status = if (gametocyte_columns) "PASS" else "FAIL"))
write.table(result, file.path(repo, "audit/phase3_preflight/AT31_stage_population_floors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
print(result, right = FALSE)
if (any(result$status != "PASS")) {
  message("AT31 FAIL"); quit(status = 1)
}
message("AT31 PASS: ", nrow(result), " checks")
