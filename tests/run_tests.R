# Tests of the metadata tables and of the shared functions. They need no
# data: run with `pixi run test` from the repository root.

Sys.setenv(PLAVISCA_WORKSPACE = Sys.getenv("PLAVISCA_WORKSPACE", tempdir()))
STAGE <- "tests"
source("R/lib/setup.R")
source("R/lib/cells.R")

failures <- 0
test <- function(name, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    message("  error: ", conditionMessage(e))
    FALSE
  })
  cat(if (ok) "ok    " else "FAIL  ", name, "\n", sep = "")
  if (!ok) failures <<- failures + 1
}
fails <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

# --- metadata ---------------------------------------------------------------

studies <- read_meta("studies.tsv")
runs <- read_meta("runs.tsv")
stages <- read_meta("stages.tsv")
columns <- read_meta("columns.tsv")
expected <- read_meta("expected_counts.tsv")

test("study labels are unique", !anyDuplicated(studies$study_label))
test("every run belongs to a known study", all(runs$study_label %in% studies$study_label))
test("every study has runs", all(studies$study_label %in% runs$study_label))
test("run accessions are unique", !anyDuplicated(runs$run_id))
test("run suffixes are unique within each study",
  all(tapply(substring(runs$run_id, nchar(runs$run_id) - 2), runs$study_label, function(x) !anyDuplicated(x))))
test("author sample names are unique within each study",
  all(tapply(runs$author_sample, runs$study_label, function(x) !anyDuplicated(x[!is.na(x)]))))
test("tissues are one of the three tissue values",
  all(runs$tissue %in% c("Host blood", "Host liver", "Vector salivary gland")))
test("sample stages are known stages", all(runs$sample_stage %in% c(NA, stages$stage)))
test("sporozoite samples are salivary gland samples",
  all(runs$tissue[runs$sample_stage %in% "Sporozoite"] == "Vector salivary gland"))
test("every study has import settings", all(studies$study_label %in% names(PARAMS$import)))
test("the merge order lists every study once", setequal(PARAMS$integration$merge_order, studies$study_label) &&
  !anyDuplicated(PARAMS$integration$merge_order))
test("broad stages are blood, liver or sporozoite",
  all(stages$broad_stage %in% c("Blood stage", "Liver stage", "Sporozoite stage")))
test("time point stages are known stages", all(read_meta("idc_timepoints.tsv")$stage %in% stages$stage))
test("gametocyte marker sets are known", all(read_meta("gametocyte_markers.tsv")$marker_set %in%
  c("female", "male", "panel_commitment", "panel_male", "panel_female", "panel_both_sexes")))
test("exported columns are unique", !anyDuplicated(columns$column))
test("expected counts: one total per study", setequal(expected$study_label[expected$check == "study_cells"], studies$study_label))
test("expected counts: stage totals add up to study totals", {
  totals <- setNames(as.integer(expected$cells[expected$check == "study_cells"]), expected$study_label[expected$check == "study_cells"])
  by_stage <- tapply(as.integer(expected$cells[expected$check == "stage_cells"]), expected$study_label[expected$check == "stage_cells"], sum)
  all(by_stage[names(totals)] == totals)
})
test("expected stages are known stages",
  all(expected$stage[expected$check == "stage_cells"] %in% c(stages$stage, "(none)")))
test("removed Mancio-Silva2022 cells are '.1' copies of kept cells", {
  d <- read_meta("mancio_silva2022_duplicate_cells.tsv")
  all(d$removed_cell_id == paste0(d$kept_cell_id, ".1"))
})

# --- shared functions -------------------------------------------------------

test("cell_id uses the last 3 digits of the run", cell_id("SRR11008269", "AAAC") == "269_AAAC")
toy_runs <- data.frame(run_id = c("SRR0000269", "SRR0000270"), author_sample = c("NIH_Sa", "NIH_CQ"))
test("run_of_cell finds the run of a cell", identical(run_of_cell(c("270_A", "269_B"), toy_runs), c("SRR0000270", "SRR0000269")))
test("run_of_cell fails on an unknown prefix", fails(run_of_cell("271_A", toy_runs)))
test("author_sample_runs maps author samples to runs", identical(unname(author_sample_runs(toy_runs)["NIH_CQ"]), "SRR0000270"))
test("Hazzard2024 author cell names are parsed", {
  f <- tempfile()
  writeLines(c("PseudoGroup", "5309_2_AAAC\tSexual", "5309_1_GGGT\tGroupA"), f)
  r <- data.frame(run_id = c("SRR27021955", "SRR27021956"), author_sample = c("5309_2", "5309_1"))
  x <- hazzard2024_author_cells(f, r)
  identical(x$cell_id, c("955_AAAC", "956_GGGT")) && identical(x$author_cell_id, c("5309_2_AAAC", "5309_1_GGGT"))
})
test("check fails when the condition is false", fails(check(FALSE, "x")))
test("check_unique fails on duplicates", fails(check_unique(c("a", "a"), "x")))
test("check_same_order fails on a different order", fails(check_same_order(c("a", "b"), c("b", "a"), "x")))

cat(sprintf("\n%d failure(s)\n", failures))
quit(status = if (failures > 0) 1 else 0)
