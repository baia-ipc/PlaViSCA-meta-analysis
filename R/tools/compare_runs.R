# Compare the results of two runs: which cells, stages and reported values
# differ. Use it after a change to see what the change did, or before a
# release against the run of the last release.
#
# Usage: Rscript R/tools/compare_runs.R <run A> <run B>
#        (run names under <workspace>/runs/, or paths to run directories)
#
# Reads results/app_data/cleaned_dataset.rds (per-cell table) and
# results/tables_and_values/values.tsv of both runs, and prints:
#   1. cells per study: in both runs, only in A, only in B
#   2. for the cells in both runs, how many changed stage and gametocyte
#      state, and the largest changes of each study (A -> B, number of cells)
#   3. the reported values that differ, or exist in only one run
# Nothing is written.

STAGE <- "tools"
source("R/lib/setup.R")

args <- commandArgs(TRUE)
check(length(args) == 2, "usage: Rscript R/tools/compare_runs.R <run A> <run B>")
run_path <- function(x) if (dir.exists(x)) x else file.path(WORKSPACE, "runs", x)
runs <- setNames(vapply(args, run_path, ""), c("A", "B"))
for (r in names(runs)) check(dir.exists(runs[[r]]), "run ", r, " not found: ", runs[[r]])
cat(sprintf("A: %s\nB: %s\n", runs[["A"]], runs[["B"]]))

read_cells <- function(run) {
  path <- file.path(run, "results", "app_data", "cleaned_dataset.rds")
  check(file.exists(path), "missing: ", path)
  cells <- readRDS(path)$mr_data
  # cell names are unique only within a study
  data.frame(key = paste(cells$study_label, cells$barcode), study_label = cells$study_label,
    stage = coalesce(as.character(cells$stage), "(none)"),
    gametocyte_state = coalesce(as.character(cells$gametocyte_state), "(none)"))
}
read_values <- function(run) {
  path <- file.path(run, "results", "tables_and_values", "values.tsv")
  check(file.exists(path), "missing: ", path)
  utils::read.delim(path, stringsAsFactors = FALSE, colClasses = "character")[, c("key", "formatted")]
}
print_table <- function(df) {
  if (nrow(df) == 0) cat("  (none)\n") else print(df, row.names = FALSE)
}

# --- 1. cells per study -----------------------------------------------------

a <- read_cells(runs[["A"]])
b <- read_cells(runs[["B"]])
cat("\n1. Cells per study\n")
studies <- sort(union(a$study_label, b$study_label))
print_table(data.frame(study_label = studies,
  A = as.integer(table(factor(a$study_label, studies))),
  B = as.integer(table(factor(b$study_label, studies))),
  only_in_A = as.integer(table(factor(a$study_label[!a$key %in% b$key], studies))),
  only_in_B = as.integer(table(factor(b$study_label[!b$key %in% a$key], studies)))))

# --- 2. stage changes -------------------------------------------------------

both <- inner_join(a, b, by = c("key", "study_label"), suffix = c("_A", "_B"))
for (column in c("stage", "gametocyte_state")) {
  from <- both[[paste0(column, "_A")]]
  to <- both[[paste0(column, "_B")]]
  cat(sprintf("\n2. %s: %d of %d cells in both runs changed\n", column, sum(from != to), nrow(both)))
  changes <- data.frame(study_label = both$study_label, A = from, B = to)[from != to, ] |>
    dplyr::count(study_label, A, B, name = "cells") |>
    arrange(study_label, desc(cells))
  # the largest changes of each study; the others summed in one line
  shown <- changes |> group_by(study_label) |> slice_head(n = 10) |> ungroup()
  others <- anti_join(changes, shown, by = c("study_label", "A", "B")) |>
    group_by(study_label) |> summarise(A = "(other changes)", B = "", cells = sum(cells), .groups = "drop")
  print_table(as.data.frame(arrange(bind_rows(shown, others), study_label)))
}

# --- 3. reported values -----------------------------------------------------

values <- full_join(read_values(runs[["A"]]), read_values(runs[["B"]]), by = "key", suffix = c("_A", "_B"))
differ <- values[is.na(values$formatted_A) | is.na(values$formatted_B) | values$formatted_A != values$formatted_B, ]
cat(sprintf("\n3. Reported values: %d of %d differ\n", nrow(differ), nrow(values)))
print_table(setNames(differ, c("key", "A", "B")))
