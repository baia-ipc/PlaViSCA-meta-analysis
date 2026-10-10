# Shared setup, sourced at the top of every stage script.
#
# Provides:
#   CODE_DIR, PARAMS                 repository root and config/params.yaml
#   INPUTS, RUN_DIR, RESULTS         workspace locations (see config/local.env.example)
#   stage_dir(name)                  output directory of a stage (created if needed)
#   read_meta(name)                  read a table from metadata/
#   write_table(df, path)            write a TSV and record it in the run's output log
#   save_object(obj, path)           saveRDS() and record it in the run's output log
#   check(), check_*()               stop with a clear message when an expectation fails

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
})

# --- locations --------------------------------------------------------------

# Stages are run from the repository root (see run.sh).
CODE_DIR <- normalizePath(".")
if (!file.exists(file.path(CODE_DIR, "config", "params.yaml"))) {
  stop("run the stages from the repository root (see run.sh)", call. = FALSE)
}

# config/local.env (KEY=VALUE lines) sets the workspace; environment
# variables of the same name take precedence.
local({
  env_file <- file.path(CODE_DIR, "config", "local.env")
  if (!file.exists(env_file)) return()
  for (line in readLines(env_file)) {
    line <- trimws(sub("#.*$", "", line))
    if (!nzchar(line) || !grepl("=", line, fixed = TRUE)) next
    key <- sub("=.*$", "", line)
    value <- gsub("^[\"']|[\"']$", "", sub("^[^=]*=", "", line))
    if (!nzchar(Sys.getenv(key))) do.call(Sys.setenv, setNames(list(value), key))
  }
})

WORKSPACE <- Sys.getenv("PLAVISCA_WORKSPACE")
if (!nzchar(WORKSPACE)) {
  stop("PLAVISCA_WORKSPACE is not set: copy config/local.env.example to config/local.env and edit it", call. = FALSE)
}
RUN_ID <- Sys.getenv("PLAVISCA_RUN", "main")
INPUTS <- file.path(WORKSPACE, "inputs")
RUN_DIR <- file.path(WORKSPACE, "runs", RUN_ID)
RESULTS <- file.path(RUN_DIR, "results")
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)

# PLAVISCA_PARAMS can name another parameter file (e.g. for a test run)
PARAMS_FILE <- Sys.getenv("PLAVISCA_PARAMS", file.path(CODE_DIR, "config", "params.yaml"))
PARAMS <- yaml::read_yaml(PARAMS_FILE)

input_path <- function(...) {
  path <- file.path(INPUTS, ...)
  check(file.exists(path), "input file missing: ", path)
  path
}

stage_dir <- function(name) {
  path <- file.path(RUN_DIR, "outs", name)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

results_dir <- function(name) {
  path <- file.path(RESULTS, name)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# --- checks -----------------------------------------------------------------

check <- function(ok, ...) {
  if (!isTRUE(all(ok))) stop("check failed: ", ..., call. = FALSE)
  invisible(TRUE)
}

check_unique <- function(x, what) {
  dups <- unique(x[duplicated(x)])
  check(length(dups) == 0, sprintf("%d duplicated %s (e.g. %s)", length(dups), what,
    paste(head(dups, 5), collapse = ", ")))
}

check_same_set <- function(a, b, what) {
  check(setequal(a, b), sprintf("%s differ: %d only in the first, %d only in the second (e.g. %s)",
    what, length(setdiff(a, b)), length(setdiff(b, a)),
    paste(head(c(setdiff(a, b), setdiff(b, a)), 5), collapse = ", ")))
}

check_same_order <- function(a, b, what) {
  check_same_set(a, b, what)
  check(identical(as.character(a), as.character(b)), what, " are the same but in a different order")
}

check_no_na <- function(x, what) {
  check(!anyNA(x), sprintf("%d missing values in %s", sum(is.na(x)), what))
}

# Cell counts are compared with metadata/expected_counts.tsv. A difference
# stops the pipeline: if it is intended, update that file and say why in the
# commit message.
check_expected_counts <- function(observed, check_name) {
  expected <- read_meta("expected_counts.tsv")
  expected <- expected[expected$check == check_name, , drop = FALSE]
  check(nrow(expected) > 0, "no rows for '", check_name, "' in metadata/expected_counts.tsv")
  key <- function(df) paste(df$study_label, ifelse(is.na(df$stage), "", df$stage), sep = " / ")
  merged <- merge(
    data.frame(key = key(expected), expected = as.integer(expected$cells)),
    data.frame(key = key(observed), observed = observed$cells),
    all = TRUE
  )
  merged[is.na(merged)] <- 0
  bad <- merged[merged$expected != merged$observed, , drop = FALSE]
  if (nrow(bad) > 0 && Sys.getenv("PLAVISCA_COUNT_CHECK") == "report") {
    # counts changed on purpose: list the differences and go on; once they
    # are reviewed, copy the new counts into metadata/expected_counts.tsv
    write_table(merged, file.path(stage_dir(STAGE), paste0("expected_counts_", check_name, ".tsv")))
    warning(sprintf("cell counts differ from metadata/expected_counts.tsv (%s) in %d row(s)", check_name, nrow(bad)), call. = FALSE)
    return(invisible(FALSE))
  }
  if (nrow(bad) > 0) {
    print(bad, row.names = FALSE)
    stop(sprintf("cell counts differ from metadata/expected_counts.tsv (%s) in %d row(s)", check_name, nrow(bad)), call. = FALSE)
  }
  message(sprintf("cell counts match metadata/expected_counts.tsv (%s, %d rows)", check_name, nrow(merged)))
  invisible(TRUE)
}

# --- tables and objects -----------------------------------------------------

# All columns are read as text; convert explicitly where a number or a
# TRUE/FALSE flag is needed.
read_meta <- function(name) {
  utils::read.delim(file.path(CODE_DIR, "metadata", name), stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = "", colClasses = "character")
}

# Every file written by a stage is recorded in <run>/outputs.tsv with its
# checksum and the code version that produced it.
record_output <- function(path) {
  log_path <- file.path(RUN_DIR, "outputs.tsv")
  commit <- tryCatch(system2("git", c("-C", CODE_DIR, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_, warning = function(w) NA_character_)
  dirty <- tryCatch(length(system2("git", c("-C", CODE_DIR, "status", "--porcelain", "--", "R", "config", "metadata"), stdout = TRUE)) > 0,
    error = function(e) NA)
  row <- data.frame(
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stage = STAGE,
    file = sub(paste0("^", RUN_DIR, "/"), "", normalizePath(path)),
    sha256 = as.character(openssl::sha256(file(path, raw = TRUE))),
    bytes = file.size(path),
    code_commit = if (length(commit)) commit[1] else NA,
    uncommitted_changes = dirty
  )
  utils::write.table(row, log_path, sep = "\t", quote = FALSE, row.names = FALSE,
    append = file.exists(log_path), col.names = !file.exists(log_path))
  invisible(path)
}

write_table <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  record_output(path)
}

save_object <- function(obj, path) {
  saveRDS(obj, path)
  record_output(path)
}

# Save the R session (package versions) of a stage next to its outputs.
save_session_info <- function(dir) {
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(dir, "session_info.txt"))
}

STUDIES <- read_meta("studies.tsv")$study_label

message(sprintf("[%s] workspace %s, run '%s'", if (exists("STAGE")) STAGE else "?", WORKSPACE, RUN_ID))
