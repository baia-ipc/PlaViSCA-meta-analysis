# Cell identifiers.
#
# A cell is named "<last 3 digits of its run accession>_<cell barcode>", e.g.
# "269_AAACCTGAGCTCAACT" for barcode AAACCTGAGCTCAACT of run SRR11008269.
# The 3-digit suffix identifies the run within a study only (e.g. 270 ends
# runs of both Sa2020 and Mancio-Silva2022), so lookups take the runs of one study.
# Mancio-Silva2022 keeps the authors' cell names, which already identify the
# sample (e.g. "D1SEQ100_GGTGCCCGGGGG").
#
# The functions below link these names to the cell names used by the authors
# in their own data, through the author sample names in metadata/runs.tsv.

cell_id <- function(run_id, barcode) paste0(substring(run_id, nchar(run_id) - 2), "_", barcode)

# runs: the rows of metadata/runs.tsv of one study
run_of_cell <- function(cell_ids, runs) {
  suffix <- sub("_.*$", "", cell_ids)
  run_suffix <- substring(runs$run_id, nchar(runs$run_id) - 2)
  check_unique(run_suffix, "run accession suffixes")
  run <- runs$run_id[match(suffix, run_suffix)]
  check_no_na(run, "run of each cell (cell name prefix not matching a run)")
  run
}

# Map author sample names (metadata/runs.tsv) to run accessions.
author_sample_runs <- function(runs) {
  runs <- runs[!is.na(runs$author_sample), ]
  check_unique(runs$author_sample, "author sample names")
  setNames(runs$run_id, runs$author_sample)
}

# Our names of author cells named "<author sample>_<barcode>"; the sample
# name may itself contain "_" (e.g. "NIH_Sa_AAAC...", "5309_2_AAAC...").
# what: the source of the names, for error messages
sample_barcode_cells <- function(author_ids, runs, what) {
  run <- author_sample_runs(runs)[sub("_[^_]+$", "", author_ids)]
  check_no_na(run, paste("run of each", what, "cell"))
  data.frame(cell_id = cell_id(run, sub("^.*_", "", author_ids)), author_cell_id = author_ids,
    stringsAsFactors = FALSE)
}

# Hazzard2024: the authors' processed table (Proccessed_Data.txt) has one row
# per retained cell, named "<author sample>_<barcode>", e.g. "5309_2_AAAC...".
hazzard2024_author_cells <- function(path, runs) {
  sample_barcode_cells(rownames(utils::read.delim(path)), runs, "Hazzard2024 author")
}

# Sa2020: metadata/sa2020_main_analysis_cells.tsv lists the cells of the
# authors' main analysis, e.g. "NIH_Sa_AAAC...".
sa2020_main_analysis_cells <- function(runs) {
  sample_barcode_cells(read_meta("sa2020_main_analysis_cells.tsv")$author_cell_id, runs, "Sa2020 main analysis")
}

# Ruberto2022_2: the authors' object combines sporozoites with blood stages;
# the sporozoites are its cells with Stage == "Sporozoite", named
# "<barcode>_<author sample>", e.g. "AAAC..._Case_909".
ruberto2022_2_author_cells <- function(path, runs) {
  author_object <- readRDS(path)
  stage <- author_object@meta.data$Stage
  author_ids <- rownames(author_object@meta.data)[!is.na(stage) & stage == "Sporozoite"]
  barcode <- sub("_.*$", "", author_ids)
  sample <- sub("^[^_]+_", "", author_ids)
  run <- author_sample_runs(runs)[sample]
  check_no_na(run, "run of each Ruberto2022_2 author cell")
  data.frame(cell_id = cell_id(run, barcode), author_cell_id = author_ids, stringsAsFactors = FALSE)
}
