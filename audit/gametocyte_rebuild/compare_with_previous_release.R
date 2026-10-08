#!/usr/bin/env Rscript
# Compare the rebuilt atlas with the previous (2026-09-21) release: cell
# membership (AT03 evidence) and harmonized_life_cycle_stage changes.
# Usage (repo root): Rscript audit/gametocyte_rebuild/compare_with_previous_release.R \
#   <previous pv_all_studies.rds> <new pv_all_studies.rds> <out_prefix>
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(SeuratObject))
args <- commandArgs(TRUE)
old <- readRDS(args[1])@meta.data
new <- readRDS(args[2])@meta.data
prefix <- args[3]
wt <- function(x, suffix) write.table(x, paste0(prefix, suffix), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

membership <- data.frame(
  check = c("total_cells_previous", "total_cells_new", "cell_id_set_equal", "cells_only_in_previous", "cells_only_in_new"),
  value = c(nrow(old), nrow(new), setequal(rownames(old), rownames(new)),
            length(setdiff(rownames(old), rownames(new))), length(setdiff(rownames(new), rownames(old))))
)
membership$status <- if (setequal(rownames(old), rownames(new)) && nrow(new) == 105949L) "PASS" else "FAIL"
wt(membership, "membership.tsv")
if (membership$status[1] != "PASS") stop("cell membership differs from the previous release")

new <- new[rownames(old), ]
lab <- function(x) ifelse(is.na(x), "NA", x)
o <- lab(old$harmonized_life_cycle_stage); n <- lab(new$harmonized_life_cycle_stage)
changed <- o != n
by_study <- as.data.frame.matrix(table(old$study_label, factor(changed, c(FALSE, TRUE), c("unchanged", "changed"))))
by_study <- cbind(study_label = rownames(by_study), by_study, row.names = NULL)
wt(by_study, "stage_changes_by_study.tsv")
transitions <- as.data.frame(table(study_label = old$study_label[changed], old_stage = o[changed], new_stage = n[changed]))
transitions <- transitions[transitions$Freq > 0, ]
names(transitions)[4] <- "n_cells"
wt(transitions[order(-transitions$n_cells), ], "stage_transitions.tsv")
totals <- as.data.frame(table(old_stage = o[changed], new_stage = n[changed]))
totals <- totals[totals$Freq > 0, ]
names(totals)[3] <- "n_cells"
wt(totals[order(-totals$n_cells), ], "stage_transitions_total.tsv")

options(width = 200)
print(membership, row.names = FALSE)
cat(sprintf("\ncells whose harmonized_life_cycle_stage changed: %d\n", sum(changed)))
print(by_study, row.names = FALSE)
print(totals[order(-totals$n_cells), ], row.names = FALSE)
