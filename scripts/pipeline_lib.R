# pipeline_lib.R
#
# Shared production-repair helpers for the PlaViSca pre-processing pipeline.
# Implements the fail-loud assertion discipline, controlled vocabularies, and
# deterministic mapping tables specified in
# audit/shared_pipeline/{metadata_schema,annotation_precedence,acceptance_tests}.tsv.
#
# Sourced by every per-study preprocessing script and by singleR.R/
# integration.R/flatten_data.R. Never itself reads or writes study data.

suppressPackageStartupMessages(library(dplyr))

# ---------------------------------------------------------------------------
# Environment defect workaround (Phase 2): rtracklayer's internal
# .parseSpeciesAsMetadata() has an upstream bug where a failed NCBI organism
# lookup returns NULL (via `warning()`, which itself returns NULL) instead
# of list(), which then crashes `metadata<-()` downstream in
# import.gff3(). This always triggers for PlasmoDB's GFF files, whose
# "##species http://.../wwwtax.cgi?id=9000000036" pragma uses a
# PlasmoDB-internal synthetic taxID that is never a real, resolvable NCBI
# taxonomy ID - so every import.gff3() call on scripts/ref/*.gff hits this
# path regardless of network availability. This patch only fixes the
# failure-path return value (list() instead of NULL); it does not change
# behavior when a species lookup succeeds. Applied once, here, so every
# per-study script that sources this file (all of which call import.gff3())
# is covered without a per-script workaround.
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  tryCatch(
    {
      ns <- asNamespace("rtracklayer")
      if (exists(".parseSpeciesAsMetadata", envir = ns, inherits = FALSE)) {
        patched_parseSpeciesAsMetadata <- function(lines) {
          species <- unique(grep("##species", lines, fixed = TRUE, value = TRUE))
          if (length(species) > 1L) stop("multiple species definitions found")
          metadata <- list()
          if (length(species) == 1L) {
            species <- sub("##species ", "", species, fixed = TRUE)
            if (get("isNCBISpeciesURL", envir = ns)(species)) {
              metadata <- tryCatch(
                get("metadataFromNCBI", envir = ns)(species),
                error = function(e) {
                  warning("failed to retrieve organism information from NCBI - continuing without species metadata")
                  list()
                }
              )
            }
          }
          metadata
        }
        environment(patched_parseSpeciesAsMetadata) <- ns
        unlockBinding(".parseSpeciesAsMetadata", ns)
        assign(".parseSpeciesAsMetadata", patched_parseSpeciesAsMetadata, envir = ns)
        lockBinding(".parseSpeciesAsMetadata", ns)
      }
    },
    error = function(e) {
      message("Note: could not apply rtracklayer .parseSpeciesAsMetadata workaround: ", conditionMessage(e))
    }
  )
}

# ---------------------------------------------------------------------------
# Controlled vocabularies (fixes D007/D056: one canonical string per study,
# never re-typed per script)
# ---------------------------------------------------------------------------

STUDY_LABELS <- c(
  mancio_silva2022 = "Mancio-Silva2022",
  sa2020 = "Sa2020",
  ruberto2022_1 = "Ruberto2022_1",
  ruberto2022_2 = "Ruberto2022_2",
  hazzard2022 = "Hazzard2022",
  hazzard2024 = "Hazzard2024"
)

# parasite_broad_stage is a deterministic, version-controlled image of
# harmonized_life_cycle_stage (metadata_schema.tsv row 41). It must never be
# independently overwritten - always derive it with map_broad_stage().
HARMONIZED_TO_BROAD_STAGE <- c(
  "EEF" = "Liver stage",
  "Hypnozoite" = "Liver stage",
  "Schizont (Liver stage)" = "Liver stage",
  "Ring" = "Blood stage",
  "Trophozoite" = "Blood stage",
  "Schizont (Blood stage)" = "Blood stage",
  "Merozoite" = "Blood stage",
  "Male gametocyte" = "Blood stage",
  "Female gametocyte" = "Blood stage",
  "Asexual (unresolved stage)" = "Blood stage",
  "Sporozoite" = "Sporozoite stage"
)

# tissue_or_sample_type canonical vocabulary (fixes the D058 vocabulary-drift
# root cause: "Host blood" vs "Mammalian host: blood" etc. must never diverge
# again). Every per-study script must emit exactly one of these values.
TISSUE_OR_SAMPLE_TYPE <- c(
  "Host blood",
  "Host liver",
  "Vector salivary gland",
  "Vector midgut"
)

# Populations that are never eligible for blood-stage IDC reference
# similarity classification (D058), regardless of expression magnitude.
NON_BLOOD_TISSUE_TYPES <- c("Host liver", "Vector salivary gland", "Vector midgut")

# ---------------------------------------------------------------------------
# Fail-loud assertion helpers (AT04, AT05: no vector-length mismatch, no
# positional recycling ever silently accepted)
# ---------------------------------------------------------------------------

#' Stop with a clearly labeled pipeline error (never a silent warning).
pipeline_fail <- function(msg) {
  stop(paste0("[PIPELINE ASSERTION FAILED] ", msg), call. = FALSE)
}

#' Assert that `x` has exactly `n` elements before it is used as a per-cell
#' metadata vector. Fixes D002,D027,D028,D033 (AT04).
assert_cardinality <- function(x, n, label = deparse(substitute(x))) {
  if (length(x) != n) {
    pipeline_fail(sprintf(
      "cardinality mismatch for '%s': length %d, expected %d (per-cell assignment must never rely on recycling)",
      label, length(x), n
    ))
  }
  invisible(TRUE)
}

#' Assert that every id in `keys` is present exactly once in `table_keys`
#' (no duplicate run keys, no unmatched runs, no unexpected extra runs).
#' Fixes AT05 and the "fail loudly" requirements in section 1 of the repair
#' spec.
assert_keyed_join <- function(keys, table_keys, label = "keys") {
  dup_keys <- table_keys[duplicated(table_keys)]
  if (length(dup_keys) > 0) {
    pipeline_fail(sprintf(
      "duplicate run keys in reference table for '%s': %s",
      label, paste(unique(dup_keys), collapse = ", ")
    ))
  }
  unmatched <- setdiff(keys, table_keys)
  if (length(unmatched) > 0) {
    pipeline_fail(sprintf(
      "unmatched runs for '%s': %s (present in data, absent from keyed reference table)",
      label, paste(unmatched, collapse = ", ")
    ))
  }
  extra <- setdiff(table_keys, keys)
  if (length(extra) > 0) {
    pipeline_fail(sprintf(
      "unexpected extra runs in keyed reference table for '%s': %s (present in table, absent from data)",
      label, paste(extra, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#' Assert no duplicated final cell IDs (AT01).
assert_unique_cell_ids <- function(ids, label = "cell IDs") {
  dups <- ids[duplicated(ids)]
  if (length(dups) > 0) {
    pipeline_fail(sprintf(
      "duplicated final cell IDs for '%s': %d duplicates (e.g. %s)",
      label, length(dups), paste(utils::head(unique(dups), 5), collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#' Assert an object has exactly `expected_n` cells (unexpected cell count).
assert_cell_count <- function(object_ncol, expected_n, label = "object") {
  if (object_ncol != expected_n) {
    pipeline_fail(sprintf(
      "unexpected cell count for '%s': got %d, expected %d",
      label, object_ncol, expected_n
    ))
  }
  invisible(TRUE)
}

#' Assert no mandatory source ID is missing (true NA disallowed for
#' mandatory identity fields).
assert_no_missing_mandatory <- function(x, label = deparse(substitute(x))) {
  if (anyNA(x) || any(x == "" & !is.na(x))) {
    pipeline_fail(sprintf(
      "missing mandatory source ID(s) in '%s': %d NA/empty of %d",
      label, sum(is.na(x) | x == ""), length(x)
    ))
  }
  invisible(TRUE)
}

#' Assert a character vector never contains the literal string "NA" (true NA
#' only). Fixes D024,D030,D062 (AT06).
assert_no_literal_na_string <- function(x, label = deparse(substitute(x))) {
  bad <- !is.na(x) & x == "NA"
  if (any(bad)) {
    pipeline_fail(sprintf(
      "literal string \"NA\" found in '%s' (%d cells) - use true NA (NA_character_), never the string \"NA\"",
      label, sum(bad)
    ))
  }
  invisible(TRUE)
}

#' Assert exact rowname-set equality between a metadata table and an
#' expression matrix before any join/cbind (AT16,AT17; fixes D053,D055).
assert_rowname_set_equal <- function(a_names, b_names, label_a = "a", label_b = "b") {
  if (!identical(sort(a_names), sort(b_names))) {
    only_a <- setdiff(a_names, b_names)
    only_b <- setdiff(b_names, a_names)
    pipeline_fail(sprintf(
      "rowname set mismatch between '%s' and '%s': %d only in %s, %d only in %s",
      label_a, label_b, length(only_a), label_a, length(only_b), label_b
    ))
  }
  invisible(TRUE)
}

#' Assert exact rowname order equality (stronger than set equality; used
#' immediately before a positional cbind, AT16).
assert_rowname_order_equal <- function(a_names, b_names, label_a = "a", label_b = "b") {
  assert_rowname_set_equal(a_names, b_names, label_a, label_b)
  if (!identical(a_names, b_names)) {
    pipeline_fail(sprintf(
      "rowname ORDER mismatch between '%s' and '%s' (sets match but order differs) - reorder explicitly before cbind",
      label_a, label_b
    ))
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Deterministic broad-stage derivation (fixes D020: makes "Sporozoite
# detailed stage + Blood broad stage" impossible by construction; AT09,AT21)
# ---------------------------------------------------------------------------

#' Deterministically derive parasite_broad_stage from
#' harmonized_life_cycle_stage. Never call this ad hoc with a hand-written
#' case_when - always route through this single mapping table so broad/
#' detailed can never silently diverge again.
map_broad_stage <- function(harmonized_life_cycle_stage) {
  unmapped <- setdiff(
    unique(harmonized_life_cycle_stage[!is.na(harmonized_life_cycle_stage)]),
    names(HARMONIZED_TO_BROAD_STAGE)
  )
  if (length(unmapped) > 0) {
    pipeline_fail(sprintf(
      "harmonized_life_cycle_stage value(s) not present in HARMONIZED_TO_BROAD_STAGE mapping table: %s",
      paste(unmapped, collapse = ", ")
    ))
  }
  unname(HARMONIZED_TO_BROAD_STAGE[harmonized_life_cycle_stage])
}

# ---------------------------------------------------------------------------
# Keyed run-table metadata application (fixes D002,D026-D028,D033: replaces
# `rep()`/positional-vector per-cell metadata assignment with an explicit
# keyed join against a committed run table, with cardinality/key assertions)
# ---------------------------------------------------------------------------

#' Split a list of per-run Seurat objects into those that ARE keyed in
#' run_table (survive to the final retained population) and those that are
#' NOT (e.g. Hazzard2024's imported-but-not-QC-retained runs, which a later
#' explicit subset() step in the merge script removes before saveRDS). This
#' makes explicit, at the earliest responsible step, exactly which runs are
#' authoritative vs transient/pre-filter - never a silent positional
#' assumption.
split_run_membership <- function(seurat_list, run_table) {
  run_ids <- names(seurat_list)
  list(
    retained = seurat_list[run_ids %in% run_table$run_id],
    pending_filter = seurat_list[!run_ids %in% run_table$run_id]
  )
}

#' Apply per-run metadata from a keyed run table to a list of per-run Seurat
#' objects, matching strictly on run_id. Fails loudly on any duplicate,
#' unmatched, or unexpected-extra run key (AT05).
#'
#' @param seurat_list named list of per-run Seurat objects, names = run_id
#' @param run_table data.frame with one row per run_id and the metadata
#'   columns to assign (must include column `run_id`)
#' @param fields character vector of column names in run_table to copy onto
#'   each Seurat object's metadata
apply_keyed_run_metadata <- function(seurat_list, run_table, fields) {
  run_ids <- names(seurat_list)
  assert_keyed_join(run_ids, run_table$run_id, label = "run_id")
  missing_fields <- setdiff(fields, colnames(run_table))
  if (length(missing_fields) > 0) {
    pipeline_fail(sprintf(
      "run_table is missing required field(s): %s",
      paste(missing_fields, collapse = ", ")
    ))
  }
  for (rid in run_ids) {
    row <- run_table[run_table$run_id == rid, , drop = FALSE]
    if (nrow(row) != 1) {
      pipeline_fail(sprintf("run_table has %d rows for run_id '%s', expected exactly 1", nrow(row), rid))
    }
    n_cells <- ncol(seurat_list[[rid]])
    for (f in fields) {
      seurat_list[[rid]][[f]] <- rep(row[[f]], n_cells)
    }
  }
  seurat_list
}

# ---------------------------------------------------------------------------
# Count provenance (Phase 2 Part I.4): every study must carry explicit,
# controlled-vocabulary provenance for its expression counts, not just
# Ruberto2022_1. Values are assigned once here so all six per-study scripts
# use the same field names/semantics.
# ---------------------------------------------------------------------------

#' Set the four count-provenance metadata fields on a Seurat object (applies
#' the same value to every current cell - call once per study object, after
#' cells are finalized, before saveRDS()).
set_count_provenance <- function(seurat_obj, source_of_counts, counting_pipeline,
                                  count_reference_version, count_provenance_status) {
  n_cells <- ncol(seurat_obj)
  # Direct @meta.data assignment (not `$<-`/`[[<-`), which is markedly
  # faster than Seurat's S4 setter path at tens-of-thousands-of-cells scale.
  seurat_obj@meta.data$source_of_counts <- rep(source_of_counts, n_cells)
  seurat_obj@meta.data$counting_pipeline <- rep(counting_pipeline, n_cells)
  seurat_obj@meta.data$count_reference_version <- rep(count_reference_version, n_cells)
  seurat_obj@meta.data$count_provenance_status <- rep(count_provenance_status, n_cells)
  seurat_obj
}

#' Set total_umi_count/near_empty_expression_flag from the adopted RNA
#' counts layer. singleR.R's IDC eligibility gate (Part 1) requires
#' total_umi_count to be present and non-NA for every cell in the merged
#' atlas, not just Ruberto2022_1 - every per-study script must call this
#' before saveRDS() so the field is never silently NA (and therefore
#' silently FALSE-like in a `>=` comparison) for the other five studies.
#' Assigns directly via @meta.data (bypasses the Seurat `$<-`/`[[<-` S4
#' dispatch path, which has been observed to be slow at tens-of-thousands-
#' of-cells scale) for speed on the largest studies.
set_qc_umi_fields <- function(seurat_obj, assay = "RNA", layer = "counts") {
  total_umi <- Matrix::colSums(SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer))
  seurat_obj@meta.data$total_umi_count <- unname(total_umi[colnames(seurat_obj)])
  seurat_obj@meta.data$near_empty_expression_flag <- seurat_obj@meta.data$total_umi_count <= 1
  seurat_obj
}

# ---------------------------------------------------------------------------
# Build manifest (Part 9 / AT23)
# ---------------------------------------------------------------------------

#' Record a build-manifest row for a produced artifact. Appends to
#' build_manifest.tsv at the repository root. Never overwrites prior rows.
record_build_manifest <- function(
  artifact_path,
  script_path,
  cell_count = NA_integer_,
  seed = NA_integer_,
  notes = "",
  manifest_path = "build_manifest.tsv"
) {
  sha256 <- function(path) {
    if (!file.exists(path)) {
      return(NA_character_)
    }
    tryCatch(
      as.character(openssl::sha256(file(path, raw = TRUE))),
      error = function(e) {
        tryCatch(
          digest::digest(file = path, algo = "sha256"),
          error = function(e2) NA_character_
        )
      }
    )
  }

  git_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  if (length(git_commit) == 0) git_commit <- NA_character_

  pkgs <- c("Seurat", "SeuratObject", "harmony", "SingleR")
  pkg_versions <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))

  row <- data.frame(
    date = as.character(Sys.time()),
    artifact_path = artifact_path,
    output_sha256 = sha256(artifact_path),
    script_path = script_path,
    git_commit = if (length(git_commit) > 0) git_commit[1] else NA_character_,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    Seurat_version = pkg_versions[["Seurat"]],
    SeuratObject_version = pkg_versions[["SeuratObject"]],
    harmony_version = pkg_versions[["harmony"]],
    SingleR_version = pkg_versions[["SingleR"]],
    cell_count = cell_count,
    seed = seed,
    notes = notes,
    stringsAsFactors = FALSE
  )

  if (file.exists(manifest_path)) {
    utils::write.table(
      row,
      manifest_path,
      append = TRUE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else {
    utils::write.table(
      row,
      manifest_path,
      append = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }
  invisible(row)
}

# ---------------------------------------------------------------------------
# AT31: stage-population floors. audit/shared_pipeline/stage_population_floors.tsv
# records, per population, the minimum number of cells that must carry a
# given harmonized_life_cycle_stage (optionally restricted to one
# source_stage_provenance; "any" = no restriction). Every source-defined
# population and the per-study gametocyte totals are listed, so a future
# precedence or gating change that drops one fails loudly instead of
# silently shrinking it.
# ---------------------------------------------------------------------------
STAGE_POPULATION_FLOORS_PATH <- "audit/shared_pipeline/stage_population_floors.tsv"

check_stage_population_floors <- function(meta, floors_path = STAGE_POPULATION_FLOORS_PATH) {
  if (!file.exists(floors_path)) {
    pipeline_fail(sprintf("Stage-population floor table %s is missing (AT31)", floors_path))
  }
  floors <- utils::read.delim(floors_path, stringsAsFactors = FALSE, check.names = FALSE)
  observed <- vapply(seq_len(nrow(floors)), function(i) {
    sel <- meta$study_label == floors$study_label[[i]] &
      meta$harmonized_life_cycle_stage %in% floors$harmonized_life_cycle_stage[[i]]
    if (floors$source_stage_provenance[[i]] != "any") {
      sel <- sel & meta$source_stage_provenance %in% floors$source_stage_provenance[[i]]
    }
    sum(sel, na.rm = TRUE)
  }, integer(1))
  floors$observed_cells <- observed
  floors$status <- ifelse(observed >= floors$min_cells, "PASS", "FAIL")
  floors
}

assert_stage_population_floors <- function(meta, floors_path = STAGE_POPULATION_FLOORS_PATH) {
  result <- check_stage_population_floors(meta, floors_path)
  failed <- result[result$status == "FAIL", , drop = FALSE]
  if (nrow(failed) > 0) {
    pipeline_fail(sprintf(
      "AT31 regression: %d stage population(s) fell below their recorded floor: %s",
      nrow(failed),
      paste(sprintf(
        "%s %s [%s] %d < %d", failed$study_label, failed$harmonized_life_cycle_stage,
        failed$source_stage_provenance, failed$observed_cells, failed$min_cells
      ), collapse = "; ")
    ))
  }
  invisible(result)
}

# ---------------------------------------------------------------------------
# rRNA features. The rRNA genes are the GFF gene records whose biotype is
# rRNA (ncRNA_gene with ebi_biotype=rRNA). Protein-coding genes whose
# *description* merely mentions rRNA (rRNA methyltransferases, pre-rRNA-
# processing proteins, fibrillarin, ...) are mRNAs and are NOT rRNA. The GFF
# is parsed directly: rtracklayer::import.gff3 queries NCBI for organism
# metadata and fails offline. Returned IDs use the Seurat feature spelling
# ("PVP01-..."), as in the per-study objects.
# ---------------------------------------------------------------------------
GFF_PATH <- "scripts/ref/PlasmoDB-68_PvivaxP01.gff"

gff_rrna_gene_ids <- function(gff_path = GFF_PATH) {
  rows <- strsplit(grep("^#", readLines(gff_path), value = TRUE, invert = TRUE), "\t", fixed = TRUE)
  rows <- rows[lengths(rows) == 9L]
  type <- vapply(rows, `[[`, character(1), 3L)
  attrs <- vapply(rows, `[[`, character(1), 9L)
  is_rrna_gene <- type == "ncRNA_gene" & grepl("(^|;)ebi_biotype=rRNA(;|$)", attrs)
  ids <- sub("^ID=", "", regmatches(attrs[is_rrna_gene], regexpr("ID=[^;]+", attrs[is_rrna_gene])))
  if (length(ids) == 0L || anyDuplicated(ids)) {
    pipeline_fail(sprintf("expected a non-empty, duplicate-free rRNA gene list from %s", gff_path))
  }
  gsub("_", "-", ids)
}

assert_no_rrna_features <- function(features, rrna_ids, label = "object") {
  hit <- intersect(gsub("_", "-", features), rrna_ids)
  if (length(hit) > 0) {
    pipeline_fail(sprintf("%s still contains %d rRNA feature(s): %s", label, length(hit), paste(utils::head(hit, 10), collapse = ", ")))
  }
  invisible(TRUE)
}
