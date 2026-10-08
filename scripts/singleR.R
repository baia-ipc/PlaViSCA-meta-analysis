library(DropletUtils) # scRNA-seq empty droplet detection
library(readxl) # data import
library(rtracklayer) # contains functions to handle gff/gtf files
library(Seurat) # contains scRNA-seq analysis functions
library(SingleCellExperiment) # single-cell data storage
library(SingleR) # single-cell data identification tool
library(scCustomize) # single-cell data handling and plotting tools
library(tidyverse)

# set current working directory
.plavisca_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PREPROCESS_ROOT", unset = NA_character_),
    getwd(),
    "/home/sopheap/pvsca_b/pre_process_data"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) {
    stop(
      "Could not locate the pre_process_data project root (looked for scripts/",
      "pipeline_lib.R under $PLAVISCA_PREPROCESS_ROOT, the current working ",
      "directory, and the historical hard-coded path). Set the ",
      "PLAVISCA_PREPROCESS_ROOT environment variable to the pre_process_data ",
      "directory's absolute path, or run this script from that directory.",
      call. = FALSE
    )
  }
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

pv.combined.all <- readRDS("./pv_all_studies.rds")

# ============================================================================
# PART 0: annotation-precedence bookkeeping (Part 4 / annotation_precedence.tsv)
# ============================================================================
# source_life_cycle_stage / source_stage_provenance are set upstream, at the
# earliest responsible per-study preprocessing step (silva.R,
# ruberto2022_1/2_pv_analysis_script.R, hazzard2022_pv_analysis_script.R),
# and are NEVER overwritten here. Every study object must already carry
# tissue_or_sample_type using the single canonical TISSUE_OR_SAMPLE_TYPE
# vocabulary (fixes the D058 vocabulary-drift root cause) before it reaches
# this script.
required_provenance_fields <- c("tissue_or_sample_type", "total_umi_count")
missing_fields <- setdiff(required_provenance_fields, colnames(pv.combined.all@meta.data))
if (length(missing_fields) > 0) {
  pipeline_fail(sprintf(
    "pv_all_studies.rds is missing required provenance field(s) %s - every per-study script must set these before the merged atlas reaches singleR.R",
    paste(missing_fields, collapse = ", ")
  ))
}

unexpected_tissue_values <- setdiff(
  unique(pv.combined.all$tissue_or_sample_type[!is.na(pv.combined.all$tissue_or_sample_type)]),
  TISSUE_OR_SAMPLE_TYPE
)
if (length(unexpected_tissue_values) > 0) {
  pipeline_fail(sprintf(
    "tissue_or_sample_type contains value(s) outside the canonical TISSUE_OR_SAMPLE_TYPE vocabulary: %s (fixes the D058 vocabulary-drift root cause - add the value to TISSUE_OR_SAMPLE_TYPE in pipeline_lib.R only after confirming it is a genuinely new population, never as a silent workaround)",
    paste(unexpected_tissue_values, collapse = ", ")
  ))
}

if (!"source_life_cycle_stage" %in% colnames(pv.combined.all@meta.data)) {
  pv.combined.all$source_life_cycle_stage <- NA_character_
}
if (!"source_stage_provenance" %in% colnames(pv.combined.all@meta.data)) {
  pv.combined.all$source_stage_provenance <- NA_character_
}

# ============================================================================
# PART 1: eligibility gating BEFORE classification (fixes D058, root cause of
# D021,D034,D037,D044) - blood-IDC reference similarity is computed ONLY for
# biologically eligible blood-stage cells with usable expression. Liver
# stages, sporozoites, mosquito-stage populations, source-defined sexual
# populations, and zero/near-zero expression cells are excluded up front,
# never classified-then-masked.
# ============================================================================
MIN_UMI_FOR_IDC_INFERENCE <- 10 # documented, explicit minimum (AT11)

idc_eligible <- (
  pv.combined.all$tissue_or_sample_type == "Host blood" &
    !is.na(pv.combined.all$tissue_or_sample_type) &
    (is.na(pv.combined.all$source_stage_provenance) |
      pv.combined.all$source_stage_provenance != "source_selection_defined") &
    pv.combined.all$total_umi_count >= MIN_UMI_FOR_IDC_INFERENCE
)

cat(sprintf(
  "[singleR.R] IDC eligibility gate: %d/%d cells eligible for blood-IDC reference similarity classification (%d excluded: non-blood tissue, source-selection-defined stage, or <%d total UMI)\n",
  sum(idc_eligible), ncol(pv.combined.all), sum(!idc_eligible), MIN_UMI_FOR_IDC_INFERENCE
))

# ============================================================================
# PART 2: Zhu SMRU1 -> PvP01 IDC reference construction (fixes D059: logged
# drop_na count, and the reference is explicitly documented as a bulk,
# single-replicate, cross-platform microarray time course, never presented
# as scRNA-seq-comparable counts)
# ============================================================================
Zhu_HiDef_PVX_IDC_timecourse_1 <- read_excel(
  "./data/Zhu_SciReps_2016.xls",
  sheet = "SMRU1"
)

PVX_PVP01_orthologs <- read_excel("./data/GenesByOrthologs_Summary.xlsx")

Zhu_HiDef_PVX_IDC_timecourse_1 <- as.data.frame(Zhu_HiDef_PVX_IDC_timecourse_1)
n_zhu_genes_before_join <- nrow(Zhu_HiDef_PVX_IDC_timecourse_1)

Zhu_HiDef_PVX_IDC_timecourse_1 <- left_join(
  Zhu_HiDef_PVX_IDC_timecourse_1,
  PVX_PVP01_orthologs,
  by = c("Gid" = "InputOrtholog")
)
n_zhu_genes_after_join_before_drop <- nrow(Zhu_HiDef_PVX_IDC_timecourse_1)
Zhu_HiDef_PVX_IDC_timecourse_1 <- drop_na(Zhu_HiDef_PVX_IDC_timecourse_1)
n_zhu_genes_dropped <- n_zhu_genes_after_join_before_drop - nrow(Zhu_HiDef_PVX_IDC_timecourse_1)

cat(sprintf(
  "[singleR.R] Zhu SMRU1 reference: %d genes before ortholog join, %d dropped by drop_na() (unmapped/incomplete orthologs), %d retained in the final reference\n",
  n_zhu_genes_before_join, n_zhu_genes_dropped, nrow(Zhu_HiDef_PVX_IDC_timecourse_1)
))

Zhu_HiDef_PVX_IDC_timecourse_1 <- Zhu_HiDef_PVX_IDC_timecourse_1[, c(2:9)]
rownames(
  Zhu_HiDef_PVX_IDC_timecourse_1
) <- Zhu_HiDef_PVX_IDC_timecourse_1$GeneID
Zhu_HiDef_PVX_IDC_timecourse_1 <- Zhu_HiDef_PVX_IDC_timecourse_1[, c(1:7)]
counts <- as.matrix(Zhu_HiDef_PVX_IDC_timecourse_1)

# D059: assay explicitly documented/named as a bulk, single-replicate,
# log2-scale microarray reference - NOT scRNA-seq counts, despite the
# assay-slot name required by SummarizedExperiment/SingleR conventions.
IDC_cycle <- SummarizedExperiment(
  assay = SimpleList(logcounts = counts),
  colData = colnames(Zhu_HiDef_PVX_IDC_timecourse_1)
)
IDC_cycle$X <- NULL
metadata(IDC_cycle)$reference_provenance <- paste(
  "Zhu et al. 2016 Sci Rep SMRU1 IDC time course:",
  "bulk, single-replicate, log2-scale microarray intensities across 7",
  "timepoints (6,12,24,32,36,42,48h), converted PVX->PvP01 via a 1:1",
  "ortholog join. NOT scRNA-seq-comparable counts."
)

cat(sprintf(
  "[singleR.R] IDC reference dimensions: %d genes x %d timepoints\n",
  nrow(IDC_cycle), ncol(IDC_cycle)
))

# ============================================================================
# PART 3: SingleR classification, restricted to eligible cells only (fixes
# D058's core defect: gating happens BEFORE, not after, classification)
# ============================================================================
sce.eligible <- as.SingleCellExperiment(DietSeurat(
  subset(pv.combined.all, cells = colnames(pv.combined.all)[idc_eligible]),
  layers = "data",
  assays = "RNA"
))
rownames(sce.eligible) <- gsub("_", "-", rownames(sce.eligible))

gene_intersection <- intersect(rownames(sce.eligible), rownames(IDC_cycle))
cat(sprintf(
  "[singleR.R] Test/reference gene intersection: %d genes (of %d reference genes, %d test genes)\n",
  length(gene_intersection), nrow(IDC_cycle), nrow(sce.eligible)
))
cat(sprintf(
  "[singleR.R] Eligible cells classified: %d; ineligible cells (gated out, NA by construction): %d\n",
  ncol(sce.eligible), ncol(pv.combined.all) - ncol(sce.eligible)
))

pred.IDC.eligible <- SingleR(
  test = sce.eligible,
  ref = IDC_cycle,
  assay.type.test = 1,
  labels = colnames(IDC_cycle)
)

# ---------------------------------------------------------------------------
# D058: renamed from the misleading "hour_post_invasion" to
# idc_reference_similarity_label - these are discrete reference-profile
# labels (which of the 7 Zhu timepoints a cell's expression best correlates
# with), never a measured infection time. NA by construction for every
# ineligible cell (never computed, never masked post hoc).
# D061: raw (unpruned) label, delta.next, and a pruned flag are all
# persisted, not just the final pruned label.
# ---------------------------------------------------------------------------
pv.combined.all$idc_reference_similarity_label <- NA_character_
pv.combined.all$idc_similarity_raw_label <- NA_character_
pv.combined.all$idc_similarity_score_delta <- NA_real_
pv.combined.all$idc_similarity_pruned_flag <- NA

eligible_cells <- colnames(sce.eligible)
pv.combined.all$idc_reference_similarity_label[eligible_cells] <- pred.IDC.eligible$pruned.labels
pv.combined.all$idc_similarity_raw_label[eligible_cells] <- pred.IDC.eligible$labels
pv.combined.all$idc_similarity_score_delta[eligible_cells] <- pred.IDC.eligible$delta.next
pv.combined.all$idc_similarity_pruned_flag[eligible_cells] <-
  is.na(pred.IDC.eligible$pruned.labels) & !is.na(pred.IDC.eligible$labels)

# ============================================================================
# PART 4: heuristic stage bins (fixes D060: explicitly labeled as a
# PlaViSca-defined heuristic, never presented as source-author annotation)
# ============================================================================
idc_hour_numeric <- suppressWarnings(as.numeric(pv.combined.all$idc_reference_similarity_label))
pv.combined.all$idc_heuristic_stage_bin <- case_when(
  is.na(idc_hour_numeric) ~ NA_character_,
  idc_hour_numeric >= 0 & idc_hour_numeric < 18 ~ "Ring",
  idc_hour_numeric >= 18 & idc_hour_numeric < 30 ~ "Trophozoite",
  # Phase 2 fix: this heuristic bin feeds directly into
  # harmonized_life_cycle_stage (Part 5a below), which must use the
  # blood/liver-disambiguated vocabulary in HARMONIZED_TO_BROAD_STAGE
  # (pipeline_lib.R) - a bare "Schizont" is not a key in that table and
  # would crash map_broad_stage() the first time a cell actually lands in
  # this bin (never previously exercised end-to-end before the Phase 2
  # atlas rebuild).
  idc_hour_numeric >= 30 & idc_hour_numeric < 46 ~ "Schizont (Blood stage)",
  idc_hour_numeric >= 46 & idc_hour_numeric <= 48 ~ "Merozoite",
  TRUE ~ NA_character_
)
# These cutoffs (0-<18 Ring, 18-<30 Trophozoite, 30-<46 Schizont, 46-48
# Merozoite) are PlaViSca-defined heuristics with no citation in Zhu et al.
# 2016 or elsewhere - never present idc_heuristic_stage_bin as an
# author-established or literature-sourced boundary (D060).

plotDeltaDistribution(pred.IDC.eligible, ncol = 3)

# ============================================================================
# PART 5: deterministic annotation precedence (Part 4 of the repair spec;
# annotation_precedence.tsv) - source stage always wins for population-type
# gating over any inference; harmonized_life_cycle_stage is assembled
# explicitly, once, from documented per-population precedence rules.
# ============================================================================
harmonized <- rep(NA_character_, ncol(pv.combined.all))

# 5a. Start from the IDC heuristic bin wherever a cell was eligible and
# received one (blood-stage asexual populations).
harmonized[!is.na(pv.combined.all$idc_heuristic_stage_bin)] <-
  pv.combined.all$idc_heuristic_stage_bin[!is.na(pv.combined.all$idc_heuristic_stage_bin)]

# 5b. Historical refine_state (Mancio-Silva2022 topcell-derived detailed
# label) wins over the IDC bin where present - preserved as its own
# source_stage_provenance value, not silently merged.
if ("refine_state" %in% colnames(pv.combined.all@meta.data)) {
  has_refine_state <- !is.na(pv.combined.all$refine_state)
  harmonized[has_refine_state] <- pv.combined.all$refine_state[has_refine_state]
  pv.combined.all$source_stage_provenance[has_refine_state] <- ifelse(
    is.na(pv.combined.all$source_stage_provenance[has_refine_state]),
    "source_state_defined",
    pv.combined.all$source_stage_provenance[has_refine_state]
  )
}

# 5c. Ruberto2022_1 liver-form source labels (Hypnozoite/Schizont) win for
# life_cycle_stage; recorded as "Schizont (Liver stage)" per the target
# vocabulary. Liver-stage cells never receive pred_gametocyte/IDC labels -
# enforced structurally below in Part 7, not just by this precedence step.
if ("liver_form" %in% colnames(pv.combined.all@meta.data)) {
  liver_hypnozoite <- pv.combined.all$liver_form %in% c("Hypnozoites", "Hypnozoite")
  liver_schizont <- pv.combined.all$liver_form %in% c("Schizonts", "Schizont")
  harmonized[liver_hypnozoite] <- "Hypnozoite"
  harmonized[liver_schizont] <- "Schizont (Liver stage)"
  is_liver <- liver_hypnozoite | liver_schizont
  pv.combined.all$source_stage_provenance[is_liver] <- "source_selection_defined"
}

# 5d. Source-selection/source-run-defined Sporozoite populations
# (Ruberto2022_2: all cells; Hazzard2022: runs 498/499) are authoritative
# and win over EVERYTHING, including any prior step above. This sets BOTH
# the detailed and (via map_broad_stage() below) the broad field atomically,
# fixing D020's Sporozoite-detailed/Blood-broad contradiction by
# construction, and is applied last among the "source wins" rules so a
# downstream inference step can never re-open it (see Part 7).
source_defined_sporozoite <- (
  !is.na(pv.combined.all$source_stage_provenance) &
    pv.combined.all$source_stage_provenance == "source_selection_defined" &
    !is.na(pv.combined.all$source_life_cycle_stage) &
    pv.combined.all$source_life_cycle_stage == "Sporozoite"
)
harmonized[source_defined_sporozoite] <- "Sporozoite"

pv.combined.all$harmonized_life_cycle_stage <- harmonized
pv.combined.all$source_stage_provenance[is.na(pv.combined.all$source_stage_provenance)] <- "singleR_inferred"

# ============================================================================
# PART 6: de novo clustering (fixes D047: use the real Seurat 5.3.0
# argument `graph.name`, not the nonexistent `graph_name`, and explicitly
# assert the intended graph exists rather than silently falling back)
# ============================================================================
intended_graph_name <- "RNA_snn"
if (!intended_graph_name %in% names(pv.combined.all@graphs)) {
  pipeline_fail(sprintf(
    "Intended clustering graph '%s' does not exist on this object (available: %s). integration.R must build this graph via FindNeighbors before singleR.R runs FindClusters.",
    intended_graph_name, paste(names(pv.combined.all@graphs), collapse = ", ")
  ))
}

pv.combined.all <- FindClusters(
  pv.combined.all,
  resolution = 0.05,
  graph.name = intended_graph_name
)

# ============================================================================
# PART 7: gametocyte marker module scores + marker-driven cluster identity
# (fixes D048: no hard-coded numeric cluster->sex mapping; DEC12 option (b)
# interim safeguard - a mandatory, versioned marker-based verification step
# runs before any cluster is ever named Female/Male gametocyte, and the
# pipeline fails loudly rather than silently assigning sex from an
# unvalidated cluster)
# ============================================================================
femaleGams <- c(
  "PVP01-1207200", "PVP01-0616100", "PVP01-1119300", "PVP01-1465500",
  "PVP01-1259400", "PVP01-1441000", "PVP01-0702600", "PVP01-1306800",
  "PVP01-0946800", "PVP01-0517400", "PVP01-1027600", "PVP01-1003000",
  "PVP01-1143200", "PVP01-1017500", "PVP01-1024300", "PVP01-0806000",
  "PVP01-1240300", "PVP01-0603600", "PVP01-0712800"
)

maleGams <- c(
  "PVP01-1262200", "PVP01-0530800", "PVP01-1412100", "PVP01-1025600",
  "PVP01-1266500", "PVP01-1229400"
)

pv.combined.all <- AddModuleScore(
  object = pv.combined.all,
  features = list(c(femaleGams)),
  name = "femaleGams_module_score"
)
pv.combined.all <- AddModuleScore(
  object = pv.combined.all,
  features = list(c(maleGams)),
  name = "maleGams_module_score"
)

# Cells not eligible for blood-stage gametocyte inference (non-blood tissue,
# source-selection-defined stage, or near-zero expression) are excluded from
# marker-based cluster scoring, mirroring the IDC eligibility gate (D037).
gametocyte_eligible <- idc_eligible

cluster_marker_summary <- pv.combined.all@meta.data %>%
  rownames_to_column("cell_id") %>%
  filter(gametocyte_eligible[cell_id]) %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    mean_female_score = mean(femaleGams_module_score1, na.rm = TRUE),
    mean_male_score = mean(maleGams_module_score1, na.rm = TRUE),
    .groups = "drop"
  )

# Versioned cluster/population marker summary (DEC12 (b): the deterministic
# evidence a human reviewer or a future automated check validates before
# trusting any cluster->sex mapping - written next to the object it
# describes, never silently trusted).
write_tsv(cluster_marker_summary, "cluster_sex_marker_summary.tsv")

# A cluster is called Female/Male gametocyte ONLY if its mean module score
# for one sex clears an explicit, documented threshold AND exceeds the other
# sex's mean score by an explicit margin. Any cluster that does not clear
# both bars is left unassigned - never silently defaulted to a numeric
# cluster identity. If NO cluster clears the bar in either direction, the
# pipeline fails loudly rather than silently proceeding with an unvalidated
# empty mapping.
FEMALE_SCORE_THRESHOLD <- 0.10
MALE_SCORE_THRESHOLD <- 0.10
SCORE_MARGIN <- 0.05

cluster_marker_summary <- cluster_marker_summary %>%
  mutate(
    call = case_when(
      mean_female_score >= FEMALE_SCORE_THRESHOLD &
        (mean_female_score - mean_male_score) >= SCORE_MARGIN ~ "Female gametocyte",
      mean_male_score >= MALE_SCORE_THRESHOLD &
        (mean_male_score - mean_female_score) >= SCORE_MARGIN ~ "Male gametocyte",
      TRUE ~ "Asexual"
    )
  )

if (!any(cluster_marker_summary$call %in% c("Female gametocyte", "Male gametocyte"))) {
  pipeline_fail(paste(
    "No cluster cleared the marker-score threshold for either Female or Male",
    "gametocyte identity (see cluster_sex_marker_summary.tsv). Refusing to",
    "silently proceed with an unvalidated cluster->sex mapping (D048/DEC12).",
    "Either the clustering resolution needs adjustment, the thresholds need",
    "reviewed adjustment, or gametocytes are genuinely absent/undetectable",
    "in this build - a human must review cluster_sex_marker_summary.tsv",
    "before proceeding."
  ))
}

cluster_to_label <- setNames(cluster_marker_summary$call, cluster_marker_summary$seurat_clusters)
pv.combined.all$pred_gametocyte_sex <- unname(cluster_to_label[as.character(pv.combined.all$seurat_clusters)])
pv.combined.all$pred_gametocyte_sex[pv.combined.all$pred_gametocyte_sex == "Asexual"] <- NA_character_

# D037,AT10,AT11: liver-stage cells and zero/near-zero-expression cells can
# never carry a gametocyte-sex call, regardless of what cluster they fall
# into - this is a structural gate, not a downstream cleanup step.
not_gametocyte_eligible <- !gametocyte_eligible
pv.combined.all$pred_gametocyte_sex[not_gametocyte_eligible] <- NA_character_

# ============================================================================
# PART 8: cluster-level gametocyte call enters harmonized_life_cycle_stage
# (precedence tier 2 of annotation_precedence.tsv: source-defined label >
# marker/cluster gametocyte call > IDC heuristic bin). The call replaces the
# IDC heuristic bin for gametocyte-eligible cells only; source-defined cells
# (source_selection_defined: sporozoites, Ruberto2022_1 liver forms;
# source_state_defined: Mancio-Silva2022 refine_state) are never overwritten.
# idc_heuristic_stage_bin is kept as its own column, so the asexual
# classification of these cells remains inspectable. Sa2020 source labels
# (tier 1) are applied on top of this in Part 9.
# ============================================================================
male_female_call <- !is.na(pv.combined.all$pred_gametocyte_sex)
source_defined_stage <- pv.combined.all$source_stage_provenance %in%
  c("source_selection_defined", "source_state_defined")
overwrite_target <- male_female_call & gametocyte_eligible & !source_defined_stage
pv.combined.all$harmonized_life_cycle_stage[overwrite_target] <-
  pv.combined.all$pred_gametocyte_sex[overwrite_target]
pv.combined.all$source_stage_provenance[overwrite_target] <- "marker_cluster_inferred"
cat(sprintf(
  "[singleR.R] Cluster-level gametocyte call written to harmonized_life_cycle_stage for %d cell(s) (%d Female, %d Male)\n",
  sum(overwrite_target),
  sum(pv.combined.all$pred_gametocyte_sex[overwrite_target] == "Female gametocyte"),
  sum(pv.combined.all$pred_gametocyte_sex[overwrite_target] == "Male gametocyte")
))

# Hazzard2024 source_pseudogroup ("Sexual" carries no sex and no stage) is a
# source field only and never enters harmonized_life_cycle_stage; its
# agreement with the cluster call is recorded in the build log.
hazzard2024_cells <- pv.combined.all$study_label == "Hazzard2024"
cat("[singleR.R] Hazzard2024 source_pseudogroup x pred_gametocyte_sex (cluster call):\n")
print(table(
  source_pseudogroup = pv.combined.all$source_pseudogroup[hazzard2024_cells],
  cluster = pv.combined.all$pred_gametocyte_sex[hazzard2024_cells],
  useNA = "ifany"
))

violation <- male_female_call &
  !is.na(pv.combined.all$source_stage_provenance) &
  pv.combined.all$source_stage_provenance == "source_selection_defined"
if (any(violation)) {
  pipeline_fail(sprintf(
    "%d source-selection-defined cell(s) received a gametocyte-sex call - this must be structurally impossible (D045/D048/D049). Cell IDs: %s",
    sum(violation), paste(utils::head(colnames(pv.combined.all)[violation], 10), collapse = ", ")
  ))
}

# ============================================================================
# PART 9: read gametocyte metadata from Sa et al. and reconcile with
# life_cycle_stage (fixes D014,D015: exact set-membership match instead of
# the buggy `filter(type == c(...))` parity recycling, and the authoritative
# NIH/PB_MACS prefix crosswalk instead of the wrong-cardinality mapping that
# collided NIH_Ao/NIH_CQ onto the same run suffix and dropped PB_MACS
# entirely)
# ============================================================================
sa2020_source_sex <- read_excel("./data/pbio.3000711.s034.xlsx") |>
  select(sample = Sample, type = Type) |>
  mutate(
    prefix = sub("(_[^_]+)$", "", sample),
    cells = sub("^[^_]+_[^_]+_", "", sample)
  ) |>
  select(prefix, cells, type) |>
  mutate(
    # D015: authoritative crosswalk (validated against Table 1 + ENA +
    # barcode overlap) - NIH_Sa=269, NIH_CQ=270, NIH_Ao=271, PB_MACS=272.
    # The historical version collided NIH_Ao/NIH_CQ onto the same suffix
    # (269) and never mapped PB_MACS at all.
    prefix = case_when(
      prefix == "AMRU_Ao" ~ "278",
      prefix == "AMRU_CQ" ~ "277",
      prefix == "AMRU_Sa" ~ "276",
      prefix == "Ches_Ao" ~ "275",
      prefix == "Ches_CQ" ~ "274",
      prefix == "Ches_Sa" ~ "273",
      prefix == "NIH_Sa" ~ "269",
      prefix == "NIH_CQ" ~ "270",
      prefix == "NIH_Ao" ~ "271",
      prefix == "PB_MACS" ~ "272",
      TRUE ~ NA_character_
    ),
    new_cells = paste0(prefix, "_", cells),
    type = case_when(
      type == "Female" ~ "Female gametocyte",
      type == "Male" ~ "Male gametocyte"
    )
  ) |>
  # D014: exact set-membership match, not the parity-recycled
  # `filter(type == c("Male gametocyte","Female gametocyte"))`.
  filter(type %in% c("Male gametocyte", "Female gametocyte"))

pv.combined.all$source_sex_annotation <- NA_character_
sa2020_match_cell <- intersect(pv.combined.all$barcode, sa2020_source_sex$new_cells)
cat(sprintf(
  "[singleR.R] Sa2020 source sex-annotation import: %d matched cells (of %d source-labeled rows)\n",
  length(sa2020_match_cell), nrow(sa2020_source_sex)
))

sa2020_barcode_idx <- match(pv.combined.all$barcode, sa2020_source_sex$new_cells)
has_sa2020_match <- !is.na(sa2020_barcode_idx)
pv.combined.all$source_sex_annotation[has_sa2020_match] <-
  sa2020_source_sex$type[sa2020_barcode_idx[has_sa2020_match]]

# Source labels are only ever imported into Sa2020 host-blood cells; a match
# anywhere else would be a crosswalk error (AT13), never a label to apply.
bad_match <- has_sa2020_match & (
  pv.combined.all$study_label != "Sa2020" |
    pv.combined.all$tissue_or_sample_type != "Host blood" |
    pv.combined.all$source_stage_provenance == "source_selection_defined"
)
if (any(bad_match)) {
  pipeline_fail(sprintf(
    "%d Sa2020 source-label match(es) fall outside Sa2020 host-blood cells (AT13 crosswalk regression)",
    sum(bad_match)
  ))
}

# Build-log record of Sa2020 source label vs cluster call, before the
# source label takes precedence.
sa2020_cells <- pv.combined.all$study_label == "Sa2020"
cat("[singleR.R] Sa2020 source_sex_annotation x pred_gametocyte_sex (cluster call):\n")
print(table(
  source = pv.combined.all$source_sex_annotation[sa2020_cells],
  cluster = pv.combined.all$pred_gametocyte_sex[sa2020_cells],
  useNA = "ifany"
))
n_opposite_sex <- sum(has_sa2020_match &
  !is.na(pv.combined.all$pred_gametocyte_sex) &
  pv.combined.all$pred_gametocyte_sex != pv.combined.all$source_sex_annotation)
n_source_no_call <- sum(has_sa2020_match & is.na(pv.combined.all$pred_gametocyte_sex))
cat(sprintf(
  "[singleR.R] Sa2020 source/cluster disagreements resolved in favour of the source label: %d opposite-sex cluster call, %d source gametocyte without cluster call (total %d of %d source-labelled cells)\n",
  n_opposite_sex, n_source_no_call, n_opposite_sex + n_source_no_call, sum(has_sa2020_match)
))

# Precedence tier 1 (annotation_precedence.tsv): the Sa2020 author labels
# are source-defined and win over the cluster call for those cells only.
pv.combined.all$harmonized_life_cycle_stage[has_sa2020_match] <-
  pv.combined.all$source_sex_annotation[has_sa2020_match]
pv.combined.all$source_stage_provenance[has_sa2020_match] <- "source_state_defined"

# adjudicated_display_stage follows the same deterministic precedence
# (source-defined > cluster call > IDC bin) and is therefore identical to
# harmonized_life_cycle_stage.
pv.combined.all$adjudicated_display_stage <- pv.combined.all$harmonized_life_cycle_stage

# ============================================================================
# PART 10: broad/detailed stage consistency (fixes D020 by construction) -
# parasite_broad_stage is deterministically derived from
# harmonized_life_cycle_stage via the single version-controlled mapping
# table in pipeline_lib.R, and must never be independently set anywhere
# else in this script.
# ============================================================================
pv.combined.all$parasite_broad_stage <- map_broad_stage(pv.combined.all$harmonized_life_cycle_stage)

# AT09/AT21 self-check: source-defined sporozoites must be broad-labeled
# Sporozoite stage, never Blood stage (the exact D020 contradiction).
sporozoite_broad_check <- pv.combined.all$parasite_broad_stage[source_defined_sporozoite]
if (length(sporozoite_broad_check) > 0 && any(sporozoite_broad_check != "Sporozoite stage", na.rm = TRUE)) {
  pipeline_fail("D020 regression: at least one source-selection-defined Sporozoite cell has parasite_broad_stage != 'Sporozoite stage'")
}

# ============================================================================
# PART 11: post-integration removal policy (fixes D045+D049 together,
# DEC10 option (b)) - gate on source-defined population type, never solely
# on downstream classifier output, and log every removed cell as part of
# build provenance. Because source-selection-defined Sporozoite cells can no
# longer receive a gametocyte-sex call (Part 8's structural assertion), the
# expected removal count is zero unless a separately approved scientific
# decision changes this.
# ============================================================================
removal_candidates <- (
  pv.combined.all$harmonized_life_cycle_stage %in% c("Male gametocyte", "Female gametocyte") &
    grepl("Anopheles", pv.combined.all$host_species) &
    (is.na(pv.combined.all$source_stage_provenance) |
      pv.combined.all$source_stage_provenance != "source_selection_defined")
)

n_removed <- sum(removal_candidates)
if (n_removed > 0) {
  removal_manifest <- pv.combined.all@meta.data[removal_candidates, , drop = FALSE] %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id,
      study_label,
      run_id,
      host_species,
      harmonized_life_cycle_stage,
      reason = "Anopheles-host cell classified as gametocyte by downstream classifier, not source-selection-defined - removed per DEC10 policy"
    )
  write_tsv(removal_manifest, "post_integration_removal_manifest.tsv")
  cat(sprintf("[singleR.R] Removing %d cell(s) per the post-integration Anopheles+gametocyte policy - see post_integration_removal_manifest.tsv\n", n_removed))
} else {
  cat("[singleR.R] Post-integration Anopheles+gametocyte removal filter: 0 cells removed (expected under the current build)\n")
}

pv.combined.all <- subset(pv.combined.all, cells = colnames(pv.combined.all)[!removal_candidates])

# ============================================================================
# Final assertions (AT07,AT08,AT09,AT10,AT21,AT22)
# ============================================================================
assert_unique_cell_ids(colnames(pv.combined.all), label = "pv_all_studies.rds (post-singleR.R)")

non_blood_idc <- pv.combined.all$tissue_or_sample_type != "Host blood" &
  !is.na(pv.combined.all$tissue_or_sample_type)
if (any(!is.na(pv.combined.all$idc_reference_similarity_label[non_blood_idc]))) {
  pipeline_fail("AT08 regression: idc_reference_similarity_label is non-NA for at least one non-blood-stage cell")
}

liver_gametocyte <- pv.combined.all$parasite_broad_stage == "Liver stage" &
  !is.na(pv.combined.all$parasite_broad_stage)
if (any(!is.na(pv.combined.all$pred_gametocyte_sex[liver_gametocyte]))) {
  pipeline_fail("AT10 regression: pred_gametocyte_sex is non-NA for at least one Liver stage cell")
}

# A gametocyte harmonized_life_cycle_stage may only come from the Sa2020
# source labels (source_state_defined, Host blood) or the cluster call
# (marker_cluster_inferred, gametocyte-eligible cells) - never on a liver,
# sporozoite or other source-selection-defined cell.
harmonized_gametocyte <- pv.combined.all$harmonized_life_cycle_stage %in%
  c("Male gametocyte", "Female gametocyte")
allowed_gametocyte_route <- (
  pv.combined.all$source_stage_provenance == "marker_cluster_inferred" &
    !is.na(pv.combined.all$pred_gametocyte_sex) &
    pv.combined.all$harmonized_life_cycle_stage == pv.combined.all$pred_gametocyte_sex
) | (
  pv.combined.all$source_stage_provenance == "source_state_defined" &
    !is.na(pv.combined.all$source_sex_annotation) &
    pv.combined.all$harmonized_life_cycle_stage == pv.combined.all$source_sex_annotation
)
if (any(harmonized_gametocyte & !allowed_gametocyte_route) ||
  any(harmonized_gametocyte & pv.combined.all$tissue_or_sample_type != "Host blood")) {
  pipeline_fail(sprintf(
    "%d cell(s) carry a gametocyte harmonized_life_cycle_stage outside the source-label/cluster-call precedence routes or outside host blood",
    sum(harmonized_gametocyte & (!allowed_gametocyte_route | pv.combined.all$tissue_or_sample_type != "Host blood"))
  ))
}

near_zero <- pv.combined.all$total_umi_count < MIN_UMI_FOR_IDC_INFERENCE
if (any(!is.na(pv.combined.all$idc_reference_similarity_label[near_zero])) ||
  any(!is.na(pv.combined.all$pred_gametocyte_sex[near_zero]))) {
  pipeline_fail("AT11 regression: a near-zero-expression cell received idc_reference_similarity_label or pred_gametocyte_sex")
}

expected_broad <- unname(map_broad_stage(pv.combined.all$harmonized_life_cycle_stage))
actual_broad <- unname(pv.combined.all$parasite_broad_stage)
# Value-wise comparison (not identical()): subset() can change attributes
# such as names/factor-vs-character on metadata columns without changing
# the underlying per-cell values - a false-positive-prone identical()
# check here would fail-loud on a purely cosmetic difference. NA must
# align position-for-position on both sides; every non-NA pair must be
# equal.
mismatch <- !(
  (is.na(expected_broad) & is.na(actual_broad)) |
    (!is.na(expected_broad) & !is.na(actual_broad) & expected_broad == actual_broad)
)
if (any(mismatch)) {
  pipeline_fail(sprintf(
    "AT21 regression: parasite_broad_stage no longer matches the deterministic mapping of harmonized_life_cycle_stage (%d cell(s) mismatched)",
    sum(mismatch)
  ))
}

# rRNA genes are removed before normalization in integration.R; none may
# reach the final object.
assert_no_rrna_features(rownames(pv.combined.all), gff_rrna_gene_ids(), label = "pv_all_studies.rds (post-singleR.R)")

# AT31: every source-defined population and each study's harmonized
# gametocyte total must keep at least its recorded cell count.
floor_check <- assert_stage_population_floors(pv.combined.all@meta.data)
cat(sprintf("[singleR.R] AT31 stage-population floors: %d/%d PASS\n", sum(floor_check$status == "PASS"), nrow(floor_check)))

# ============================================================================
# Cleanup and save
# ============================================================================
features <- gsub("-", "_", rownames(pv.combined.all[["RNA"]]@features))
rownames(pv.combined.all[["RNA"]]@features) <- features

pv.combined.all$orig.ident <- NULL
pv.combined.all$femaleGams_module_score1 <- NULL
pv.combined.all$maleGams_module_score1 <- NULL
pv.combined.all@project.name <- "PlaViSca"

saveRDS(pv.combined.all, file = "./pv_all_studies.rds")
record_build_manifest(
  artifact_path = "pv_all_studies.rds",
  script_path = "scripts/singleR.R",
  cell_count = ncol(pv.combined.all),
  notes = sprintf(
    "Fixes D014,D015,D020,D021,D034,D037,D044,D045,D047,D048,D049,D058,D059,D060,D061; %d cells removed per DEC10 post-integration policy",
    n_removed
  )
)
