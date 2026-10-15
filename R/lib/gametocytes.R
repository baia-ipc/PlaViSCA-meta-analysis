# Gametocyte call per blood-stage cluster, from the detection of marker genes.
#
# For each cluster and marker set, the enrichment is the share of the
# cluster's cells detecting the marker genes (counts > 0) minus the share
# detecting control genes of similar expression, in percentage points. The
# control genes of a marker gene are n_controls random genes (not markers)
# of the same expression bin (mean normalised expression over all cells,
# n_bins bins of equal size). Comparing with genes of similar expression
# makes the enrichment independent of sequencing depth: deeper clusters
# detect more of every gene, markers or not.

# counts, data: genes x cells (counts and normalised); cluster: one value per
# cell; marker_sets: named list of gene vectors. Returns clusters x sets.
marker_enrichment <- function(counts, data, cluster, marker_sets, n_bins, n_controls, seed) {
  cluster <- factor(cluster)
  membership <- Matrix::sparse.model.matrix(~ 0 + cluster)
  colnames(membership) <- levels(cluster)
  detection <- as.matrix((counts > 0) %*% membership) /
    matrix(Matrix::colSums(membership), nrow(counts), ncol(membership), byrow = TRUE)
  mean_expression <- Matrix::rowMeans(data)
  bin <- setNames(cut(rank(mean_expression, ties.method = "first"), n_bins, labels = FALSE), rownames(data))
  candidates <- setdiff(rownames(data), unlist(marker_sets))
  set.seed(seed)
  sapply(marker_sets, function(genes) {
    controls <- unlist(lapply(genes, function(g) sample(candidates[bin[candidates] == bin[g]], n_controls)))
    100 * (colMeans(detection[genes, , drop = FALSE]) - colMeans(detection[controls, , drop = FALSE]))
  })
}

# enrichment: data frame with blood_cluster, female, male, male_mature
# (percentage points) and relative_umi (cluster median UMIs / median over the
# clusters). p: the annotation parameters. Adds call and state.
call_gametocyte_clusters <- function(enrichment, p) {
  enrichment |>
    dplyr::mutate(
      call = dplyr::case_when(
        female >= p$gametocyte_min_enrichment & female > male ~ "Female gametocyte",
        male >= p$gametocyte_min_enrichment & male > female ~ "Male gametocyte",
        TRUE ~ NA_character_),
      state = dplyr::case_when(
        call == "Male gametocyte" & male_mature >= p$late_male_min_enrichment ~ "Late male",
        call == "Male gametocyte" ~ "Early male",
        call == "Female gametocyte" & relative_umi < p$developing_max_rel_umi ~ "Developing gametocyte",
        call == "Female gametocyte" & relative_umi >= p$late_female_min_rel_umi ~ "Late female",
        call == "Female gametocyte" ~ "Early female",
        TRUE ~ NA_character_),
      # developing gametocytes: on the female branch by the markers, but too
      # early (few UMIs) for the sex to be called
      call = ifelse(state %in% "Developing gametocyte", "Developing gametocyte", call))
}
