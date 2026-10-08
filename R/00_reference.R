# Stage 00: gene reference.
#
# Reads the PlasmoDB-68 PvivaxP01 GFF and the STARsolo gene lists, and writes:
#   genes.tsv            one row per gene (ID, name, description, type, rRNA flag)
#   starsolo_genes.tsv   the genes counted by STARsolo (identical in every run)
#   gene_annotation.rds  gene ID, name and description, as used by the app
#
# The GFF is parsed directly (tab-separated text) rather than with
# rtracklayer, which tries to contact NCBI for the species and fails offline.

STAGE <- "00_reference"
source("R/lib/setup.R")

out <- stage_dir(STAGE)

# --- genes from the GFF -----------------------------------------------------

gff_path <- input_path(PARAMS$reference$gff)
lines <- readLines(gff_path)
fields <- strsplit(lines[!startsWith(lines, "#")], "\t", fixed = TRUE)
fields <- fields[lengths(fields) == 9]
type <- vapply(fields, `[[`, "", 3)
attributes <- vapply(fields, `[[`, "", 9)

gene_types <- c("protein_coding_gene", "ncRNA_gene", "pseudogene")
is_gene <- type %in% gene_types

attribute <- function(attrs, key) {
  value <- ifelse(grepl(paste0("(^|;)", key, "="), attrs),
    sub(paste0("^.*?(^|;)", key, "=([^;]*).*$"), "\\2", attrs, perl = TRUE), NA_character_)
  # GFF3 percent-encodes reserved characters (%2C = ",", %3B = ";")
  vapply(value, function(v) if (is.na(v)) NA_character_ else utils::URLdecode(v), "", USE.NAMES = FALSE)
}

genes <- data.frame(
  gene_id = attribute(attributes[is_gene], "ID"),
  name = attribute(attributes[is_gene], "Name"),
  description = attribute(attributes[is_gene], "description"),
  type = type[is_gene],
  biotype = attribute(attributes[is_gene], "ebi_biotype"),
  stringsAsFactors = FALSE
)
check_unique(genes$gene_id, "gene IDs in the GFF")
# Seurat replaces "_" with "-" in feature names
genes$seurat_id <- gsub("_", "-", genes$gene_id)
# rRNA genes: ncRNA genes with biotype rRNA (protein-coding genes whose
# description mentions rRNA, e.g. rRNA methyltransferases, are not rRNA)
genes$is_rrna <- genes$type == "ncRNA_gene" & genes$biotype %in% "rRNA"
message(sprintf("%d genes in the GFF (%s), %d rRNA genes",
  nrow(genes), paste(names(table(genes$type)), table(genes$type), sep = ": ", collapse = ", "), sum(genes$is_rrna)))
check(sum(genes$is_rrna) > 0, "no rRNA genes found in ", gff_path)

write_table(genes, file.path(out, "genes.tsv"))

# Gene annotation for the app: genes with a standard PVP01_ + 7 digit ID
annotation <- genes |>
  filter(grepl("^PVP01_[0-9]{7}$", gene_id)) |>
  transmute(id = gene_id, name, description)
save_object(annotation, file.path(out, "gene_annotation.rds"))

# --- genes counted by STARsolo ----------------------------------------------

feature_files <- list.files(file.path(INPUTS, PARAMS$import$starsolo_counts), pattern = "^features.tsv$",
  recursive = TRUE, full.names = TRUE)
feature_files <- grep("/GeneFull/(raw|filtered)/features.tsv$", feature_files, value = TRUE)
check(length(feature_files) > 0, "no STARsolo features.tsv found under ", file.path(INPUTS, PARAMS$import$starsolo_counts))
feature_lists <- lapply(feature_files, function(f) utils::read.delim(f, header = FALSE)$V1)
check(all(vapply(feature_lists, identical, TRUE, feature_lists[[1]])),
  "the STARsolo runs were not all counted against the same gene list")
starsolo_genes <- data.frame(gene_id = feature_lists[[1]])
check(all(starsolo_genes$gene_id %in% genes$gene_id), "STARsolo counts genes that are not in the GFF")
message(sprintf("%d STARsolo gene lists, all identical: %d genes", length(feature_files), nrow(starsolo_genes)))
write_table(starsolo_genes, file.path(out, "starsolo_genes.tsv"))
