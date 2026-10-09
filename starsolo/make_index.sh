#!/usr/bin/env bash
# Build the STARsolo reference in <workspace>/inputs/reference/:
#   PlasmoDB-68_PvivaxP01_Genome.fasta   genome (ENA GCA_900093555.2, PlasmoDB names)
#   PlasmoDB-68_PvivaxP01_utr2exon.gtf   gene models (from the PlasmoDB-68 GFF)
#   sequence_names.tsv                   PlasmoDB and ENA sequence names
#   star_index_overhang74, _overhang90   STAR indices (10x v2 and v3 studies)
# Needs PlasmoDB-68_PvivaxP01.gff and GCA_900093555.2.fasta in that folder.
set -euo pipefail
cd "$(dirname "$0")"
[ -f ../config/local.env ] && . ../config/local.env
ref="${STARSOLO_REFERENCE:-$PLAVISCA_WORKSPACE/inputs/reference}"

pixi run python -W ignore make_reference.py -v "$ref/PlasmoDB-68_PvivaxP01.gff" "$ref/GCA_900093555.2.fasta" contig_pairs.tsv \
  "$ref/PlasmoDB-68_PvivaxP01_Genome.fasta" "$ref/PlasmoDB-68_PvivaxP01_utr2exon.gtf" "$ref/sequence_names.tsv"
for overhang in 74 90; do
  index="$ref/star_index_overhang$overhang"
  rm -rf "$index"; mkdir -p "$index"
  pixi run STAR --runThreadN 8 --runMode genomeGenerate --genomeDir "$index" \
    --genomeFastaFiles "$ref/PlasmoDB-68_PvivaxP01_Genome.fasta" --genomeSAindexNbases 11 \
    --sjdbGTFfile "$ref/PlasmoDB-68_PvivaxP01_utr2exon.gtf" --sjdbOverhang $overhang \
    --outFileNamePrefix "$index/"
done
