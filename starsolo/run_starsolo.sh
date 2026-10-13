#!/usr/bin/env bash
# Count the reads of one run with STARsolo (the input of stage 10).
#
# Usage: starsolo/run_starsolo.sh [--genome-reads] <study_pmid> <run_id> <fastq_dir> <out_dir> [threads]
#
#   <fastq_dir>/<study_pmid>/<run_id>_1.fastq.gz, _2.fastq.gz   reads (R1: barcode+UMI, R2: cDNA)
#   <out_dir>/<study_pmid>/<run_id>_solo_out/                   STARsolo output
#   <out_dir>/<study_pmid>/<run_id>_solo_out/SHA256SUMS         checksums of the files stage 10 reads
#
# --genome-reads: also write <run_id>_solo_out/genome_reads_per_barcode.tsv,
#   the unique reads of each barcode over the whole genome (count_genome_reads.py),
#   from a temporary coordinate-sorted BAM. Stage 10 uses it for Hazzard2022.
#   The count matrices are the same with or without it.
#
# Reference files (made by starsolo/make_reference.py and STAR genomeGenerate,
# see starsolo/README.md) are taken from <workspace>/inputs/reference/.
#
# Settings by 10x chemistry (Sa2020, PMID 32365102, is 10x v2; the others v3):
#   v2: whitelist 737K-august-2016, UMI 10 nt, reads trimmed to 75 nt
#       (seqtk trimfq), index with sjdbOverhang 74
#   v3: whitelist 3M-february-2018, UMI 12 nt, index with sjdbOverhang 90
#
# Only R2 (cDNA) is aligned. R1 gives the barcode (first 16 nt) and the UMI
# (next 10 or 12 nt); the rest of R1, when present (Sa2020 and the
# Hazzard2022 blood runs have 75-nt R1), is poly-T and adapter and is ignored.
set -euo pipefail
cd "$(dirname "$0")"

genome_reads=
if [ "${1:-}" = "--genome-reads" ]; then genome_reads=1; shift; fi
pmid=$1; run=$2; fastq_dir=$3; out_dir=$4; threads=${5:-15}
[ -f ../config/local.env ] && . ../config/local.env
# STARSOLO_REFERENCE: another reference folder (default: the workspace one)
ref="${STARSOLO_REFERENCE:-$PLAVISCA_WORKSPACE/inputs/reference}"

if [ "$pmid" = 32365102 ]; then
  whitelist=737K-august-2016.txt; umi=10; index=star_index_overhang74; keep=75
else
  whitelist=3M-february-2018.txt; umi=12; index=star_index_overhang90; keep=
fi

mkdir -p "$out_dir/$pmid"
work=$(mktemp -d "$out_dir/$pmid/.$run.XXXX")
trap 'rm -rf "$work"' EXIT

# reads already trimmed to $keep bases are named <run>_trimmed_<1|2>.fastq.gz
reads=$run
[ -f "$fastq_dir/$pmid/${run}_1.fastq.gz" ] || reads=${run}_trimmed
r1="$fastq_dir/$pmid/${reads}_1.fastq.gz"
r2="$fastq_dir/$pmid/${reads}_2.fastq.gz"
if [ -n "$keep" ]; then
  # keep the first $keep bases of each read (3' trimming)
  for i in 1 2; do
    fq="$fastq_dir/$pmid/${reads}_$i.fastq.gz"
    len=$( (zcat "$fq" || true) | awk 'NR == 2 { print length($0); exit }')
    if [ "$len" -gt "$keep" ]; then
      pixi run seqtk trimfq -e $((len - keep)) "$fq" | gzip -1 > "$work/${run}_$i.fastq.gz"
    else
      ln -s "$fq" "$work/${run}_$i.fastq.gz"
    fi
  done
  r1="$work/${run}_1.fastq.gz"; r2="$work/${run}_2.fastq.gz"
fi

zcat -f "$ref/$whitelist.gz" > "$work/$whitelist"
alignments=(--outSAMtype None)
[ -n "$genome_reads" ] && alignments=(--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM CB UB \
  --limitBAMsortRAM 20000000000 --outBAMsortingThreadN "$threads")
pixi run STAR \
  --runThreadN "$threads" \
  --genomeDir "$ref/$index" \
  --readFilesIn "$r2" "$r1" \
  --soloBarcodeReadLength 0 \
  --sjdbGTFfile "$ref/PlasmoDB-68_PvivaxP01_utr2exon.gtf" \
  --readFilesCommand zcat \
  --soloCBwhitelist "$work/$whitelist" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 \
  --soloUMIlen $umi \
  --soloUMIfiltering MultiGeneUMI \
  --alignIntronMin 1 \
  --alignIntronMax 2756 \
  "${alignments[@]}" \
  --soloFeatures Gene GeneFull \
  --outStd Log \
  --outTmpDir "$work/star_tmp" \
  --outFileNamePrefix "$work/${run}_solo_out/"

if [ -n "$genome_reads" ]; then
  bam="$work/${run}_solo_out/Aligned.sortedByCoord.out.bam"
  pixi run samtools view "$bam" | pixi run python count_genome_reads.py "$work/${run}_solo_out/genome_reads_per_barcode.tsv"
  rm "$bam"
fi

# checksums of the files stage 10 reads, checked by ./run.sh check-inputs
(cd "$work/${run}_solo_out" &&
  checked=(Log.out Log.final.out Solo.out/GeneFull) &&
  { [ -z "$genome_reads" ] || checked+=(genome_reads_per_barcode.tsv); } &&
  find "${checked[@]}" -type f | LC_ALL=C sort | xargs sha256sum > SHA256SUMS)

rm -rf "$out_dir/$pmid/${run}_solo_out"
mv "$work/${run}_solo_out" "$out_dir/$pmid/${run}_solo_out"
echo "done $pmid $run"
