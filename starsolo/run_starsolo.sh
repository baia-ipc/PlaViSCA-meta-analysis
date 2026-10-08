#!/usr/bin/env bash
# Count the reads of one run with STARsolo (the input of stage 10).
#
# Usage: starsolo/run_starsolo.sh <study_pmid> <run_id> <fastq_dir> <out_dir> [threads]
#
#   <fastq_dir>/<study_pmid>/<run_id>_1.fastq.gz, _2.fastq.gz   reads (R1: barcode+UMI, R2: cDNA)
#   <out_dir>/<study_pmid>/<run_id>_solo_out/                   STARsolo output
#
# Reference files (made by starsolo/make_reference.py and STAR genomeGenerate,
# see starsolo/README.md) are taken from <workspace>/inputs/reference/.
#
# Settings by 10x chemistry (Sa2020, PMID 32365102, is 10x v2; the others v3):
#   v2: whitelist 737K-august-2016, UMI 10 nt, clip 26 nt of R1, reads
#       trimmed to 75 nt (seqtk trimfq), index with sjdbOverhang 74
#   v3: whitelist 3M-february-2018, UMI 12 nt, clip 28 nt of R1, index with
#       sjdbOverhang 90
set -euo pipefail
cd "$(dirname "$0")"

pmid=$1; run=$2; fastq_dir=$3; out_dir=$4; threads=${5:-15}
[ -f ../config/local.env ] && . ../config/local.env
ref="$PLAVISCA_WORKSPACE/inputs/reference"

if [ "$pmid" = 32365102 ]; then
  whitelist=737K-august-2016.txt; umi=10; clip="26 0"; index=star_index_overhang74; keep=75
else
  whitelist=3M-february-2018.txt; umi=12; clip="28 0"; index=star_index_overhang90; keep=
fi

mkdir -p "$out_dir/$pmid"
work=$(mktemp -d "$out_dir/$pmid/.$run.XXXX")
trap 'rm -rf "$work"' EXIT

r1="$fastq_dir/$pmid/${run}_1.fastq.gz"
r2="$fastq_dir/$pmid/${run}_2.fastq.gz"
if [ -n "$keep" ]; then
  # keep the first $keep bases of each read (3' trimming)
  for i in 1 2; do
    fq="$fastq_dir/$pmid/${run}_$i.fastq.gz"
    len=$(zcat "$fq" | awk 'NR == 2 { print length($0); exit }')
    if [ "$len" -gt "$keep" ]; then
      pixi run seqtk trimfq -e $((len - keep)) "$fq" | gzip -1 > "$work/${run}_$i.fastq.gz"
    else
      ln -s "$fq" "$work/${run}_$i.fastq.gz"
    fi
  done
  r1="$work/${run}_1.fastq.gz"; r2="$work/${run}_2.fastq.gz"
fi

zcat -f "$ref/$whitelist.gz" > "$work/$whitelist"
pixi run STAR \
  --runThreadN "$threads" \
  --genomeDir "$ref/$index" \
  --readFilesIn "$r1" "$r2" \
  --soloBarcodeMate 1 \
  --sjdbGTFfile "$ref/PlasmoDB-68_PvivaxP01_utr2exon.gtf" \
  --readFilesCommand zcat \
  --soloCBwhitelist "$work/$whitelist" \
  --soloType CB_UMI_Simple \
  --clip5pNbases $clip \
  --soloCBlen 16 \
  --soloUMIlen $umi \
  --soloUMIfiltering MultiGeneUMI \
  --alignIntronMin 1 \
  --alignIntronMax 2756 \
  --outSAMtype None \
  --soloFeatures Gene GeneFull \
  --outStd Log \
  --outTmpDir "$work/star_tmp" \
  --outFileNamePrefix "$work/${run}_solo_out/"

rm -rf "$out_dir/$pmid/${run}_solo_out"
mv "$work/${run}_solo_out" "$out_dir/$pmid/${run}_solo_out"
echo "done $pmid $run"
