# STARsolo: from reads to count matrices

These scripts make the input of stage 10 (`<workspace>/inputs/starsolo_counts/`)
from the sequencing reads of the four studies counted with STARsolo (Sa2020,
Ruberto2022_2, Hazzard2022, Hazzard2024). They have their own software
environment (`pixi.toml` here: STAR 2.7.11b, seqtk, Python).

## Reference

`make_index.sh` builds, in `<workspace>/inputs/reference/`:

| File | What it is |
|---|---|
| `PlasmoDB-68_PvivaxP01_Genome.fasta` | the P. vivax P01 genome |
| `PlasmoDB-68_PvivaxP01_utr2exon.gtf` | the gene models of the PlasmoDB-68 GFF |
| `sequence_names.tsv` | PlasmoDB and ENA name of each sequence |
| `star_index_overhang74`, `star_index_overhang90` | STAR indices |

Starting files:
- `PlasmoDB-68_PvivaxP01.gff`, the PlasmoDB-68 annotation;
- `GCA_900093555.2.fasta`, the ENA assembly of P01:
  `https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_900093555.2`.

**Genome.** The genome sequences come from ENA, because PlasmoDB no longer
serves release-68 files.
- They are renamed to the PlasmoDB names used in the GFF
  (`make_reference.py`): chromosomes and organelle genomes by name,
  unplaced contigs by length.
- In 8 cases two or three unplaced contigs share the same length; each of
  these 8 contigs carries one gene. They are paired as listed in
  `contig_pairs.tsv`. For 7 of them the pairing is fixed by the splice
  motifs of the junctions found in the original STARsolo runs: only one
  candidate sequence has the right motif at every junction. Contig
  `Transfer.PvP01_00_208.final` has no junctions; of its two candidates,
  only one gives counts identical to the original (run ERR5087440).
- The 64 ENA contigs without genes are kept, under their ENA names, so that
  reads mapping there are not counted as unique.

**GTF.** One transcript per gene. Its exon lines are the GFF `exon`,
`five_prime_UTR` and `three_prime_UTR` records, so the UTRs are repeated as
exons.

**STAR indices.** `--genomeSAindexNbases 11`, with `--sjdbOverhang 74` for
the 10x v2 study (Sa2020) and 90 for the 10x v3 studies.

## Counting

```bash
starsolo/run_starsolo.sh <study_pmid> <run_id> <fastq_dir> <out_dir> [threads]
```

This counts one run (`<fastq_dir>/<pmid>/<run>_1.fastq.gz` and `_2.fastq.gz`)
into `<out_dir>/<pmid>/<run>_solo_out/`. The STAR options are in the script.
The folder gets a `SHA256SUMS` with the checksums of the files stage 10
reads (logs, `Solo.out/GeneFull/`, `genome_reads_per_barcode.tsv`), which
`./run.sh check-inputs` verifies.
Only read 2 (cDNA) is aligned; read 1 gives the cell barcode and the UMI.
Sa2020 reads are first trimmed to 75 bases (seqtk). The runs to count are
those of `metadata/runs.tsv`.

## Unique reads per barcode over the genome

```bash
starsolo/run_starsolo.sh --genome-reads 36525464 SRR20710500 <fastq_dir> <out_dir>
```

With `--genome-reads`, the run's folder also gets
`genome_reads_per_barcode.tsv` (barcode, unique reads). Stage 10 uses it
for Hazzard2022, whose authors define a cell as a barcode with at least
5,000 (blood) or 250 (sporozoites) unique reads over the whole genome (the
paper says "more than"; their counting script uses `>= 5000`).

`count_genome_reads.py` counts, from a temporary coordinate-sorted BAM, the
reads that are primary alignments, map to one place only and have a
corrected barcode and UMI; reads with the same barcode, UMI, strand and
position are counted once, wherever they map. The count matrices are the
same with or without the option. The four Hazzard2022 runs of
`inputs/starsolo_counts/` were made with it.

## Why only read 2 is aligned

In 10x 3' libraries, read 2 is the cDNA, on the same strand as the RNA;
read 1 holds the barcode and the UMI, and in some runs (Sa2020, the
Hazzard2022 blood runs: 75 bases) also poly-T and adapter. STARsolo is
therefore given read 2 as the read to align and read 1 only as the barcode
read (`--readFilesIn R2 R1 --soloBarcodeReadLength 0`), and counts reads on
the gene strand (STARsolo default).

Aligning read 1 as a mate of read 2 instead (`--soloBarcodeMate 1` with
`--clip5pNbases`) gives wrong counts: the pair takes the strand of read 1, so
reads are counted on the opposite strand of the genes, and where read 1 is
longer than barcode + UMI most pairs are discarded as too short.

Check: in run ERR5087440 (Ruberto2022_2), for the 2,570 cells of the
authors' analysis, the gene counts correlate with the authors' counts with
a median per-cell Spearman correlation of 0.88.

The reads are the FASTQ files of the runs in `metadata/runs.tsv`, downloaded
from the SRA/ENA into `<fastq_dir>/<PMID>/`, with an `md5sum.txt` of the
files.

The barcode whitelists `3M-february-2018.txt.gz` and
`737K-august-2016.txt.gz` come from
<https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/>. Each is
identical to the barcode list of the raw matrices of its chemistry.

The order of the genes in the GTF, and therefore of the rows of the
matrices, follows the GFF.
