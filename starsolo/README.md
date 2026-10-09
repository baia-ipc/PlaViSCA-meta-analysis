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
Sa2020 reads are first trimmed to 75 bases (seqtk). The runs to count are
those of `metadata/runs.tsv`.

## How this reproduces the counts used

The STARsolo outputs in `inputs/starsolo_counts/` were made with the same
STAR version and options. A first recomputation, with the same-length
contigs paired in numeric order, differed from them mainly in the genes of
the wrongly paired contigs (and in a few other genes by 1 to 7 UMIs); one
run (SRR20710501) was identical. This led to `contig_pairs.tsv`. The
comparison with the corrected reference is in progress.

The barcode whitelists `3M-february-2018.txt.gz` and
`737K-august-2016.txt.gz` come from
<https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/>. Each is
identical to the barcode list of the raw matrices of its chemistry.

The order of the genes in the GTF (and therefore the order of the rows of
the matrices) follows the GFF. It differs from the order in
`inputs/starsolo_counts/`; the counts do not.
