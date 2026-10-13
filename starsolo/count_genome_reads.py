#!/usr/bin/env python3
"""
Count the unique reads of each cell barcode over the whole genome, from the
coordinate-sorted alignments of one STARsolo run.

Reads the alignments as SAM on standard input:
    samtools view <bam> | count_genome_reads.py <out_tsv>

Usage:
    count_genome_reads.py [options] <out_tsv>

Arguments:
    <out_tsv>     table: barcode, unique_reads (sorted by unique_reads)

Options:
    -h --help     Show this help message and exit.
    -V --version  Show the script version and exit.
    -v --verbose  Enable verbose logging.
    -q --quiet    Suppress log output.

A read counts if it is a primary alignment, maps to one place only (NH:i:1)
and has a corrected cell barcode (CB) and UMI (UB). Reads with the same
barcode, UMI, strand and position are counted once (PCR duplicates). Reads are
counted wherever they map, inside or outside genes.

This is the cell definition of Hazzard et al. 2022 (PMID 36525464): a barcode
is a cell if it has enough unique reads over the genome.
"""
import sys

from docopt import docopt
from loguru import logger


def count_unique_reads(sam) -> dict:
    counts = {}
    seen = set()
    position = None
    for line in sam:
        fields = line.rstrip("\n").split("\t")
        flag = int(fields[1])
        if flag & 0x904:  # unmapped, secondary or supplementary
            continue
        tags = dict(t.split(":", 2)[0::2] for t in fields[11:])
        barcode, umi = tags.get("CB", "-"), tags.get("UB", "-")
        if tags.get("NH") != "1" or barcode == "-" or umi == "-":
            continue
        if (fields[2], fields[3]) != position:  # sorted input: duplicates share a position
            position = (fields[2], fields[3])
            seen.clear()
        key = (barcode, umi, flag & 16)
        if key in seen:
            continue
        seen.add(key)
        counts[barcode] = counts.get(barcode, 0) + 1
    return counts


def set_logger(verbose: bool, quiet: bool) -> None:
    logger.remove()
    if quiet:
        logger.add(lambda msg: None, level="DEBUG")
    elif verbose:
        logger.add(lambda msg: sys.stderr.write(msg), level="INFO")
    else:
        logger.add(lambda msg: sys.stderr.write(msg), level="WARNING")


def count_genome_reads(out_tsv: str) -> None:
    counts = count_unique_reads(sys.stdin)
    with open(out_tsv, "w") as out:
        out.write("barcode\tunique_reads\n")
        for barcode, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            out.write(f"{barcode}\t{n}\n")
    logger.info(f"{len(counts)} barcodes, {sum(counts.values())} unique reads")


if __name__ == "__main__":
    args = docopt(__doc__, version="count_genome_reads.py 1.0")
    set_logger(args["--verbose"], args["--quiet"])
    count_genome_reads(args["<out_tsv>"])
