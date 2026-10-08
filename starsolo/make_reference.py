#!/usr/bin/env python3
"""
Build the STARsolo reference files from the PlasmoDB-68 PvivaxP01 GFF and the
ENA assembly of P. vivax P01 (GCA_900093555.2).

Usage:
    make_reference.py [options] <gff> <ena_fasta> <out_fasta> <out_gtf> <out_names>

Arguments:
    <gff>        PlasmoDB-68_PvivaxP01.gff
    <ena_fasta>  ENA assembly GCA_900093555.2 (FASTA)
    <out_fasta>  genome FASTA with the PlasmoDB sequence names
    <out_gtf>    GTF for STAR
    <out_names>  table: PlasmoDB sequence name, ENA sequence name

Options:
    -h --help     Show this help message and exit.
    -V --version  Show the script version and exit.
    -v --verbose  Enable verbose logging.
    -q --quiet    Suppress log output.

Genome: the ENA sequences are renamed to the PlasmoDB names (PvP01_01_v2,
..., Transfer.PvP01_00_1.final) and written in the order of the GFF
##sequence-region lines, followed by the ENA contigs without genes (which are
not in the GFF; they keep their ENA name). Chromosomes, apicoplast and mitochondrion are
matched by name, the unplaced contigs by length; when several ENA contigs
have the same length they are paired in numeric order. Every sequence must
have the length given in the GFF.

GTF: one transcript per gene (mRNA, ncRNA types, pseudogenic transcripts).
Its exon lines are the GFF exon, five_prime_UTR and three_prime_UTR records
(the UTRs are repeated as exons, as in the GTF used to make the counts).
"""
import collections
import re
import sys
from urllib.parse import unquote

from docopt import docopt
from loguru import logger


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(chunks)
            m = re.search(r"(?:chromosome|contig|scaffold): ([^,\s]+)", line)
            name, chunks = m.group(1).upper(), []
        else:
            chunks.append(line.strip())
    seqs[name] = "".join(chunks)
    return seqs


def sequence_regions(gff):
    return [(f[1], int(f[3])) for f in (l.split() for l in open(gff) if l.startswith("##sequence-region"))]


def ena_name(plasmodb_name):
    m = re.match(r"PvP01_(\d\d)_v\d", plasmodb_name)
    if m:
        return str(int(m.group(1)))
    for organelle in ("API", "MIT"):
        if plasmodb_name.startswith("PvP01_" + organelle):
            return organelle
    return None


def match_sequences(regions, ena):
    matched, used = {}, set()
    for name, length in regions:
        e = ena_name(name)
        if e:
            matched[name] = e
            used.add(e)
    contig_number = lambda n: int(re.search(r"(\d+)$", n.replace(".final", "")).group(1))
    by_length = collections.defaultdict(list)
    for e, s in ena.items():
        if e not in used:
            by_length[len(s)].append(e)
    for length in by_length:
        by_length[length].sort(key=contig_number)
    for name, length in sorted((r for r in regions if r[0] not in matched), key=lambda r: contig_number(r[0])):
        candidates = by_length.get(length)
        if not candidates:
            sys.exit(f"no ENA sequence of length {length} for {name}")
        matched[name] = candidates.pop(0)
    return matched


def write_fasta(regions, ena, matched, path):
    with open(path, "w") as out:
        records = [(name, ena[matched[name]], length) for name, length in regions]
        # contigs without genes (absent from the GFF) are kept, under their
        # ENA name: reads mapping there must not be counted as unique
        unused = sorted(set(ena) - set(matched.values()))
        records += [(name, ena[name], len(ena[name])) for name in unused]
        for name, seq, length in records:
            if len(seq) != length:
                sys.exit(f"{name}: length {len(seq)} in ENA, {length} in the GFF")
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def write_gtf(gff, path):
    transcript_types = {"mRNA", "ncRNA", "rRNA", "tRNA", "snoRNA", "snRNA", "pseudogenic_transcript"}
    exon_types = {"exon", "five_prime_UTR", "three_prime_UTR"}
    transcript_gene, rows = {}, []
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        attrs = dict(a.split("=", 1) for a in f[8].split(";") if "=" in a)
        if f[2] in transcript_types:
            transcript_gene[attrs["ID"]] = attrs["Parent"]
        elif f[2] in exon_types:
            rows.append((f, unquote(attrs["Parent"])))
    with open(path, "w") as out:
        for f, transcript in rows:
            gene = transcript_gene[transcript]
            out.write("\t".join(f[:2] + ["exon"] + f[3:8] +
                [f'gene_id "{gene}"; transcript_id "{transcript}";']) + "\n")


def set_logger(verbose: bool, quiet: bool) -> None:
    logger.remove()
    if quiet:
        logger.add(lambda msg: None, level="DEBUG")
    elif verbose:
        logger.add(lambda msg: sys.stderr.write(msg), level="INFO")
    else:
        logger.add(lambda msg: sys.stderr.write(msg), level="WARNING")


def make_reference(gff: str, ena_fasta: str, out_fasta: str, out_gtf: str, out_names: str) -> None:
    regions = sequence_regions(gff)
    ena = read_fasta(ena_fasta)
    matched = match_sequences(regions, ena)
    logger.info(f"{len(regions)} GFF sequences matched to ENA; {len(ena) - len(regions)} ENA contigs without genes added")
    write_fasta(regions, ena, matched, out_fasta)
    write_gtf(gff, out_gtf)
    with open(out_names, "w") as out:
        out.write("plasmodb_name\tena_name\n")
        for name, _ in regions:
            out.write(f"{name}\t{matched[name]}\n")


if __name__ == "__main__":
    args = docopt(__doc__, version="make_reference.py 1.0")
    set_logger(args["--verbose"], args["--quiet"])
    make_reference(args["<gff>"], args["<ena_fasta>"], args["<out_fasta>"], args["<out_gtf>"], args["<out_names>"])
