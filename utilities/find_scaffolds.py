import argparse
import gzip
import re
from Bio import SeqIO

"""
Usage:
    python find_scaffolds.py --fasta_file <fasta_file_path> --bed_file <out_bed_path>

Example:
    python find_scaffolds.py --fasta_file genome.fasta --bed_file scaffolds.bed

Parse FASTA(.gz) file and outpu BED file which represent stretches of 'N' characters.
"""


################################################################################
##                           Function Definitions                             ##
################################################################################

def find_scaffolds(sequence, chrom_name):
    """
    Identifies the positions of 'N' characters in a genomic sequence.

    Args:
        sequence: The genomic sequence as a Biopython Seq object.
        chrom_name: The name of the chromosome or contig.

    Returns:
        A list of tuples, each representing the start and end positions of a scaffold.
    """
    scaffolds = []
    pattern = re.compile('N+', re.IGNORECASE)  

    sequence_str = str(sequence)

    for match in pattern.finditer(sequence_str):
        start, end = match.start(), match.end()
        scaffolds.append((chrom_name, start, end))

    return scaffolds

def read_fasta_file(filename):
    """
    Reads a FASTA file using Biopython and processes each sequence to find scaffolds.

    Args:
        filename: Path to the FASTA file.

    Returns:
        A list of all scaffolds found in all sequences within the file.
    """
    all_scaffolds = []

    open_func = gzip.open if filename.endswith('.gz') else open
    with open_func(filename, 'rt') as f:
        for record in SeqIO.parse(f, "fasta"):
            chrom_name = record.id
            sequence = record.seq
            scaffolds = find_scaffolds(sequence, chrom_name)
            all_scaffolds.extend(scaffolds)

    return all_scaffolds

def write_bed_file(scaffolds, output_file):
    """
    Writes the scaffold positions to a BED file.

    Args:
        scaffolds: A list of tuples with scaffold positions.
        output_file: Path to the output BED file.
    """
    with open(output_file, 'w') as f:
        for chrom, start, end in scaffolds:
            f.write(f"{chrom}\t{start}\t{end}\tN_stretch\n")

################################################################################
##                                  MAIN                                      ##
################################################################################

def main(args):

    scaffolds = read_fasta_file(args.fasta_file)
    write_bed_file(scaffolds, args.bed_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find scaffolds in FASTA sequences using Biopython and output their positions in a BED file.")
    parser.add_argument("--fasta_file", type=str, required=True, help="Path to the input FASTA file.")
    parser.add_argument("--bed_file", type=str, required=True, help="Path to the output BED file.")

    args = parser.parse_args()

    main(args)