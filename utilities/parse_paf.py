import argparse

"""
This script processes a PAF file to filter and print alignment information based on a minimum length criterion.
It takes input from a PAF file and prints the filtered results to a specified output file.

Usage:
    python parse_paf.py --paf_file <path_to_paf_file> [--min_length <minimum_alignment_length>] --output_file <path_to_output_file>

Example:
    python parse_paf.py --paf_file alignments.paf --min_length 100000 --output_file filtered_alignments.txt
"""

################################################################################
##                           Function Definitions                             ##
################################################################################

def parse_paf_line(line, min_length=100000):
    """
    Parses a single line from a PAF file, filtering based on minimum length.

    Args:
        line: String containing a single line from the PAF file.
        min_length: Minimum alignment length (default: 100000).

    Returns:
        A dictionary containing relevant information from the PAF line if it meets the
        minimum length criteria, otherwise None.
    """
    fields = line.strip().split("\t")
    if len(fields) < 12 or int(fields[9]) < min_length:
        return None

    # Extracting fields and optional tags for further processing
    # See PAF spec @ https://github.com/lh3/miniasm/blob/master/PAF.md
    parsed_line = {
        'query_name': fields[0],
        'query_start': fields[2],
        'query_end': fields[3],
        'target_name': fields[5],
        'ref_pos_start': fields[7],
        'ref_pos_end': fields[8],
        'ref_seq_len': fields[6],
        'tp': None,  # Type of alignment
        'de': None,  # Sequence divergence
    }

    # Parsing optional fields for 'tp' and 'de'
    for field in fields[12:]:
        tag, _, value = field.split(':')
        if tag == 'tp':
            parsed_line['tp'] = value
        elif tag == 'de':
            parsed_line['de'] = value

    # Filtering non-primary/secondary alignments
    if parsed_line['tp'] not in ['P', 'I'] or parsed_line['de'] is None:
        return None

    return parsed_line


def read_paf(filename, min_length=100000):
    """
    Reads a PAF file and returns a list of parsed alignments.

    Args:
        filename: Path to the PAF file.
        min_length: Minimum alignment length (default: 100000).

    Returns:
        A list of dictionaries, where each dictionary contains parsed information from
        a valid alignment line in the PAF file.
    """
    alignments = []
    with open(filename, 'r') as f:
        for line in f:
            parsed_line = parse_paf_line(line, min_length)
            if parsed_line:
                alignments.append(parsed_line)
    return alignments


def print_lines_to_file(alignments, output_file):
    """
    Prints all elements of alignments along with a header line to a specified file.

    Args:
        alignments: List of dictionaries containing parsed information from alignments.
        output_file: Path to the file where the output will be written.
    """

    with open(output_file, 'w') as f:
        if alignments:
            # Generate header from dictionary keys
            header = '\t'.join(list(alignments[0].keys()))  + '\n'
            f.write(header)

            # Iterate over each alignment and print its values along with the label
            for alignment in alignments:
                line = '\t'.join(str(value) for value in alignment.values()) + '\n'
                f.write(line)



################################################################################
##                                  MAIN                                      ##
################################################################################

def main(args):
    alignments = read_paf(args.paf_file, args.min_length)
    print_lines_to_file(alignments, args.output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process phased alignment data.")
    parser.add_argument("--paf_file", type=str, required=True, help="Path to the PAF file.")
    parser.add_argument("--min_length", type=int, default=100000, help="Minimum alignment length (default: 100000).")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output file where results will be written.")
    args = parser.parse_args()

    main(args)
