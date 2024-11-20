"""
Usage:
    python asm_plot.py --input_file <input_file_path> --bed_file <bed_file_path> --output_file <output_file_path>

Example:
    python asm_plot.py --input_file data.txt --bed_file regions.bed --output_file plot.png

This script reads a tabulated input file and a BED file with genomic regions, parses the data, and generates a graphical representation using matplotlib. The plot visualizes genomic alignments and specific regions from the BED file on corresponding chromosome tracks.
"""

import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import csv
import re

ref_height  = 0.10
asm_height  = 0.05
asm_spacing = 0.15

################################################################################
##                           Function Definitions                             ##
################################################################################

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]

def parse_tab_line(line):
    """
    Parses a single line from the input tabulated file.

    Args:
        line: String containing a single line from the file.

    Returns:
        A dictionary with parsed values including target_name, ref_pos_start, ref_pos_end, and query_name.
    """
    fields = line.strip().split("\t")
    target_name, ref_pos_start, ref_pos_end, query_name = fields[3], int(fields[4]), int(fields[5]), fields[0]
    return {'chrom': target_name, 'start': ref_pos_start, 'end': ref_pos_end, 'query': query_name}

def read_tab_file(filename):
    """
    Reads the input file and parses each line.

    Args:
        filename: Path to the input file.

    Returns:
        A list of dictionaries, each containing the parsed line information.
    """
    data = []
    with open(filename, 'r') as f:
        next(f)  # Skip header line
        for line in f:
            if line.strip() and not line.startswith('#'):
                parsed_line = parse_tab_line(line)
                data.append(parsed_line)
    return data

def read_tsv_to_dict(filename, key_column=0, value_column=1):
    """
    Reads the input tsv and turns the first two columns into a dictionary

    Args:
        filename: Path to the input file.

    Returns:
        a dictionary where the first column is now the key (and the second column is the value)
    """

    with open(filename, mode='r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        result_dict = {rows[key_column]: rows[value_column] for rows in reader if len(rows) >= value_column}
    return result_dict

def parse_bed_line(line):
    """
    Parses a single line from the input BED file.

    Args:
        line: String containing a single line from the BED file.

    Returns:
        A dictionary with parsed values including chrom, start, end, type, and color.
    """
    fields = line.strip().split("\t")
    chrom, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
    # Calculate the size of the region and filter out if less than 100,000 bp
    size = end - start
    if size < 100000:
        return None  # Skip regions smaller than 100,000 bp

    type = re.split('_|\(', name)[0]  # Extract the type from the fourth column
    color_values = fields[8].split(',')  # Expecting RGB color in the ninth column
    if len(color_values) == 3:
        color = tuple(int(c) / 255.0 for c in color_values)  # Normalize to [0, 1]
    else:
        color = 'grey'  # Default color if parsing fails
    
    return {'chrom': chrom, 'start': start, 'end': end, 'type': type, 'color': color}



def read_bed_file(filename):
    """
    Reads the BED file and parses each line.

    Args:
        filename: Path to the BED file.

    Returns:
        A list of dictionaries, each containing the parsed line information.
    """
    regions = []
    with open(filename, 'r') as f:
        next(f)  # Skip header line
        for line in f:
            if line.strip() and not line.startswith('#'):
                parsed_line = parse_bed_line(line)
                if parsed_line:
                    regions.append(parsed_line)
    return regions

# Predefined dictionary for chromosome names and their sizes in base pairs
REFERENCE_CHROMOSOMES = {
    "chr1": 248387328, "chr2": 242696752, "chr3": 201105948, "chr4": 193574945,
    "chr5": 182045439, "chr6": 172126628, "chr7": 160567428, "chr8": 146259331,
    "chr9": 150617247, "chr10": 134758134, "chr11": 135127769, "chr12": 133324548,
    "chr13": 113566686, "chr14": 101161492, "chr15": 99753195, "chr16": 96330374,
    "chr17": 84276897, "chr18": 80542538, "chr19": 61707364, "chr20": 66210255,
    "chr21": 45090682, "chr22": 51324926, "chrX": 154259566, "chrY": 62460029,
    "chrM": 16569
}

def plot_data(alignment_files, color_files, regions, output_file):
    """
    Plots the parsed data and regions using matplotlib, saving the plot to an output file.

    Args:
        alignment_files: List of paths to the input tabulated files.
        regions: List of dictionaries containing parsed line information from the BED file.
        output_file: Path to the output file.
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    # Sort chromosome names from the predefined reference set
    chroms = sorted(REFERENCE_CHROMOSOMES.keys(), key=natural_sort_key, reverse=True)
    chrom_y_positions = {chrom: i for i, chrom in enumerate(chroms, start=1)}

    # Plot chromosomes as black lines based on their known sizes
    for chrom, y_pos in chrom_y_positions.items():
        if chrom in REFERENCE_CHROMOSOMES:
            chrom_length = REFERENCE_CHROMOSOMES[chrom]
            ax.add_patch(patches.Rectangle((0, y_pos - ref_height/2), chrom_length, ref_height, color='black', alpha=0.6))

    # Plot regions from the BED file on top of the chromosomes
    for r in regions:
        if r['chrom'] in REFERENCE_CHROMOSOMES:
            x_position = r['start']
            width = r['end'] - r['start']
            y_position = chrom_y_positions[r['chrom']]
            rect = patches.Rectangle((x_position, y_position - ref_height/2), width, ref_height, color=r['color'])
            ax.add_patch(rect)

    # Initialize a vertical offset for the first set of alignments
    vertical_offset = asm_spacing

    # Plot alignments from each file at a unique height
    for i in range(len(alignment_files)):
        data         = read_tab_file(alignment_files[i])
        color_lookup = read_tsv_to_dict(color_files[i])

        for d in data:
            if d['chrom'] in REFERENCE_CHROMOSOMES:
                query_name = d['query']
                color_val  = color_lookup[query_name]

                x_position = d['start']
                width = d['end'] - d['start']
                y_position = chrom_y_positions[d['chrom']] + vertical_offset
                rect = patches.Rectangle((x_position, y_position), width, asm_height, color=color_val)
                ax.add_patch(rect)

        # Increment the vertical offset for the next file's alignments
        vertical_offset += asm_spacing

    # Set plot parameters
    max_bp_length = max(REFERENCE_CHROMOSOMES.values())
    ax.set_xlim(0, max_bp_length)
    ax.set_ylim(0, len(chroms) + 1)
    ax.set_yticks(range(1, len(chroms) + 1))
    ax.set_yticklabels(chroms)
    ax.set_xlabel('Position')
    ax.set_ylabel('Chromosome')
    plt.tight_layout()

    # Save the plot to file
    plt.savefig(output_file, dpi=600)


################################################################################
##                                  MAIN                                      ##
################################################################################

def main(args):

    regions = read_bed_file(args.bed_file) if args.bed_file else []
    plot_data(args.input_files, args.color_files, regions, args.output_file)    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize multiple genomic alignments with BED file regions.")
    parser.add_argument("--input_files", nargs='+', type=str, required=True, help="Paths to the input tabulated files (up to 10).")
    parser.add_argument("--color_files", nargs="+", type=str, required=True, help="Paths to tsv file with colors for query/contig names")
    parser.add_argument("--bed_file", type=str, help="Path to the BED file with genomic regions (optional).")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output plot file.")
    args = parser.parse_args()

    main(args)
