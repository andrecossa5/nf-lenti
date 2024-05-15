#!/usr/bin/python

import dnaio
import argparse
import pandas as pd
from scipy.spatial.distance import hamming


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='get_triplets',
    description=
    """
    Get grepped file. 
    """
)

my_parser.add_argument(
    '--output_path', 
    type=str,
    default=None,
    help='Where to save the files'
)

my_parser.add_argument(
    '--read1_path', 
    type=str,
    default=None,
    help='Read1'
)

my_parser.add_argument(
    '--read2_path', 
    type=str,
    default=None,
    help='read2'
)

my_parser.add_argument(
    '--barcodes_path', 
    type=str,
    default=None,
    help='path of the barcodes'
)

my_parser.add_argument(
    '--gbc',
    action='store_true',
    help='filtering for lentiviral. Default: False.'
)

my_parser.add_argument(
    '--anchor',
    type=str,
    default=None,
    help='anchor for lentiviralbarcoding'
)

my_parser.add_argument(
    '--threshold',
    type=int,
    default=1,
    help='hamming distance from the anchor'
)

# Parse arguments
args = my_parser.parse_args()
input_read_1 = args.read1_path
input_read_2 = args.read2_path
barcodes_path = args.barcodes_path
output_read = args.output_path



##



solo_CBCs = pd.read_csv(
    barcodes_path, 
    header=None, index_col=0
)
solo_CBCs = list(solo_CBCs.index.unique())
        
with dnaio.open(input_read_1, input_read_2) as reader, open(output_read,'w') as writer:
    writer.write(f"Name\tcbc\tumi\tfeature\n")
    for record_1, record_2 in reader:
        cbc = record_1.sequence[:16]
        umi = record_1.sequence[16:]
        feature = record_2.sequence
        if args.gbc==True:
                    a_ = feature[:33]
                    gbc = feature[33:33+18]
                    h = hamming(list(a_), list(args.anchor)) * 33
                    if h <= int(args.threshold) and cbc in solo_CBCs:
                        feature = gbc
                        writer.write(f"{record_1.name}\t{cbc}\t{umi}\t{feature}\n")
        else:
            if cbc in solo_CBCs:
                 writer.write(f"{record_1.name}\t{cbc}\t{umi}\t{feature}\n")
