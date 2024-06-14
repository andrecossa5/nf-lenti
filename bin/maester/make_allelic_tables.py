#!/usr/bin/python

"""
Make allelic tables.
"""
 

##


# Libraries
import os
import sys
import argparse


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='make_allelic_tables',
    description=
    """
    Create cell-specific allelic tables.
    """
)

# Input
# my_parser.add_argument(
#     '--sample', 
#     type=str,
#     default=None,
#     help='Sample name. Default: None.'
# )




##


# Parse arguments
args = my_parser.parse_args()
input_bam = args.input_bam
cell = args.cell
fasta_MT = args.fasta_MT
min_base_qual = args.min_base_qual
min_alignment_quality = args.min_alignment_quality


##


def main():





    






# Run 
if __name__ == '__main__':
    main()
