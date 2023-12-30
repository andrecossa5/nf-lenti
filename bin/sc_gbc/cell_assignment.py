#!/usr/bin/python

# Cell assignment script

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='cell_assignment',
    description=
    """
    Script for clone calling and cell assignment.
    """
)

# Input
my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)

# Output
my_parser.add_argument(
    '--path_bulk', 
    type=str,
    default=None,
    help='Path to input bulk reference. Default: None.'
)

# treshold
my_parser.add_argument(
    '--path_sc',
    type=str,
    default=None,
    help='Path to input sc GBC reads elements. Default: None.'
)

# treshold
my_parser.add_argument(
    '--sample_map',
    type=str,
    default=None,
    help='Path to sample_map. Default: None.'
)

# treshold
my_parser.add_argument(
    '--ncores',
    type=int,
    default=8,
    help='n cores for pairwise distances calculation. Default: 8.'
)

# Spikeins
my_parser.add_argument(
    '--bulk_correction_treshold',
    type=int,
    default=1,
    help='''
    Hamming distance treshold to consider a sc GBC a "degenerate" sequence with respect 
    to a bulk reference sequence. Default: 1.
    '''
)

# Spikeins
my_parser.add_argument(
    '--umi_treshold',
    type=int,
    default=5,
    help='Min number of UMIs to consider a CBC-GBC combination supported. Default: 5.'
)

# Spikeins
my_parser.add_argument(
    '--coverage_treshold',
    type=int,
    default=None,
    help='Min coverage (i.e., n reads) to consider a CBC-GBC-UMI combination supported. Default: None.'
)

# p_treshold
my_parser.add_argument(
    '--p_treshold',
    type=int,
    default=.001,
    help='Max p_poisson treshold to consider a CBC-GBC combination supported. Default: .001.'
)

# ratio_to_most_abundant_treshold
my_parser.add_argument(
    '--ratio_to_most_abundant_treshold',
    type=float,
    default=.3,
    help='Min coverage (nUMIs / nreads) to consider a CB-GBC combination supported. Default: .3.'
)


##


# Parse arguments
args = my_parser.parse_args()
sample = args.sample
path_bulk = args.path_bulk
path_sample_map = args.sample_map
path_sc = args.path_sc
ncores = args.ncores
bulk_correction_treshold = args.bulk_correction_treshold
umi_treshold = args.umi_treshold
coverage_treshold = args.coverage_treshold
p_treshold = args.p_treshold
ratio_to_most_abundant_treshold = args.ratio_to_most_abundant_treshold

# sample = 'AML_clones'
# path_bulk = '/Users/IEO5505/Desktop/example_mito/scratch_data/bulk_GBC_reference.csv'
# path_sample_map = '/Users/IEO5505/Desktop/example_mito/scratch_data/sample_map.csv'
# path_sc = '/Users/IEO5505/Desktop/example_mito/scratch_data/GBC_read_elements.tsv.gz'
# ncores = 8
# bulk_correction_treshold = 1
# umi_treshold = 5
# coverage_treshold = 50
# p_treshold = 0.001
# ratio_to_most_abundant_treshold = .3

# Import code
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from helpers import *


##


########################################################################

def main():


    # GBC correction, clone calling and cell assignment
    try: 
        # Custom workflow. Brings together filtering strategies from difference works.
        custom_workflow(
            path_bulk, 
            path_sample_map, 
            path_sc, 
            sample, 
            ncores=ncores,
            bulk_correction_treshold=bulk_correction_treshold,
            coverage_treshold=coverage_treshold,
            umi_treshold=umi_treshold, 
            p_treshold=p_treshold,
            ratio_to_most_abundant_treshold=ratio_to_most_abundant_treshold
        )
    except:
        raise Exception(
            f'''
            Some problem has been encoutered with the custom_workflow for the {sample} sample...
            '''
        )
    

    ##


########################################################################
    
# Run
if __name__ == '__main__':
    main()