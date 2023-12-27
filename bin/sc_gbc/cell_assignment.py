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
    '--sc_correction_treshold',
    type=int,
    default=3,
    help='''
    Hamming distance treshold to consider a sc GBC a "degenerate" sequence with respect 
    to another found in single-cell data. Default: 3.
    '''
)

# Spikeins
my_parser.add_argument(
    '--umi_treshold',
    type=int,
    default=5,
    help='Min number of UMIs to consider a CB-GBC combination supported. Default: 5.'
)

# Spikeins
my_parser.add_argument(
    '--read_treshold',
    type=int,
    default=30,
    help='Min number of reads to consider a CB-GBC combination supported. Default: 30.'
)

# Spikeins
my_parser.add_argument(
    '--coverage_treshold',
    type=int,
    default=10,
    help='Min coverage (nUMIs / nreads) to consider a CB-GBC combination supported. Default: 10.'
)

# ratio_to_most_abundant_treshold
my_parser.add_argument(
    '--ratio_to_most_abundant_treshold',
    type=float,
    default=.5,
    help='Min coverage (nUMIs / nreads) to consider a CB-GBC combination supported. Default: .5.'
)

# Method
my_parser.add_argument(
    '--method',
    type=str,
    default='CR_2023',
    help='Method for cell assignment. Default: CR_2023.'
)


##


# Parse arguments
args = my_parser.parse_args()
sample = args.sample
path_bulk = args.path_bulk
path_sample_map = args.sample_map
path_sc = args.path_sc
method = args.method
ncores = args.ncores
bulk_correction_treshold = args.bulk_correction_treshold
sc_correction_treshold = args.sc_correction_treshold
read_treshold = args.read_treshold
umi_treshold = args.umi_treshold
coverage_treshold = args.coverage_treshold
ratio_to_most_abundant_treshold = args.ratio_to_most_abundant_treshold

# sample = 'AC_NT_mets_2'
# path_bulk = '/Users/IEO5505/Desktop/example_mito/scratch_data/bulk_GBC_reference.csv'
# path_sample_map = '/Users/IEO5505/Desktop/example_mito/scratch_data/sample_map.csv'
# path_sc = '/Users/IEO5505/Desktop/example_mito/scratch_data/GBC_read_elements.tsv.gz'
# method = 'LARRY_2020'
# ncores = 8
# bulk_correction_treshold = 1
# sc_correction_treshold = 3
# read_treshold = 30
# umi_treshold = 5
# coverage_treshold = 10
# ratio_to_most_abundant_treshold = .3

# Import code
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from CR_2023 import *
from LARRY_2020 import *


##


########################################################################

def main():


    # GBC correction, clone calling and cell assignment
    if method == 'CR_2023':

        try: 
            # Roda and Cossa et al. 2023, Cancer Research (Dixit and Adamson 2016, Cell)
            CR_2023_workflow(
                path_bulk, 
                path_sample_map, 
                path_sc, 
                sample, 
                ncores, 
                bulk_correction_treshold=bulk_correction_treshold,
                read_treshold=read_treshold, 
                umi_treshold=umi_treshold, 
                coverage_treshold=coverage_treshold,
                ratio_to_most_abundant_treshold=ratio_to_most_abundant_treshold
            )
        except:
            raise Exception(
                f'''
                Some problem has been encoutered with the CR_2023_workflow for the {sample} sample...
                '''
            )
    
    elif method == 'LARRY_2020':

        try: 
            # Weinreb et al. 2020, Science
            LARRY_2020_workflow(
                path_bulk, 
                path_sample_map, 
                path_sc, 
                sample, 
                sc_correction_treshold=sc_correction_treshold,
                read_treshold=read_treshold, 
                umi_treshold=umi_treshold, 
                coverage_treshold=coverage_treshold
            )
        except:
            raise Exception(
                f'''
                Some problem has been encoutered with the LARRY_2020_workflow for the {sample} sample...
                '''
            )

    else:
        raise ValueError(f'The {method} method is not supported...')
    

    ##


########################################################################
    
# Run
if __name__ == '__main__':
    main()