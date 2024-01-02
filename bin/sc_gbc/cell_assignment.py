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
    default=3,
    help='''
    Hamming distance treshold to consider a sc GBC a "degenerate" sequence with respect 
    to a bulk reference sequence. Default: 3.
    '''
)

# Spikeins
my_parser.add_argument(
    '--sc_correction_treshold',
    type=int,
    default=3,
    help='''
    Hamming distance treshold to consider a sc GBC a "degenerate" sequence with respect 
    to another in single-cell data. Default: 3.
    '''
)

# correction_type
my_parser.add_argument(
    '--correction_type',
    type=str,
    default='reference-free',
    help='Method to correct observed GBC sequences. Default: reference-free.'
)

# filtering_method
my_parser.add_argument(
    '--filtering_method',
    type=str,
    default='GMM',
    help='Method to discard noisy UMIs. Default: GMM.'
)

# Spikeins
my_parser.add_argument(
    '--coverage_treshold',
    type=int,
    default=10,
    help='''
        Min number of reads to consider a CBC-GBC-UMI combination "true". Default: 10. This
        is used only if the --filtering method is not GMM (default) or medstd.
        '''
)

# Spikeins
my_parser.add_argument(
    '--umi_treshold',
    type=int,
    default=5,
    help='Min number of UMIs to consider a CBC-GBC combination supported. Default: 5.'
)

my_parser.add_argument(
    '--from_subset', 
    action='store_true',
    help='If the qc needs to be done of a subsetted .h5ad matrix. Default: False.'
)

# p_treshold
my_parser.add_argument(
    '--p_treshold',
    type=float,
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
sc_correction_treshold = args.sc_correction_treshold
correction_type = args.correction_type
filtering_method = args.filtering_method
coverage_treshold = args.coverage_treshold
umi_treshold = args.umi_treshold
p_treshold = args.p_treshold
ratio_to_most_abundant_treshold = args.ratio_to_most_abundant_treshold

# sample = 'AC_AC_mets_3'
# path_bulk = '/Users/IEO5505/Desktop/example_mito/scratch_data/'
# path_sample_map = '/Users/IEO5505/Desktop/example_mito/scratch_data/sample_map.csv'
# path_sc = '/Users/IEO5505/Desktop/example_mito/scratch_data/GBC_read_elements.tsv.gz'
# ncores = 8
# bulk_correction_treshold = 3
# sc_correction_treshold = 3
# correction_type = 'reference-free'
# filtering_method = 'GMM'
# coverage_treshold = 10
# umi_treshold = 5
# p_treshold = 0.01
# ratio_to_most_abundant_treshold = .3


# Import code
# sys.path.append('/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/sc_gbc')
# os.chdir('/Users/IEO5505/Desktop/example_mito/scratch_results')
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from helpers import *


##


########################################################################

def main():

    try: 

        """
        Custom workflow. 
        Brings together filtering and correction strategies from difference works:
        * Adamson et al., Dixit et al., Cell 2016
        * Weinreb et al., Science 2020
        * Roda and Cossa et al., Cancer Research 2023 
        * Nadalin et al., pre-print on biorxiv 2023
        """

        cell_assignment_workflow(
            path_bulk, 
            path_sample_map, 
            path_sc, 
            sample, 
            correction_type=correction_type, 
            sc_correction_treshold=sc_correction_treshold,
            bulk_correction_treshold=bulk_correction_treshold,
            ncores=ncores,
            filtering_method=filtering_method,
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