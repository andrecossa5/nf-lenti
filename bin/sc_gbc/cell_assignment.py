#!/usr/bin/python

"""
Cell assignment script.
"""
 

##


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
    '--umi_treshold',
    type=int,
    default=5,
    help='Min number of UMIs to consider a CBC-GBC combination supported. Default: 5.'
)

# p_treshold
my_parser.add_argument(
    '--p_treshold',
    type=float,
    default=.5,
    help='Max p_poisson treshold to consider a CBC-GBC combination supported. Default: .001.'
)

# ratio_to_most_abundant_treshold
my_parser.add_argument(
    '--max_ratio_treshold',
    type=float,
    default=.8,
    help='Min ratio between a GBC nUMIs and the most abundant (nUMIs) GBC found for a given CBC. Default: .8.'
)

# Normalized abundance
my_parser.add_argument(
    '--normalized_abundance_treshold',
    type=float,
    default=.8,
    help='Min abundance (nUMIs fraction within a cell) of a CBC-GBC combination. Default: .8.'
)

# sample_params
my_parser.add_argument(
    '--sample_params',
    type=str,
    default="NULL",
    help='Path to sample_specific filtering parameters. Default: NULL.'
)


##


# Parse arguments
args = my_parser.parse_args()
sample = args.sample
path_bulk = args.path_bulk
path_sample_map = args.sample_map
path_sc = args.path_sc
bulk_correction_treshold = args.bulk_correction_treshold
umi_treshold = args.umi_treshold
p_treshold = args.p_treshold
max_ratio_treshold = args.max_ratio_treshold
normalized_abundance_treshold = args.normalized_abundance_treshold
sample_params = args.sample_params


# Import code
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from helpers import *


##


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

        if sample_params != "NULL":
            params = pd.read_csv(sample_params, index_col=0)
            params = params.loc[sample].to_dict()
        else:
            params = None
            
        cell_assignment_workflow(
            path_sc, 
            sample=sample, 
            path_bulk=path_bulk, 
            path_sample_map=path_sample_map, 
            umi_treshold=umi_treshold, 
            p_treshold=p_treshold,
            max_ratio_treshold=max_ratio_treshold,
            normalized_abundance_treshold=normalized_abundance_treshold,
            sample_params=params,
            bulk_correction_treshold=bulk_correction_treshold
        )

    except:

        raise Exception(
            f'''
            Some problem has been encoutered with the custom_workflow for the {sample} sample...
            '''
        )
    

    ##


# Run
if __name__ == '__main__':
    main()
