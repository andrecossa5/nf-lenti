#!/usr/bin/python

# maker_AFM_matrix script

import os
import argparse


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='maker_AFM_matrix',
    description=
    """
    Prepare AFM.
    """
)

# Add arguments
my_parser.add_argument(
    '--path_ch_matrix', 
    type=str,
    default=None,
    help='Path to ch_matrix. Default: None'
)

my_parser.add_argument(
    '--path_meta', 
    type=str,
    default=None,
    help='Path to cells_meta.csv. Default: None'
)

my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None'
)

my_parser.add_argument(
    '--pp_method', 
    type=str,
    default=None,
    help='Pipeline used for preprocessing. Default: mito_preprocessing'
)


##


# Parse arguments
args = my_parser.parse_args()
path_ch_matrix = args.path_ch_matrix
path_meta = args.path_meta
sample = args.sample
pp_method = args.pp_method


##


def main():
    
    from mito_utils.make_afm import make_afm as make

    afm = make(path_ch_matrix, path_meta, sample=sample, pp_method=pp_method)
    afm.write(os.path.join(path_ch_matrix, 'afm.h5ad'))


if __name__ == '__main__':
    main()