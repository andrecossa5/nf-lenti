#!/usr/bin/python

# maker_AFM_matrix script

import os
import argparse


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='make AFM matrix',
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

my_parser.add_argument(
    '--scLT_system', 
    type=str,
    default=None,
    help='scLT system. System for scLT. Default: MAESTER.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_ch_matrix = args.path_ch_matrix
path_meta = args.path_meta
sample = args.sample
scLT_system = args.scLT_system
pp_method = args.pp_method


##


def main():
    
    import mito as mt

    afm = mt.io.make_afm(
        path_ch_matrix, 
        path_meta=path_meta, 
        sample=sample, 
        pp_method=pp_method, 
        scLT_system=scLT_system
    )
    
    if os.path.isdir(path_ch_matrix):
        afm.write(os.path.join(path_ch_matrix, 'afm.h5ad'))
    else:
        afm.write(os.path.join(os.path.dirname(path_ch_matrix), 'afm.h5ad'))


if __name__ == '__main__':
    main()