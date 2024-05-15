#!/usr/bin/python

import gzip
import os
import pandas as pd
from scipy.spatial.distance import hamming
import argparse

##


my_parser = argparse.ArgumentParser(
    prog='Get CBC-UMI-GBC GUido versione',
    description=
    """
    Get CBC-UMI-GBC GUido versione.
    """
)

my_parser.add_argument(
    '--grepped_path', 
    type=str,
    default=None,
    help='path of the grepped file'
)

my_parser.add_argument(
    '--barcodes_path', 
    type=str,
    default=None,
    help='path of the barcodes'
)

my_parser.add_argument(
    '--type', 
    type=str,
    default=None,
    help='type of experiment: lentiviral, GEX, mitocondrial'
)

my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Datasample: AML_clones, MDA_clones, ...'
)

# Paths
args = my_parser.parse_args()
type =args.type
sample = args.sample
grepped_path = args.grepped_path
#path_barcodes = args.barcodes_path
#path_barcodes = os.path.join(path_barcodes, 'barcodes.tsv.gz')
path_data = os.path.join(grepped_path, 'grepped.txt.gz')


##


#type = 'sc'
#sample = 'AML_clones'
#grepped_path = '/Users/ieo6943/Documents/data/AML_clones/'
#path_barcodes = os.path.join(grepped_path, 'barcodes.tsv.gz')
#path_data = '/Users/ieo6943/Documents/data/AML_clones/grepped_scratch.txt.gz'
#path_file_GBC = '/Users/ieo6943/Documents/data/AML_clones/GBC_read_elements.tsv.gz'

##


# Parameters
anchor = 'TAGCAAACTGGGGCACAAGCTTAATTAAGAATT'
treshold = 1

##


# Read Solo-filtered CBCs 
#solo_CBCs = pd.read_csv(
#    path_barcodes, 
#    header=None, index_col=0
#)

# From unsorted bam, grepped
file_in = gzip.open(path_data, 'rb')
next(file_in)
file_out = gzip.open(f'/hpcnfs/home/ieo6943/results/{type}/{sample}/GBC_read_elements_scratch.tsv.gz', 'wb')
#file_out = gzip.open(path_file_GBC, 'wb')
#i=0
for line in file_in:
    fields = line.decode('utf-8').strip().split('\t')
    if len(fields)>=4: 
        read_name = fields[0]
        cr = fields[1]
        ur = fields[2]
        a_ = fields[3][:33]
        gbc = fields[3][33:33+18]
        h = hamming(list(a_), list(anchor)) * 33
        if h <= int(treshold):
            r = f'@{read_name}\t{cr}\t{ur}\t{gbc}\n'
            file_out.write(r.encode('utf-8'))




file_in.close()
file_out.close()

