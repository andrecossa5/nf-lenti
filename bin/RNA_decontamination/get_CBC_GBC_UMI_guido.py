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
    '--grepped_path', 
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
path_data = args.path_data
path_barcodes = args.path_barcodes
path_barcodes = os.path.join(path_data, 'barcodes.tsv.gz')
path_grepped = os.path.join(path_data, 'grepped.txt.gz')


##


# Parameters
anchor = 'TAGCAAACTGGGGCACAAGCTTAATTAAGAATT'
treshold = 1

##


# Read Solo-filtered CBCs 
solo_CBCs = pd.read_csv(
    path_barcodes, 
    header=None, index_col=0
)

# From unsorted bam, grepped
file_in = gzip.open(path_grepped, 'rb')
next(file_in)
file_out = gzip.open(f'/hpcnfs/home/ieo6943/results/{type}/{sample}/GBC_read_elements.tsv.gz', 'wb')
anchor = 'TAGCAAACTGGGGCACAAGCTTAATTAAGAATT'


for line in file_in:
    fields = line.decode('utf-8').strip().split('\t')
    #print(f'fields:{fields}')
    #print('len fields=', len(fields))
    if len(fields)>=4:
        
        read_name = fields[0]
        #print(f'read_name:{read_name}')
        
        cr = fields[1]
        #print(f'cr:{cr}')
        ur = fields[2]
        #print(f'ur:{ur}')
        a_ = fields[3][:33]
        #print(f'a_:{a_}')
        gbc = fields[3][33:33+18]
        #print(f'gbc:{gbc}')
        h = hamming(list(a_), list(anchor)) * 33

        if h <= int(treshold) and cr in solo_CBCs.index:
            r = f'@{read_name}\t{cr}\t{ur}\t{gbc}\n'
            #file_out.write(r.encode('utf-8'))



file_in.close()
file_out.close()