#!/usr/bin/python

import os
import sys
import gzip
from io import TextIOWrapper
import pandas as pd
from scipy.spatial.distance import hamming
import io

##


# Paths
path_main = '/Users/ieo6943/Documents/data/8_clones'
#path_grepped = sys.argv[1]
#path_filtered = sys.argv[2]
#anchor = sys.argv[3]
#treshold = sys.argv[4]
path_data = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data'
# path_bam = os.path.join(path_data, 'lentibam.bam')
path_filtered = os.path.join(path_data, 'barcodes.tsv.gz')
path_grepped = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/grepped.txt.gz'


##


#Parameters
anchor = 'TAGCAAACTGGGGCACAAGCTTAATTAAGAATT'
treshold = 1

##


# Read Solo-filtered CBCs 
solo_CBCs = pd.read_csv(
    os.path.join(path_data, 'barcodes.tsv.gz'), 
    header=None, index_col=0
)

# From unsorted bam, grepped
file_in = gzip.open(path_grepped, 'rb')
file_out = gzip.open('GBC_read_elements.tsv.gz', 'wb')


for line in file_in:
    fields = line.decode('utf-8').strip().split('\t')
    if len(fields)>=10:
        
        read_name = fields[0]
        cr = fields[-2].split(':')[-1]
        ur = fields[-1].split(':')[-1]
        read_sequence = fields[9]
        a_ = read_sequence[:33]
        gbc = read_sequence[33:33+18]
        h = hamming(list(a_), list(anchor)) * 33

        if h <= int(treshold) and cr not in solo_CBCs.index:

            r = f'@{read_name}\t{cr}\t{ur}\t{gbc}\n'
            file_out.write(r.encode('utf-8'))


file_in.close()
file_out.close()

