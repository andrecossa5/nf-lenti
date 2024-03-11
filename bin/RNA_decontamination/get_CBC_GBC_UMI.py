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
#path_grepped = sys.argv[1]
#path_filtered = sys.argv[2]
#anchor = sys.argv[3]
#treshold = sys.argv[4]


path_data = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data'
# path_bam = os.path.join(path_data, 'lentibam.bam')
path_filtered = os.path.join(path_data, 'barcodes.tsv.gz')
anchor = 'TAGCAAACTGGGGCACAAGCTTAATTAAGAATT'
treshold = 1
path_grepped = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/grepped.txt.gz'


##




# Read Solo-filtered CBCs 
#solo_CBCs = pd.read_csv(
 #   os.path.join(path_filtered, 'barcodes.tsv.gz'), 
  #  header=None, index_col=0
#)
solo_CBCs = pd.read_csv(
    os.path.join(path_data, 'barcodes.tsv.gz'), 
    header=None, index_col=0
)

# From unsorted bam, grepped
file_in = gzip.open(path_grepped, 'rb')
file_out_not = gzip.open('GBC_read_elements_not_1M.tsv.gz', 'wb')
file_out = gzip.open('GBC_read_elements_in_1M.tsv.gz', 'wb')



tot_lines = 0
tot_fields = 0
tot_h = 0
tot_h_in = 0
for line in file_in:
    tot_lines += 1
    fields = line.decode('utf-8').strip().split('\t')
    if len(fields)>=10:
        tot_fields += 1
        read_name = fields[0]
        cr = fields[-2].split(':')[-1]
        ur = fields[-1].split(':')[-1]
        read_sequence = fields[9]
        a_ = read_sequence[:33]
        gbc = read_sequence[33:33+18]
        print('gbc = ', gbc)
        h = hamming(list(a_), list(anchor)) * 33
        print('h =', h)
        if h <= int(treshold) and cr not in solo_CBCs.index:
            tot_h += 1
            r = f'@{read_name}\t{cr}\t{ur}\t{gbc}\n'
            file_out_not.write(r.encode('utf-8'))

        elif h <= int(treshold) and cr in solo_CBCs.index:
            tot_h_in += 1
            r = f'@{read_name}\t{cr}\t{ur}\t{gbc}\n'
            file_out.write(r.encode('utf-8'))

    print('tot_lines = ', tot_lines)
    if tot_lines >1000000: break


file_in.close()
file_out.close()
file_out_not.close()

