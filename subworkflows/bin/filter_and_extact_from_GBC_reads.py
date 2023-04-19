#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
import dnaio 


##

# Paths
path_R1_in = sys.argv[1]
path_R2_in = sys.argv[2]
path_filtered = sys.argv[3]

# Read Solo-filtered CBCs 
solo_CBCs = pd.read_csv(path_filtered + '/barcodes.tsv.gz', header=None, index_col=0)

# Open stream fqs and .txt files
fqs_in = dnaio.open(path_R1_in, file2=path_R2_in, mode='r')
fqs_out = dnaio.open('filtered_R1.fq.gz', file2='filtered_R2.fq.gz', mode='w')
CBCs = open('CBCs_by_read.tsv', 'w')
UMIs = open('UMIs_by_read.tsv', 'w')
GBCs = open('GBCs_by_read.tsv', 'w')

# Parse and write, with a single for loop
for r1,r2 in fqs_in:
    name = r1.name
    cbc = r1.sequence[:16]
    umi = r1.sequence[16:16+12]
    gbc = r2.sequence[33:33+18]
    if cbc in solo_CBCs.index:
        CBCs.write(f'@{name}\t{cbc}\n')
        UMIs.write(f'@{name}\t{umi}\n')
        GBCs.write(f'@{name}\t{gbc}\n')
        fqs_out.write(r1,r2)

# Close streams
fqs_in.close()
fqs_out.close()
CBCs.close()
UMIs.close()
GBCs.close()


##





