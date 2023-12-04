#!/usr/bin/python

import os
import sys
import pandas as pd
import dnaio 
from scipy.spatial.distance import hamming


##


# Paths
path_R1_in = sys.argv[1]
path_R2_in = sys.argv[2]
path_filtered = sys.argv[3]
anchor = sys.argv[4]
treshold = sys.argv[5]


##


# Read Solo-filtered CBCs 
solo_CBCs = pd.read_csv(
    os.path.join(path_filtered, 'barcodes.tsv.gz'), 
    header=None, index_col=0
)

# Open stream fqs and .txt files
fqs_in = dnaio.open(path_R1_in, file2=path_R2_in, mode='r')
el = open('GBC_read_elements.tsv', 'w')

# Parse danio stream and write with a single for loop
for r1,r2 in fqs_in:
    name = r1.name
    cbc = r1.sequence[:16]
    umi = r1.sequence[16:16+12]
    a_ = r2.sequence[:33]
    gbc = r2.sequence[33:33+18]
    h = hamming(list(a_), list(anchor)) * 33
    if h <= int(treshold) and cbc in solo_CBCs.index:
        el.write(f'@{name}\t{cbc}\t{umi}\t{gbc}\n')

# Close streams
fqs_in.close()
el.close()


##





