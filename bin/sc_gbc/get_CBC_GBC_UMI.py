#!/usr/bin/python

import os
import sys
import pysam
import pandas as pd
from scipy.spatial.distance import hamming


##


# Paths
path_bam = sys.argv[1]
path_filtered = sys.argv[2]
anchor = sys.argv[3]
treshold = sys.argv[4]


# path_ = '/Users/IEO5505/Desktop/example_mito/scratch_data'
# os.listdir(path_)
# path_bam = os.path.join(path_, 'lentibam.bam')
# path_filtered = os.path.join(path_, 'barcodes.tsv.gz')
# anchor = 'TAGCAAACTGGGGCACAAGCTTAATTAAGAATT'
# treshold = 1


##


# Read Solo-filtered CBCs 
solo_CBCs = pd.read_csv(
    os.path.join(path_filtered, 'barcodes.tsv.gz'), 
    header=None, index_col=0
)

# Read bam and parse records
bam_in = pysam.AlignmentFile(path_bam, 'rb')
el = open('GBC_read_elements.tsv', 'w')

for r in bam_in:
    name = r.query_name
    a_ = r.seq[:33]
    gbc = r.seq[33:33+18]
    umi = r.get_tag('UB')
    cbc = r.get_tag('CB')
    h = hamming(list(a_), list(anchor)) * 33
    if h <= int(treshold) and cbc in solo_CBCs.index:
        el.write(f'@{name}\t{cbc}\t{umi}\t{gbc}\n')

# Close streams
bam_in.close()
el.close()


##
