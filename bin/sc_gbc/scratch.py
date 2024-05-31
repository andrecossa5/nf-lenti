#!/usr/bin/python

import pysam
import numpy as np
import pandas as pd
from mito_utils.plotting_base import *



##


# Paths
bam_file = '/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/consensus.bam'
output_bam = '/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/filter_local.bam'

# Params
min_quality = 30
base_consensus_error = 0.1
read_max_consensus_error = 0.1
read_max_N = 0.2
GBC_max_N = 0
consensus_filtering_mode = 'GBC'


##


# Open files
bam_in = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

L = []
for read in bam_in:

    # Consensus read info
    umi = read.get_tag('MI:')  
    sequence = np.array(list(read.query_sequence))
    n_consensus_read = read.get_tag('cM:')
    qualities = np.array(read.query_qualities) 
    mean_quality = qualities.mean()
    base_consensus_error = np.array(read.get_tag('ce:'))/n_consensus_read
    mean_consensus_error = base_consensus_error.mean()

    # Mask
    bad_bases_test = (base_consensus_error>base_consensus_error) | (qualities<min_quality)
    sequence[bad_bases_test] = 'N'

    # Stats
    GBC = sequence[33:33+18]
    read_N_perc = np.sum(sequence=='N')/sequence.size
    GBC_N_perc = np.sum(GBC=='N')/18
    
    # Test
    read_status = 'non_supported'
    consensus_read_test = mean_consensus_error <= read_max_consensus_error and \
                          mean_quality >= min_quality and \
                          read_N_perc <= read_max_N
    if consensus_read_test:
        if consensus_filtering_mode:
            if GBC_max_N <= GBC_max_N:
                read_status = 'supported'
        read_status = 'supported'
    
    # Add to sample df
    final_seq = ''.join(sequence)
    GBC = ''.join(GBC)
    
    # Write if supported read
    if read_status == 'supported':
        read.query_sequence = final_seq
        bam_out.write(read)

    # Append to L
    d = { 
        'UMI':umi, 'GBC':GBC,
        'read':final_seq, 'mean_consensus_error':mean_consensus_error, 
        'mean_quality':mean_quality, 'read_N_perc':read_N_perc, 
        'GBC_N_perc':GBC_N_perc, 'read_status':read_status
    }
    L.append(d)


##


# Close files and create dataframe results
bam_in.close()
bam_out.close()
df = pd.DataFrame(L)


##


df.query('read_status=="supported"').groupby('UMI')['GBC'].nunique().describe()


df['UMI'].value_counts().head(50)

np.sum(df['UMI'].value_counts()==1)
df['UMI'].value_counts().size

