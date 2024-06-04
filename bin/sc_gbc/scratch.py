#!/usr/bin/python

import pysam
import numpy as np
import pandas as pd
from mito_utils.plotting_base import *
import sys


##


#ARG
bam_file = sys.argv[1]


##


def main(bam_file):

    # Open files
    bam_in = pysam.AlignmentFile(bam_file, "rb", check_sq=False)

    L = []
    for read in bam_in:

        # Consensus read info
        umi = read.get_tag('MI:')  
        sequence = np.array(list(read.query_sequence))
        n_consensus_read = read.get_tag('cD:')
        qualities = np.array(read.query_qualities) 
        mean_quality = qualities.mean()
        mean_consensus_error = read.get_tag('cE:')
        alignment_start = read.reference_start

        # Stats
        GBC = sequence[33:33+18]
        read_N_perc = np.sum(sequence=='N')/sequence.size
        GBC_N_perc = np.sum(GBC=='N')/18

        read_status = 'supported' 

        # Add to sample df
        final_seq = ''.join(sequence)
        GBC = ''.join(GBC)

        d = { 
            'UMI':umi, 'GBC':GBC,
            'read':final_seq, 'mean_consensus_error':mean_consensus_error, 
            'mean_quality':mean_quality, 'read_N_perc':read_N_perc, 
            'GBC_N_perc':GBC_N_perc, 'read_status':read_status,
            'alignment_start':alignment_start, 'n_consensus_read':n_consensus_read

        }
        L.append(d)

    # Close files and create dataframe results
    bam_in.close()
    df = pd.DataFrame(L)
    df.to_csv('filtered_consensus_read.csv')


##


if __name__ == '__main__':
    main()