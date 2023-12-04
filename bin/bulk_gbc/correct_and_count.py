#!/usr/bin/python

# Correct and count genomic GBC

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='correct_and_count',
    description=
    """
    Script for GBC filtering, correction and counting.
    """
)


# Input
my_parser.add_argument(
    '-i', 
    '--input', 
    type=str,
    default=None,
    help='''
        The path to extracted GBCs file. Must be a .csv file, 
        with one column. Default: . .
        '''
)

# Output
my_parser.add_argument(
    '-o', 
    '--output', 
    type=str,
    default=os.getcwd(),
    help='Output folder path. Default: . .'
)

# treshold
my_parser.add_argument(
    '-t', 
    '--threshold',
    type=int,
    default=3,
    help='Maximum Hamming distance for which sequences can be clustered together.'
)

# treshold
my_parser.add_argument(
    '--method',
    type=str,
    default='directional',
    help='UMI-tools method. Default: directional.'
)

# treshold
my_parser.add_argument(
    '--min_n_reads',
    type=int,
    default=1000,
    help='''
        min_n_reads that a "correct" GBCs must have to be considered to entry 
        the final whitelist. Default: 1000.'
        '''
)

# Spikeins
my_parser.add_argument(
    '--spikeins',
    type=str,
    default=None,
    help='Path to spikeins table. Default: None.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_i = args.input
path_o = args.output
method = args.method
threshold = args.threshold
min_n_reads = args.min_n_reads
spikeins = args.spikeins


##


# Import code
import numpy as np
import pandas as pd
import dask.dataframe as dd
from umi_tools import UMIClusterer


##


########################################################################

def main():

    # Remove short GBCs (len == 18bp), retain only GBC >1 read
    GBCs = dd.read_csv(path_i, header=None, sep='\t')[0]
    is_18bp = GBCs.map(lambda x: len(x) == 18)
    GBCs = GBCs.loc[is_18bp]
    GBC_counts = GBCs.value_counts().compute().astype(np.int32)
    # Save all raw counts
    GBC_counts.to_frame('read_counts').to_csv(os.path.join(path_o, 'GBC_raw_counts.csv.gz'))

    # Cluster unique GBCs
    GBC_counts = GBC_counts.loc[lambda x: x>1].copy() # Only >1 reads
    GBC_counts.index = GBC_counts.index.map(lambda x : x.encode('UTF-8'))
    clusterer = UMIClusterer(cluster_method=method)
    groups = clusterer(GBC_counts.to_dict(), threshold=threshold)
    GBC_counts.index = GBC_counts.index.map(lambda x : x.decode('UTF-8'))

    # Get non-singletons and reformat
    groups = [ g for g in groups if len(g)>1 ]
    df_l = [
        pd.DataFrame(
            [ ( g[0].decode("utf-8"), x.decode("utf-8") ) for x in g[1:] ],
            columns=['correct', 'degenerated']
        )
        for g in groups
    ]
    del groups                  # Save memory
    df = pd.concat(df_l)
    del df_l                    # Save memory

    # Compute correction_df
    df['n_reads_correct'] = GBC_counts.loc[df['correct']].values
    df['n_reads_degenerated'] = GBC_counts.loc[df['degenerated']].values
    corrected_counts = (
        df.groupby('correct')
        .apply(
            lambda x: x['n_reads_correct'].unique()[0] + x['n_reads_degenerated'].sum()
        )
    )
    df.to_csv(os.path.join(path_o, 'correction_df.csv'), index=False)

    # Remove spikeins
    if spikeins is not None:
        spikes = pd.read_csv(spikeins, index_col=0).index
        GBC_counts = GBC_counts.loc[~GBC_counts.index.isin(spikes)]

    # Get final GBC pool
    GBC_counts = GBC_counts.to_frame('read_count')
    correct = df['correct'].unique()
    degenerated = df['degenerated'].unique()
    non_clustered = GBC_counts.index[
        ~GBC_counts.index.isin(np.concatenate([degenerated, correct]))
    ]
    final_pool = GBC_counts.index[
        GBC_counts.index.isin(correct) | (GBC_counts.index.isin(non_clustered))
    ]
    assert final_pool.size == (correct.size + non_clustered.size)

    # Final counts
    final_counts = pd.DataFrame(0, index=final_pool, columns=['read_count'])
    final_counts.loc[non_clustered, 'read_count'] = GBC_counts.loc[non_clustered, 'read_count']
    final_counts.loc[correct, 'read_count'] = corrected_counts.loc[correct]

    # Filter min_n_reads, and save
    (
        final_counts
        .query('read_count>=@min_n_reads')
        .to_csv(os.path.join(path_o, 'GBC_counts_corrected.csv'))
    )


##


########################################################################

# Run
if __name__ == '__main__':
    main()