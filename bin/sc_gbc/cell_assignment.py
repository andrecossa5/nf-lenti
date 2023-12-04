#!/usr/bin/python

# Cell assignment script

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='cell_assignment',
    description=
    """
    Script for clone calling and cell assignment.
    """
)


# Input
my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)

# Output
my_parser.add_argument(
    '--path_bulk', 
    type=str,
    default=None,
    help='Path to input bulk reference. Default: None.'
)

# treshold
my_parser.add_argument(
    '--path_sc',
    type=str,
    default=None,
    help='Path to input sc GBC reads elements. Default: None.'
)

# treshold
my_parser.add_argument(
    '--sample_map',
    type=str,
    default=None,
    help='Path to sample_map. Default: None.'
)

# treshold
my_parser.add_argument(
    '--method',
    type=str,
    default='unique_combos',
    help='Method for clone calling and cell assignment. Default: unique_combos.'
)

# treshold
my_parser.add_argument(
    '--ncores',
    type=int,
    default=8,
    help='n cores for pairwise distances calculation. Default: 8.'
)

# Spikeins
my_parser.add_argument(
    '--bulk_sc_treshold',
    type=int,
    default=1,
    help='''
    Hamming distance treshold to consider a sc GBC a "degenerate" sequence with respect 
    to a bulk reference sequence. Default: 1.
    '''
)

# Spikeins
my_parser.add_argument(
    '--umi_treshold',
    type=int,
    default=5,
    help='Min number of UMIs to consider a CB-GBC combination supported. Default: 5.'
)

# Spikeins
my_parser.add_argument(
    '--read_treshold',
    type=int,
    default=15,
    help='Min number of reads to consider a CB-GBC combination supported. Default: 15.'
)

# Spikeins
my_parser.add_argument(
    '--coverage_treshold',
    type=int,
    default=3,
    help='Min coverage (nUMIs / nreads) to consider a CB-GBC combination supported. Default: 3.'
)

# ratio_to_most_abundant_treshold
my_parser.add_argument(
    '--ratio_to_most_abundant_treshold',
    type=float,
    default=3,
    help='Min coverage (nUMIs / nreads) to consider a CB-GBC combination supported. Default: 3.'
)


##


# Parse arguments
args = my_parser.parse_args()
sample = args.sample
path_bulk = args.path_bulk
path_sample_map = args.sample_map
path_sc = args.path_sc
method = args.method
ncores = args.ncores
bulk_sc_treshold = args.bulk_sc_treshold
umi_treshold = args.umi_treshold
read_treshold = args.read_treshold
coverage_treshold = args.coverage_treshold
ratio_to_most_abundant_treshold = args.ratio_to_most_abundant_treshold


##

# path_sc = '/Users/IEO5505/Desktop/PD/tmp_sc/GBC_read_elements.tsv.gz'
# sample = 'X'
# path_bulk = '/Users/IEO5505/Desktop/PD/tmp_bulk/bulk_GBC_reference.csv'
# path_sample_map = None
# method = 'unique_combos'
# ncores = 8
# bulk_sc_treshold = 1
# umi_treshold = 5
# read_treshold = 30
# coverage_treshold = 10
# ratio_to_most_abundant_treshold = .3
# os.chdir('/Users/IEO5505/Desktop/PD/tmp_sc')

##


# Utils
def to_numeric(X):
    return np.select([X=='A', X=='T', X=='C', X=='G'], [1,2,3,4], default=0)


##


def call_clones(
    df_combos, sample=sample,
    read_treshold=read_treshold, umi_treshold=umi_treshold,
    coverage_treshold=coverage_treshold, 
    ratio_to_most_abundant_treshold=ratio_to_most_abundant_treshold
    ):
    """
    Call clones as unique GBC sequence sets, obtained from a
    table of CBC-UMI-GBC combinations, filtered properly. 
    """    
    
    # Supported and unsupported CBC-GBC combinations
    test = (df_combos['read_counts'] >= read_treshold) & \
        (df_combos['umi_counts'] >= umi_treshold) & \
        (df_combos['coverage'] >= coverage_treshold) & \
        (df_combos['ratio_to_most_abundant'] >= ratio_to_most_abundant_treshold)
    df_combos['status'] = np.where(test, 'supported', 'unsupported')

    # Pivot
    M = (
        df_combos
        .query('status=="supported"')
        .pivot_table(index='CBC', columns='GBC', values='umi_counts')
    )
    M[M.isna()] = 0

    # Get unique GBC sets
    cells_with_unique_combos = M.apply(lambda x: frozenset(M.columns[x>0]), axis=1)
    cells_with_unique_combos = cells_with_unique_combos.map(lambda x: ';'.join(x))
    clones = (
        cells_with_unique_combos
        .reset_index(drop=True)
        .value_counts(normalize=False)
        .to_frame('n cells').reset_index()
        .rename(columns={'index':'GBC_set'})
        .assign(clone=lambda x: [ f'C{i}_{sample}' for i in range(x.shape[0]) ])
        .set_index('clone')
    )

    return clones, cells_with_unique_combos


##


# Import code
import pandas as pd
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from itertools import chain
import matplotlib.pyplot as plt


##


########################################################################

def main():

    # Get the right reference sequences from bulk_GBC_reference
    bulk = pd.read_csv(
        os.path.join(path_bulk, 'summary', 'bulk_GBC_reference.csv'),
        index_col=0
    )

    if path_sample_map is not None:
        sample_map = pd.read_csv(path_sample_map, index_col=0)
        if sample in sample_map.index:
            ref = sample_map.loc[sample, 'reference']
            bulk = bulk.query('sample==@ref')
            assert bulk.shape[0]>0
            print(f'Found bulk GBC sequences for the {sample} sample, from ref {ref}.')
        else:
            raise KeyError(
                f'{sample} is not present in sample_map.csv index. Check errors.'
            )

    # Reverse complement and value_counts
    sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
    sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
    d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))
    sc = sc_df['GBC'].value_counts()

 
    ##

    
    # Map sc GBCs to bulk GBCs

    # Calculate hamming distance single-cell GBCs vs bulk reference GBCs
    sc_numeric = to_numeric(np.vstack(sc.index.map(lambda x: np.array(list(x)))))
    bulk_numeric = to_numeric(np.vstack(bulk.index.map(lambda x: np.array(list(x)))))
    D = pairwise_distances(
        sc_numeric, bulk_numeric, metric='hamming', n_jobs=int(ncores)
    ) * sc_numeric.shape[1]

    # Build a correction dict for sc GBCs at hamming distance <= bulk_sc_treshold
    # from a bulk one.
    d_corr = (
        sc.to_frame('read_count')
        .assign(
            correct_GBC=[ bulk.index[i] for i in D.argmin(axis=1) ],
            hamming=D.min(axis=1),
        )
        .query('hamming<=@bulk_sc_treshold')
        ['correct_GBC'].to_dict()
    )
    # Correct GBC sequences and remove not found ones
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: d_corr[x] if x in d_corr else 'not_found' )
    grouped = sc_df.query('GBC!="not_found"').groupby(['CBC', 'GBC'])

    
    ##


    # Compute sc CBC-GBCs combos and related stats
    read_counts = grouped.size().to_frame('read_counts')
    umi_counts = grouped['UMI'].nunique().to_frame('umi_counts')
    df_combos = (
        read_counts
        .join(umi_counts)
        .reset_index()
        .assign(coverage=lambda x: x['read_counts']/x['umi_counts'])
    )
    df_combos = df_combos.join(
        df_combos
        .groupby('CBC')
        .apply(lambda x: x['umi_counts'] / x['umi_counts'].max())
        .droplevel(0)
        .to_frame('ratio_to_most_abundant')
    )
    df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')

        
    ##


    # Cell assignment

    # Call clones
    clones, cells_with_unique_combos = call_clones(
        df_combos, sample=sample,
        read_treshold=read_treshold, 
        umi_treshold=umi_treshold, 
        coverage_treshold=coverage_treshold, 
        ratio_to_most_abundant_treshold=ratio_to_most_abundant_treshold
    )
    combos = clones['GBC_set'].map(lambda x: x.split(';')).to_list()
    unique_GBC = np.unique(list(chain.from_iterable(combos)))

    # Calculate GBC sequences occurrences across these sets
    occurrences = np.zeros((unique_GBC.size, unique_GBC.size))
    for i,x in enumerate(unique_GBC):
        for j,y in enumerate(unique_GBC):
            if any([ x in c and y in c for c in combos ]):
                occurrences[i,j] = 1
    occurrences = pd.DataFrame(occurrences, index=unique_GBC, columns=unique_GBC)

    # Assign cells to a GBC_set (i.e., clone)
    cells = (
        cells_with_unique_combos
        .map(lambda x: clones.index[clones['GBC_set']==x][0])
        .to_frame('clone')
    )


    ##
    

    # Save clones, cells and GBC co-occurrences
    clones.to_csv('clones_summary_table.csv')
    cells.to_csv('cells_summary_table.csv')
    occurrences.to_csv('occurrences.csv')

    # Combinations support plot
    fig, ax = plt.subplots(figsize=(6,5))
    x = np.log10(df_combos['read_counts'])
    y = np.log10(df_combos['umi_counts'])
    ax.plot(x[df_combos['status'] == 'supported'], 
            y[df_combos['status'] == 'supported'], '.', 
            label='assigned', color='blue', markersize=.5, zorder=10)
    ax.plot(x[df_combos['status'] == 'unsupported'],
            y[df_combos['status'] == 'unsupported'], '.', 
            label='not-assigned', color='grey', markersize=.3)
    ax.set(
        title='CBC-GBC combination status', 
        xlabel='log10_read_counts', 
        ylabel='log10_umi_counts'
    )
    ax.legend()
    fig.savefig('CBC_GBC_combo_status.png')


    ##


########################################################################
    
# Run
if __name__ == '__main__':
    main()