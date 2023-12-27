"""
Utils for clone calling and cell assignment as in the LARRY paper, weinberg et al., 2020.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting_utils._utils import Timer
from plotting_utils._plotting_base import *
from itertools import chain


##


# Utils
def hamming(bc1, bc2):
    h = np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])
    return h


##


def get_combos_LARRY2020(path_bulk, path_sample_map, sample, path_sc, 
                        sc_correction_treshold=3, coverage_treshold=10
    ):

    # Get the right reference sequences from bulk_GBC_reference 
    if path_bulk is not None:
        bulk = pd.read_csv(
            os.path.join(path_bulk),#, 'summary', 'bulk_GBC_reference.csv'),
            index_col=0
        )
    else:
        raise ValueError('Specify a valid path_bulk path!')

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
    else:
        raise ValueError('Specify a valid path_sample_map path!')


    ##


    # Reverse complement GBCs from reads_elements
    sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
    sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
    d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))  

    # Count CBC-GBC-UMI combinations
    counts = sc_df.groupby(['CBC', 'GBC', 'UMI']).size().reset_index(name='count')

    # Viz distribution n reads
    fig, ax = plt.subplots(figsize=(5,5))
    counts['log'] = np.log(counts['count'])/np.log(10)
    hist(counts, 'log', c='k', ax=ax, n=50)
    format_ax(
        xticks=np.logspace(0,4,5), ax=ax, xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
        title='CBC-GBC-UMI combination n reads distribution'
    )
    ax.set_yscale('log')
    ax.axvline(np.log(10)/np.log(10),c='r')
    fig.tight_layout()

    # Filter CBC-GBC-UMI combos for min n reads
    filtered = counts.query('count>=@coverage_treshold')                                                      
    sc_df = pd.merge(
        sc_df, 
        filtered[['CBC', 'GBC', 'UMI']], 
        on=['CBC', 'GBC', 'UMI'], 
        how='inner'
    )                                                                                                               

    ##

    # Build correction map
    all_gfp_bcs = sorted(set([x for x in sc_df['GBC'].unique()]))
    good_gfp_bcs = []
    bc_map = {}
    for i,bc1 in enumerate(all_gfp_bcs):
        mapped = False
        for bc2 in good_gfp_bcs:
            h = hamming(bc1,bc2)
            if h <= sc_correction_treshold:
                mapped = True
                bc_map[bc1] = bc2
                break
        if not mapped:
            good_gfp_bcs.append(bc1)
    for bc in good_gfp_bcs: 
        bc_map[bc] = bc

    # Correct GBCs
    sc_df['GBC'] = sc_df['GBC'].map(bc_map)                               

    # Add stats
    grouped = sc_df.groupby(['CBC', 'GBC'])
    df_combos = (
        grouped.size().to_frame('read_counts').join(
        grouped['UMI'].nunique().to_frame('umi_counts')
        ).reset_index()
    )

    return df_combos, bulk, fig


##


def filter_and_pivot_LARRY2020(df_combos, read_treshold=30, umi_treshold=5):
    """
    Filter a CBC-UMI-GBC table and return it in large format.
    """
    # Supported and unsupported CBC-GBC combinations
    test = (df_combos['read_counts'] >= read_treshold) & (df_combos['umi_counts'] >= umi_treshold)
    df_combos['status'] = np.where(test, 'supported', 'unsupported')

    # Pivot
    M = (
        df_combos
        .query('status=="supported"')
        .pivot_table(index='CBC', columns='GBC', values='umi_counts')
    )
    M[M.isna()] = 0

    return M, df_combos


##


def get_clones(M):
    """
    Get the clone dataframe for each sample.
    """
    cells_with_unique_combos = M.apply(lambda x: frozenset(M.columns[x>0]), axis=1)
    cells_with_unique_combos = cells_with_unique_combos.map(lambda x: ';'.join(x))
    clones = (
        cells_with_unique_combos
        .reset_index(drop=True)
        .value_counts(normalize=False)
        .to_frame('n cells').reset_index()
        .rename(columns={'index':'GBC_set'})
        .assign(clone=lambda x: [ f'C{i}' for i in range(x.shape[0]) ])
        .set_index('clone')
    )
    clones['prevalence'] = clones['n cells'] / clones['n cells'].sum()

    return clones


##


def LARRY_2020_workflow(path_bulk, path_sample_map, path_sc, sample, 
    sc_correction_treshold=3, read_treshold=30, umi_treshold=5, coverage_treshold=10
    ):
    """
    Complete clone calling workflow.
    Read bulk and single-cell data, correct sc GBCs with the bulk reference, obtain a df with CBC-GBC combos.
    """

    T = Timer()

    T.start()
    f = open('clone_calling_summary.txt', 'w')
    f.write(f'# LARRY_2020 clone calling and cell assignment workflow: \n')
    f.write(f'\n')
    f.write(f'Input params:\n')
    f.write(f'  * sc_correction_treshold: {sc_correction_treshold}\n')
    f.write(f'  * read_treshold: {read_treshold}\n')
    f.write(f'  * umi_treshold: {umi_treshold}\n')
    f.write(f'  * coverage_treshold: {coverage_treshold}\n')
    f.write(f'\n')

    # Get CBC-GBC combos
    df_combos, df_bulk, fig = get_combos_LARRY2020(
        path_bulk, path_sample_map, sample, path_sc, 
        sc_correction_treshold=sc_correction_treshold, 
        coverage_treshold=coverage_treshold
    )
    fig.savefig('CBC_GBC_UMI_read_distribution.png')

    # Filter combos and pivot
    M, df_combos = filter_and_pivot_LARRY2020(df_combos, read_treshold=read_treshold, umi_treshold=umi_treshold)

    # General checks
    f.write('# General checks \n')
    f.write(f'- Unique, "good" GBCs, post correction (without bulk reference whitelist): {df_combos["GBC"].unique().size}\n')
    f.write(f'- n unsupported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="unsupported").sum()}\n')
    f.write(f'- n supported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="supported").sum()}\n')
    f.write(f'- Observed MOI (from supported CBC-GBC combos only): {(M>0).sum(axis=1).median():.2f} median, (+-{(M>0).sum(axis=1).std():.2f})\n')
    f.write(f'\n')

    # Relationship with bulk (GBC sequences) checks
    pseudobulk_sc = M.sum(axis=0) / M.sum(axis=0).sum()
    pseudobulk_sc = pseudobulk_sc.sort_values(ascending=False)
    common = set(pseudobulk_sc.index.to_list()) & set(df_bulk.index.to_list())
    pseudobulk_sc = pseudobulk_sc.loc[list(common)]
    bulk = df_bulk.loc[list(common)]['read_count'] / df_bulk.loc[list(common)]['read_count'].sum()
    corr = np.corrcoef(pseudobulk_sc, bulk)[0,1]
    f.write('# Individual "good" GBC sequences checks \n')
    f.write(f'- n "good" GBC sequences, bulk ({df_bulk.shape[0]}) vs sc ({df_combos["GBC"].unique().size})\n')
    f.write(f'- n common GBC: {len(common)}\n')
    f.write(f'- Fraction of bulk GBC found in sc: {df_bulk.index.isin(df_combos["GBC"].unique()).sum()/df_bulk.shape[0]:.2f}\n')
    f.write(f'- Common GBCs abundance (normalized nUMIs from pseudobulk scRNA-seq vs normalized read counts from bulk DNA-seq) correlation: {corr:.2f}\n')
    f.write(f'\n')

    # GBC sets checks
    sets = get_clones(M)
    GBC_set = list(chain.from_iterable(sets['GBC_set'].map(lambda x: x.split(';')).to_list()))
    redundancy = 1-np.unique(GBC_set).size/len(GBC_set)
    occurrences = pd.Series(GBC_set).value_counts().sort_values(ascending=False)
    f.write('# GBC sets (i.e. "good" GBCs combinations which are supported in a given cell subset) checks \n')
    f.write(f'- Unique GBCs sets: {sets.shape[0]}\n')
    f.write(f'- Unique GBCs in these sets: {np.unique(GBC_set).size}\n')
    f.write(f'- GBCs redundancy across sets: {redundancy:.2f}\n')
    f.write(f'- GBCs occurrences across sets: {occurrences.median():.2f} (+- {occurrences.std():.2f})\n')
    f.write('- Ratio between cell counts assigned to a single GBC (top 10 highly-recurrent GBCs, or "frequent flyers") \n')
    f.write('  and the total cell counts of other "ambiguous", possibly related GBC_sets (i.e., GBC sets containing the given GBC combined with other ones).\n')
    n = min([occurrences.size, sets.shape[0], 10])
    for i,x in enumerate(occurrences.index[:n]):
        clone = sets.iloc[i,1]
        total = sets.loc[sets['GBC_set'].str.contains(x)]['n cells'].sum()
        f.write(    f'  -> GBC {i+1} ({x}): occurences {occurrences[i]}; cell_count ratio: {clone/total:.2f}\n')
    f.write(f'\n')

    # Get clones (unique GBC only here) and cells tables
    unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
    filtered_M = M.loc[unique_cells]
    clones_df = get_clones(filtered_M)
    cells_df = (
        filtered_M
        .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
        .to_frame('GBC_set')
    )

    # Final clones checks
    common = set(clones_df.set_index('GBC_set').index.to_list()) & \
             set(df_bulk.index.to_list())
    corr = np.corrcoef(
        clones_df.set_index('GBC_set').loc[list(common), 'prevalence'], 
        bulk.loc[list(common)]
    )[0,1]
    f.write('# Final clones (i.e., distinct populations of uniquely barcoded cells only) checks \n')
    f.write(f'- n starting CBC (STARSolo): {df_combos["CBC"].unique().size}\n')
    f.write(f'- n uniquely barcoded cells: {cells_df.shape[0]}\n')
    f.write(f'- n clones: {clones_df.shape[0]}\n')
    f.write(f'- n clones>=10 cells: {clones_df["n cells"].loc[lambda x:x>=10].size}\n')
    f.write(f'- Common sc-bulk clones (N.B., one barcode only, to calculate prevalences correlation): {len(common)}\n')
    f.write(f'- sc vs bulk prevalence correlation: {corr:.2f}\n')
    f.write(f'\n')

    # Save
    df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')
    clones_df.to_csv('clones_summary_table.csv')
    cells_df.to_csv('cells_summary_table.csv')


    ##


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


    f.write(f'Finished workflow in time {T.stop()}.')
    f.close()


##
    