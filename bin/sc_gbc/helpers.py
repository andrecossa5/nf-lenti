"""
Utils for custom cell_assignement. Brings together elements from LARRY_2020, CR_2023 and
other (pre-)pubblication works.
"""

import os
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
from scipy.stats import poisson
from itertools import chain
from plotting_utils._plotting_base import *
from plotting_utils._utils import Timer


##


def to_numeric(X):
    return np.select([X=='A', X=='T', X=='C', X=='G'], [1,2,3,4], default=0)


##


def get_CBC_GBC_combos(path_bulk, path_sample_map, path_sc, sample, ncores=8, 
                        bulk_correction_treshold=1, coverage_treshold="auto"):
    """
    Create a table of CBC-UMI-GBC combinations from single-cell data, after correcting 
    GBCs with a bulk reference.
    """

    # Get the right reference sequences from bulk_GBC_reference 
    if path_bulk is not None:
        bulk = pd.read_csv(
            os.path.join(path_bulk, 'summary', 'bulk_GBC_reference.csv'),
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

    # Reverse complement and value_counts
    sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
    sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
    d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))

    # Count CBC-GBC-UMI combinations
    counts = sc_df.groupby(['CBC', 'GBC', 'UMI']).size().reset_index(name='count')
    medstd = counts['count'].median() + counts['count'].std()
    coverage_treshold = coverage_treshold if coverage_treshold != "auto" else medstd

    # Viz distribution n reads
    fig, axs = plt.subplots(1,2,figsize=(10,5))

    counts['log'] = np.log(counts['count'])/np.log(10)
    hist(counts, 'log', c='k', ax=axs[0], n=50, a=.7)
    format_ax(
        xticks=np.logspace(0,4,5), ax=axs[0], xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
        title='CBC-GBC-UMI combination \n n reads distribution before GBC correction'
    )
    axs[0].set_yscale('log')
    axs[0].axvline(np.log(coverage_treshold)/np.log(10), c='r', linewidth=3, linestyle='--')
 

    ##

    
    # Map sc GBCs to bulk GBCs

    # Calculate hamming distance single-cell GBCs vs bulk reference GBCs
    sc = sc_df['GBC'].value_counts()
    sc_numeric = to_numeric(np.vstack(sc.index.map(lambda x: np.array(list(x)))))
    bulk_numeric = to_numeric(np.vstack(bulk.index.map(lambda x: np.array(list(x)))))
    D = pairwise_distances(
        sc_numeric, bulk_numeric, metric='hamming', n_jobs=int(ncores)
    ) * sc_numeric.shape[1]

    # Build a correction dict for sc GBCs at hamming distance <= bulk_sc_treshold from a bulk one.
    d_corr = (
        sc.to_frame('read_count')
        .assign(
            correct_GBC=[ bulk.index[i] for i in D.argmin(axis=1) ],
            hamming=D.min(axis=1),
        )
        .query('hamming<=@bulk_correction_treshold')
        ['correct_GBC'].to_dict()
    )
    # Correct GBC sequences and remove not found ones
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: d_corr[x] if x in d_corr else 'not_found')
    sc_df = sc_df.query('GBC!="not_found"')

    # Count again
    del counts
    counts = sc_df.groupby(['CBC', 'GBC', 'UMI']).size().reset_index(name='count')
    medstd = counts['count'].median() + counts['count'].std()
    coverage_treshold = coverage_treshold if coverage_treshold != "auto" else medstd

    # Viz distribution n reads, after correction
    counts['log'] = np.log(counts['count'])/np.log(10)
    hist(counts, 'log', c='k', ax=axs[1], n=50, a=.7)
    format_ax(
        xticks=np.logspace(0,4,5), ax=axs[1], xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
        title='CBC-GBC-UMI combination \n n reads distribution after GBC correction'
    )
    axs[1].set_yscale('log')
    axs[1].axvline(np.log(coverage_treshold)/np.log(10), c='r', linewidth=3, linestyle='--')
    fig.tight_layout()

    # Compute CBC-GBC combos with filtered UMIs, and related stats
    filtered = counts.query('count>=@coverage_treshold')
    df_combos = (
        pd.merge(
            sc_df, filtered[['CBC', 'GBC', 'UMI']],
            on=['CBC', 'GBC', 'UMI'], how='inner'
        )
        .groupby(['CBC', 'GBC'])['UMI']
        .nunique().to_frame('umi')
        .reset_index()
        .assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )
    )

    return df_combos, bulk, fig


##


def filter_and_pivot(df_combos, umi_treshold=5, p_treshold=.001, ratio_to_most_abundant_treshold=.3):
    """
    Filter a CBC-GBC-UMI table and return it in large format.
    """

    # p of observing higher counts than the observed ones, assuming UMI counts are Poisson distributed
    avg_gbc = df_combos.groupby('GBC')['umi'].median()
    p_poisson = []
    for i in range(df_combos.shape[0]):
        x = df_combos.iloc[i,:]
        gbc = x['GBC']
        obs = x['umi']
        exp = avg_gbc.loc[gbc]
        p = 1-poisson.cdf(obs-1, mu=exp)
        p_poisson.append(p)
    df_combos['p'] = p_poisson
    df_combos['logp'] = -np.log10(p_poisson)
    df_combos['logumi'] = np.log10(df_combos['umi'])

    ##

    # Supported and unsupported CBC-GBC combinations
    test = (df_combos['umi'] >= umi_treshold) & \
        (df_combos['p'] >= p_treshold) & \
        (df_combos['max_ratio'] >= ratio_to_most_abundant_treshold)
    df_combos['status'] = np.where(test, 'supported', 'unsupported')

    # Filter and pivot
    M = (
        df_combos
        .query('status=="supported"')
        .pivot_table(index='CBC', columns='GBC', values='umi')
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


def custom_workflow(
    path_bulk, path_sample_map, path_sc, sample, ncores=8, bulk_correction_treshold=1, 
    coverage_treshold=None, umi_treshold=5, p_treshold=.001, ratio_to_most_abundant_treshold=.3
    ):
    """
    Complete clone calling workflow.
    Read bulk and single-cell data, correct sc GBCs with the bulk reference, obtain a df with CBC-GBC combinations,
    filter them and assign cells to their clonal labels.
    """

    T = Timer()

    T.start()
    f = open('clone_calling_summary.txt', 'w')
    f.write(f'# Custom clone calling and cell assignment workflow: \n')
    f.write(f'\n')
    f.write(f'Input params:\n')
    f.write(f'  * bulk_correction_treshold: {bulk_correction_treshold}\n')
    f.write(f'  * coverage_treshold: {coverage_treshold}\n')
    f.write(f'  * umi_treshold: {umi_treshold}\n')
    f.write(f'  * p_treshold: {p_treshold}\n')
    f.write(f'  * ratio_to_most_abundant_treshold: {ratio_to_most_abundant_treshold}\n')
    f.write(f'\n')

    df_combos, df_bulk, fig = get_CBC_GBC_combos(
        path_bulk, 
        path_sample_map, 
        path_sc, 
        sample, 
        ncores=ncores, 
        bulk_correction_treshold=bulk_correction_treshold,
        coverage_treshold=coverage_treshold
    )
    fig.savefig('CBC_GBC_UMI_read_distribution.png')
        
    # Filter CBC_GBC combos and pivot 
    M, df_combos = filter_and_pivot(
        df_combos,
        umi_treshold=umi_treshold,
        p_treshold=p_treshold,
        ratio_to_most_abundant_treshold=ratio_to_most_abundant_treshold
    )

    # General checks
    f.write('# General checks \n')
    f.write(f'- Unique, "good" GBCs, post correction (with bulk reference whitelist): {df_combos["GBC"].unique().size}\n')
    f.write(f'- n unsupported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="unsupported").sum()}\n')
    f.write(f'- n supported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="supported").sum()}\n')
    f.write(f'- Observed MOI (from supported CBC-GBC combos only): {(M>0).sum(axis=1).median():.2f} median, (+-{(M>0).sum(axis=1).std():.2f})\n')
    f.write(f'\n')

    # Relationship with bulk (GBC sequences) checks
    pseudobulk_sc = M.sum(axis=0) / M.sum(axis=0).sum()
    pseudobulk_sc = pseudobulk_sc.sort_values(ascending=False)
    bulk = df_bulk.loc[pseudobulk_sc.index]['read_count'] / df_bulk.loc[pseudobulk_sc.index]['read_count'].sum()
    corr = np.corrcoef(pseudobulk_sc, bulk)[0,1]
    f.write('# Individual "good" GBC sequences checks\n')
    f.write(f'- n "good" GBC sequences, bulk ({df_bulk.shape[0]}) vs sc ({df_combos["GBC"].unique().size})\n')
    f.write(f'- Fraction of bulk GBCs found in sc: {df_bulk.index.isin(df_combos["GBC"].unique()).sum()/df_bulk.shape[0]:.2f}\n')
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

    ## 

    # Get clones (from cell with 1 single GBC only here) and cells tables

    # Get 1-GBC CBCs
    unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
    filtered_M = M.loc[unique_cells]
    clones_df = get_clones(filtered_M)
    cells_df = (
        filtered_M
        .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
        .to_frame('GBC_set')
    )

    # Final clones checks
    corr = np.corrcoef(clones_df.set_index('GBC_set')['prevalence'], bulk.loc[clones_df['GBC_set']])[0,1]
    f.write('# Final clones (i.e., distinct populations of uniquely barcoded cells only) checks \n')
    f.write(f'- n starting CBC (STARSolo): {df_combos["CBC"].unique().size}\n')
    f.write(f'- n uniquely barcoded cells: {cells_df.shape[0]}\n')
    f.write(f'- n clones: {clones_df.shape[0]}\n')
    f.write(f'- n clones>=10 cells: {clones_df["n cells"].loc[lambda x:x>=10].size}\n')
    f.write(f'- sc vs bulk prevalence correlation: {corr:.2f}\n')
    f.write(f'\n')

    # Save
    df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')
    clones_df.to_csv('clones_summary_table.csv')
    cells_df.to_csv('cells_summary_table.csv')

    # Final CBC-GBC combo support plot

    # Viz p_poisson vs nUMIs
    fig, ax = plt.subplots(figsize=(5,5))
    scatter(df_combos, 'umi', 'logp', by='max_ratio', marker='o', s=10, vmin=.2, vmax=.8, ax=ax, c='Spectral_r')
    format_ax(
        ax, title='p Poisson vs nUMIs, all CBC-GBC combinations', 
        xlabel='nUMIs', ylabel='-log10(p_poisson)', reduce_spines=True
    )
    ax.axhline(y=-np.log10(p_treshold), color='k', linestyle='--')
    ax.text(.05, .9, f'Total CBC-GBC combo: {df_combos.shape[0]}', transform=ax.transAxes)
    n_filtered = df_combos.query('status=="supported"').shape[0]
    ax.text(.05, .86, 
        f'n CBC-GBC combo retained: {n_filtered} ({n_filtered/df_combos.shape[0]*100:.2f}%)',
        transform=ax.transAxes
    )
    fig.tight_layout()
    fig.savefig('CBC_GBC_combo_status.png', dpi=300)


    ##


    f.write(f'- Execution time: {T.stop()}\n')
    f.close()


##