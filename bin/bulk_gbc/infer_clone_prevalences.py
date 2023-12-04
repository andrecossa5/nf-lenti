#!/usr/bin/python

# Infer clone prevalences script

########################################################################

# Parsing CLI args 

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='infer_clone_prevalences',
    description=
    """
    Script for final inference of clone prevalence from each sample bulk GBC counts.
    """
)

# Add arguments

# Input
my_parser.add_argument(
    '-i', 
    '--input', 
    type=str,
    default=None,
    help='The path to GBC_counts.csv file. Default: . .'
)

# Output
my_parser.add_argument(
    '-o', 
    '--output', 
    type=str,
    default=os.getcwd(),
    help='Output folder path. Default: . .'
)

# Path spikeins table
my_parser.add_argument( 
    '--path_spikeins_table', 
    type=str,
    default=None,
    help='''
    Path to spikeins table. Default: None. If different from None, the script will try to use spikeins for 
    clonal inference.
    '''
)

# Path spikeins table
my_parser.add_argument( 
    '--n_reads', 
    type=int,
    default=1000,
    help='Min n_reads to retain a GBC.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_i = args.input
path_o = args.output
path_spikeins_table = args.path_spikeins_table
n_reads = args.n_reads


##


########################################################################


# Import code
import numpy as np
import pandas as pd
import seaborn as sns
import dask.dataframe as dd
from scipy.interpolate import interp1d
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt


##


# Define helper functions
def rev_complement(x):
    d = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
    x = list(x)
    rev_x = []
    for i in range(len(x)):
        rev_x.append(d[x[i]])
    return ''.join(rev_x)[::-1]


##


def calc_FC(s):
     
    FC = []
    for i in range(s.size-1):
        fc = s[i] / s[i+1]
        FC.append(fc)
    FC.append(1)

    return np.array(FC)


##


def fit_trend(df):
    """
    Interpolation read_count, n_cells (log scale) and visualization.
    """
    x = df['log10_read_count']
    y = df['log10_n_cells']
    r_counts_cells = np.corrcoef(x, y)[0,1]
    r_FCs = np.corrcoef(df['FC_from_read_count'], df['FC_from_n_cells'])[0,1]
    f = interp1d(x, y, assume_sorted=False, fill_value='extrapolate')

    fig, ax = plt.subplots(figsize=(5.5,5))
    ax.plot(x, y, 'ko')
    sns.regplot(data=df, x='log10_read_count', y='log10_n_cells', scatter=False, ax=ax)
    title = f'Read_counts vs n_cells: r={r_counts_cells:.2f} \n FC from read_counts vs FC from n_cells: r={r_FCs:.2f}'
    ax.set(xlabel='log10 read_count', ylabel='log10 n_cells', title=title)
    ax.spines[['right', 'top']].set_visible(False)
    fig.tight_layout()
    
    return f, fig


##


def filter_spikeins_n_reads(df_counts, df_spike, n_reads=1):
    """
    Filter all spikeins from GBC counts and retain GBCs with at least n_reads.
    Than, recalculate log10_read_count and obs_frequency.
    """
    is_not_spike = ~df_counts.index.isin(df_spike.index)
    df = df_counts.loc[is_not_spike].copy()

    if n_reads is not None:
        has_enough_reads = df_counts['read_count']>=n_reads
        df = df_counts[has_enough_reads].copy()
    else:
        raise ValueError('n_reads must be an integer >=1')
    
    df['log10_read_count'] = np.log10(df['read_count'])
    df['obs_frequency'] = df['read_count'] / df['read_count'].sum()

    return df


##


def plot_distributions(df_counts, df_spike, n_reads=1):
    """
    Plot distribution of log10 read counts and observed frequencies of GBCs, 
    after spikeins removal. GBCs can be filtered further to have at least n_reads reads.
    """

    # Filter all GBCs counts
    df = filter_spikeins_n_reads(df_counts, df_spike, n_reads=n_reads)

    # Read counts / GBC frequency distributions
    fig, axs = plt.subplots(1,2,figsize=(9,5))
    
    s = df['status'].value_counts().reset_index()
    s.columns = ['status', 'counts']
    counts_txt = s['status'].astype('str') + ': ' + s['counts'].astype('str')
    sns.kdeplot(data=df, x='log10_read_count', fill=True, hue='status', ax=axs[0])
    median_ = df["log10_read_count"].median()
    std_ = df["log10_read_count"].std()
    t = f'log10 read counts: {median_:.2f} (+-{std_:.2f})'
    axs[0].set(title=t)
    axs[0].legend(labels=counts_txt.to_list(), title='status', frameon=False)
    axs[0].spines[['right', 'top']].set_visible(False)
    
    sns.kdeplot(data=df, x='obs_frequency', fill=True, hue='status', ax=axs[1])
    median_ = df["obs_frequency"].median()
    std_ = df["obs_frequency"].std()
    t = f'Observed frequency: {median_:.2f} (+-{std_:.2f})'
    axs[1].set(title=t)
    axs[1].spines[['right', 'top']].set_visible(False)
    axs[1].legend_ = None
    
    fig.suptitle(f'Read counts and frequency distributions. All GBCs > {n_reads} read (n={df.shape[0]})')
    fig.tight_layout()

    return fig
    

##


def filter_GBCs(df_counts, df_spike=None, n_reads=1, f=None, use_spike=True):
    """
    Function to filter a read count table with and without spikeins info.
    """
        
    # Calculate log10 read_counts
    
    # W/i
    if use_spike:

        # Filter to a minimum amount of reads
        df = filter_spikeins_n_reads(df_counts, df_spike, n_reads=n_reads)

        # Use interpolation to infer cellular_prevalence
        df['log10_n_cells'] = df['log10_read_count'].map(lambda x: f(x))
        df['n_cells'] = np.round(10 ** df['log10_n_cells'])
        df['cellular_prevalence'] = df['n_cells'] / df['n_cells'].sum()

        print(df.head())
        pho = np.corrcoef(df["obs_frequency"], df["cellular_prevalence"])[0,1]
        print(f'Pearson\'s r: {pho}')

        df = df.loc[:, ~df.columns.str.contains('log10|cells|obs')].copy()

    else:
        
        # Filter out only spikeins
        df = filter_spikeins_n_reads(df_counts, df_spike, n_reads=n_reads/10)

        # Use a GMM to assign GBCs to 2 mixture components and retain the ones from the top one
        X = df['read_count'].values.reshape(-1,1)
        gmm = GaussianMixture(n_components=2, random_state=1234)  
        gmm.fit(X)
        top_component_idx = np.argsort(gmm.means_.flatten())[-1]
        test = gmm.predict_proba(X)[:,top_component_idx]>.85    # Retain only GBCs assigned to the top component
        df = df.loc[test]                                       # Filter
        df['cellular_prevalence'] = df['read_count'] / df['read_count'].sum()

        print(df.head())
        pho = np.corrcoef(df["obs_frequency"], df["cellular_prevalence"])[0,1]
        print(f'Pearson\'s r: {pho}')

        df = df.loc[:, ~df.columns.str.contains('log_10|cells|obs')].copy()

    return df


##


def plot_prevalences(df):
    """
    Plot prevalences.
    """

    # Filter for common
    df_common = df.loc[(df['found_wi'])&(df['found_wo'])]
    x = df_common['cellular_prevalence_wo'].values
    y = df_common['cellular_prevalence_wi'].values
    r = np.corrcoef(x, y)[0,1]

    fig, ax = plt.subplots(figsize=(5.5,5))
    ax.plot(x, y, 'ko')
    ax.set(
        xlabel='Cellular prevalence w/o spikeins', 
        ylabel='Cellular prevalence w/i spikeins', 
        title=f'Cellular prevalences (w/i, w/o spikeins, r={r:.2f})'
    )
    ax.text(.1, .9, f'n clones w/i: {df["found_wi"].sum()}', transform=ax.transAxes)
    ax.text(.1, .85, f'n clones w/o: {df["found_wo"].sum()}', transform=ax.transAxes)
    ax.text(.1, .8, f'n common: {df_common.shape[0]}', transform=ax.transAxes)
    ax.spines[['right', 'top']].set_visible(False)
    fig.tight_layout()
    
    return fig


##


########################################################################


def main():

    # GBCs counts
    df_counts = dd.read_csv(path_i)
    df_counts.columns = ['GBC', 'read_count', 'status']

    # Filter GBCs according to their correction status
    df_counts = df_counts.query('status != "degenerated_to_remove"').compute()
    df_counts = df_counts.set_index('GBC')
    df_counts['status'] = pd.Categorical(
        df_counts['status'], 
        categories=['corrected', 'not_whitelisted', 'degenerated_not_to_remove']
    )

    # With spikeins
    if path_spikeins_table is not None:
        if os.path.exists(path_spikeins_table):
            print(f'Spikeins detected... Use this info for clonal inference.')
        else:
            raise ValueError(f'Provided path {path_spikeins_table} does not exists!')

    # Read spikeins table and check the detection of the 'centrl ones'
    df_spike = pd.read_csv(path_spikeins_table, index_col=0)
    df_spike.index = df_spike.index.map(lambda x: rev_complement(x))             
    # 2 orders dynamic range only, 3 points. 2 needed, bare minimum
    central_spikeins = df_spike.query('n_cells>=1000 and n_cells <=100000').index
    assert np.sum(df_counts.index.isin(central_spikeins)) >= 2

    # Fit spikeins trend
    df_spike = df_spike.join(df_counts)
    df_spike['log10_n_cells'] = np.log10(df_spike['n_cells'])
    df_spike['log10_read_count'] = np.log10(df_spike['read_count'])
    df_spike['FC_from_n_cells'] = calc_FC(df_spike['log10_n_cells'])
    df_spike['FC_from_read_count'] = calc_FC(df_spike['log10_read_count'])
    f, fig = fit_trend(df_spike.loc[central_spikeins])  # Trend only to central spikeins

    # Save spikeins table and fit
    (
        df_spike
        .assign(used_for_fit=lambda x: x.index.isin(central_spikeins))
        .to_csv(os.path.join(path_o, 'df_spikeins.csv'))
    )
    fig.savefig(os.path.join(path_o, 'spikeins_fit.png'), dpi=300)

    ##

    # Plot GBCs distributions, after spikeins (all of them) removal, and retaining GBCs with>n_reads
    fig = plot_distributions(df_counts, df_spike)
    fig.savefig(os.path.join(path_o, 'GBC_1read_distributions.png'), dpi=300)
    fig = plot_distributions(df_counts, df_spike, n_reads=n_reads)
    fig.savefig(os.path.join(path_o, f'GBC_morereads_distributions.png'), dpi=300)

    # W/i
    df_with = filter_GBCs(df_counts, df_spike=df_spike, n_reads=n_reads, f=f, use_spike=True)

    # W/o
    df_without = filter_GBCs(df_counts, df_spike=df_spike, n_reads=n_reads, use_spike=False)

    # Merge info
    df = (
        df_counts[['read_count', 'status']]
        .join(
            df_without[['cellular_prevalence']]
            .join(df_with[['cellular_prevalence']], lsuffix='_wo', rsuffix='_wi', how='outer')
            .assign(
                found_wi=lambda x: ~x['cellular_prevalence_wi'].isna(),
                found_wo=lambda x: ~x['cellular_prevalence_wo'].isna()
            ),
            how='right'
        )
        .sort_values('read_count', ascending=False)
    )
    df.to_csv(os.path.join(path_o, 'clonal_prevalences.csv'))

    # Last vizualization
    fig = plot_prevalences(df)
    fig.savefig(os.path.join(path_o, 'final_prevalences.png'), dpi=300)

##


########################################################################

# Run
if __name__ == '__main__':
    main()


