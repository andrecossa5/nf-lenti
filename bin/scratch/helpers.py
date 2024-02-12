"""
Utils for custom cell_assignement. Brings together elements from LARRY_2020, CR_2023 and
other (pre-)pubblication works.
g
"""

import os
import pickle
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances, r2_score, mean_squared_error
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from scipy.stats import poisson
from itertools import chain
from plotting_utils._plotting_base import *
from plotting_utils._utils import Timer






#def read_data(path_bulk, path_sample_map, path_sc, sample):
def read_data(path_bulk, path_sample_map, path_sc, sample):   
    print()

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

    sc = sc_df['GBC'].value_counts().sort_values(ascending=False)
    n_common = len(set(sc.index.to_list()) & set(bulk.index.to_list()))
    print(f'{n_common} common GBC between bulk ({bulk.index.size}) and sc ({sc_df["GBC"].unique().size})')
    n_top = len(set(sc.index[:10].to_list()) & set(bulk.index[:10].to_list()))
    print(f'{n_top} common GBC between bulk and sc top 10 clones.')
    
    return sc_df, bulk
    



##


def to_numeric(X):
    return np.select([X=='A', X=='T', X=='C', X=='G', X=='N'], [1,2,3,4,5], default=0)


##


def hamming(bc1, bc2): 
    return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])


##

 
def map_GBCs(sc_df, bulk=None, use_bulk=True, bulk_correction_treshold=3, sc_correction_treshold=3, ncores=8):
    """
    Correct sc GBC sequences with bulk-DNA reference. Outputs a sequence map: sc:bulk.
    """

    # Correct with bulk reference
    if use_bulk and bulk is not None:

        # Calculate hamming distance single-cell GBCs vs bulk reference GBCs
        print('Correct with bulk reference...')
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
        )
        d_corr = d_corr.query('hamming<=@bulk_correction_treshold')['correct_GBC'].to_dict()
    
    # Correct without bulk reference
    else:
        
        # Code readapted from https://github.com/AllonKleinLab/LARRY/blob/master/LARRY_for_10X.ipynb
        print('Use direct sc correction...')
        sc = sc_df['GBC'].value_counts()
        all_gbcs = sc.index.to_list()       # Ordered for n_reads. 
        good_gbcs = []
        d_corr = {}
        for bc1 in all_gbcs:
            mapped = False
            for bc2 in good_gbcs:
                if hamming(bc1,bc2) <= sc_correction_treshold:
                    mapped = True
                    d_corr[bc1] = bc2
                    break
            if not mapped:
                good_gbcs.append(bc1)
        for bc in good_gbcs: 
            d_corr[bc] = bc

    return d_corr


##


def count_UMIs(sc_df, gbc_col='GBC'):
    counts = (
        sc_df.groupby(['CBC', gbc_col, 'UMI'])
        .size()
        .reset_index(name='count')
    )
    counts['log'] = np.log(counts['count'])/np.log(10) # LARRY log transformation
    return counts


##


def mark_UMIs(counts, coverage_treshold=None, n_mixtures=2, nbins=50):
    """
    Mark noisy UMI read_counts.
    """

    if isinstance(coverage_treshold, int):
        coverage_treshold = coverage_treshold
        counts['status'] = np.where(counts['count']>=coverage_treshold, 'Retain', 'Filter out')
        return counts

    elif coverage_treshold == 'medstd':
        coverage_treshold = counts['count'].median() + counts['count'].std()
        counts['status'] = np.where(counts['count']>=coverage_treshold, 'Retain', 'Filter out')
        return counts

    elif coverage_treshold == 'GMM':
        x = counts['log']
        counts['bin'] = pd.cut(x, bins=nbins)
        x = counts.groupby('bin')['log'].median().values
        test = ~np.isnan(x)
        x = x[test].reshape(-1,1)
        gmm = GaussianMixture(n_components=n_mixtures, random_state=1234)  
        gmm.fit(x)
        top_component_idx = np.argsort(gmm.means_.flatten())[-1]
        counts['status'] = np.where(gmm.predict_proba(x)[:,top_component_idx]>.85, 'Retain', 'Filter out')
        return counts
    
    elif coverage_treshold == 'polyfit':

        # Prep data
        value_type = 'log'
        x = counts[value_type]
        counts['bin'] = pd.cut(x, bins=nbins)
        x = counts.groupby('bin')['log'].median().values
        test = ~np.isnan(x)
        x = x[test]
        y = counts.groupby('bin').size().values
        y = np.log10(y[test])

        # Fit polynomial of degree 3
        degree = 3
        a,b,c,d = np.polyfit(x, y, degree)
        fitted_curve = np.poly1d([a,b,c,d])
        y_fitted = fitted_curve(x)
        r_squared = r2_score(y, y_fitted)
        mse = mean_squared_error(y, y_fitted)

        # Obtain roots
        x0 = ( -2*b + np.sqrt((4*b**2)-(12*a*c)) ) / 6*a
        x1 = ( -2*b - np.sqrt((4*b**2)-(12*a*c)) ) / 6*a
        xmin = [x0, x1][np.argmin([x0,x1])]
        ymin = a*xmin**3 + b*xmin**2 + c*xmin + d

        # Define UMI status
        counts['status'] = np.where(counts['count']>=10**xmin, 'Retain', 'Filter out')
        d = {
            'counts':counts, 'x':x, 'y':y, 'y_fitted':y_fitted, 
            'r_squared':r_squared, 'mse':mse, 'xmin':xmin, 'ymin':ymin

        }

        return d


## 
    

def viz_polyfit(x, y, y_fitted, xmin, ymin, r_squared, mse, ax=None):
    """
    Visualize polynomial fit used to find UMI nreads cutoff.
    """
    ax.plot(x, y, 'ko', markersize=5)
    ax.plot(x, y_fitted)
    ax.plot(xmin, ymin, 'rx')
    ax.axvline(x=xmin, c='r')
    format_ax(ax, title='Find polynomial root', xlabel='n reads (log10)', ylabel='n CBC-GBC-UMI combinations (log10)')
    ax.text(.6, .9, f'R2: {r_squared:.2f}, MSE: {mse:.2f}', transform=ax.transAxes)
    ax.text(.6, .86, f'n_reads treshold: {round(10**xmin)}', transform=ax.transAxes)

    return ax


##


def sturges(x):   
    return round(1 + 2 * 3.322 * np.log(len(x))) 


##


def viz_UMIs(counts, ax, log=True, by=None, nbins='sturges'):
    """
    Plot UMI n reads distribution.
    """
    value_type = 'log' if log else 'count'
    nbins = nbins if nbins != 'sturges' else sturges(counts[value_type])

    if by is not None:
        colors = {'Filter out':'k', 'Retain':'r'}
        hist(counts, value_type, by=by, c=colors, ax=ax, n=nbins, a=.7)
        add_legend('UMI status', colors=colors, ax=ax, bbox_to_anchor=(1,1), loc='upper right',
                   label_size=10, ticks_size=9)
    else:
        hist(counts, value_type, c='k', ax=ax, n=nbins, a=.7)
    format_ax(
        xticks=np.logspace(0,4,5), ax=ax, 
        xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
        title='CBC-GBC-UMI combination'
    )
    if log:
        ax.set_yscale('log')

    return ax


##


def get_combos(counts, gbc_col='GBC_reference-free'):
    """
    Compute CBC-GBC stats from filtered UMIs.
    """
    filtered = counts.loc[lambda x: x['status']=='Retain']
    df_combos = (
        filtered
        .rename(columns={gbc_col:'GBC'})
        .groupby(['CBC', 'GBC'])['UMI'].nunique()
        .to_frame('umi').reset_index()
        .assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )
    )
    return df_combos


##


def filter_and_pivot(df_combos, umi_treshold=5, p_treshold=.01, max_ratio_treshold=.8, normalized_abundance_treshold=.8):
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
        (df_combos['p'] <= p_treshold) & \
        (df_combos['max_ratio'] >= max_ratio_treshold) & \
        (df_combos['normalized_abundance'] >= normalized_abundance_treshold)
    df_combos['status'] = np.where(test, 'supported', 'unsupported')

    # Filter and pivot
    M = (
        df_combos
        .query('status=="supported"')
        .pivot_table(index='CBC', columns='GBC', values='umi')
    )
    M[M.isna()] = 0

    # Print n 1-GBC cells
    print('Estimated MOI:')
    print((M>0).sum(axis=1).value_counts())

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



def cell_assignment_workflow(
    path_bulk, path_sample_map,path_sc , sample, 
    sample_params=None,
    correction_type='reference-free', sc_correction_treshold=3, 
    bulk_correction_treshold=3, ncores=8,  
    filtering_method="medstd", coverage_treshold=15, 
    umi_treshold=5, p_treshold=.5, max_ratio_treshold=.8, normalized_abundance_treshold=.8
    ):
    """
    Complete clone calling and cell assignment workflow.
    Read bulk and single-cell data, correct sc GBCs with the bulk reference, obtain a df with CBC-GBC combinations,
    filter them and assign cells to their clonal labels.
    """
    
    T = Timer()

    T.start()
    f = open('clone_calling_summary.txt', 'w')
    
    if sample_params is not None:
        correction_type = sample_params['correction_type']
        coverage_treshold = sample_params['coverage_treshold']
        umi_treshold = sample_params['umi_treshold']
        p_treshold = sample_params['p_treshold']
        max_ratio_treshold = sample_params['max_ratio_treshold']
        normalized_abundance_treshold = sample_params['normalized_abundance_treshold']
    elif filtering_method in ["medstd", "GMM"]:
        coverage_treshold = filtering_method
    else:
        print(f'Use fixed coverage treshold...{coverage_treshold}')

    f.write(f'# Custom clone calling and cell assignment workflow, sample {sample}: \n')
    f.write('\n')
    f.write('Input params:\n')
    f.write(f'  * bulk_correction_treshold: {bulk_correction_treshold}\n')
    f.write(f'  * coverage_treshold: {coverage_treshold} \n')
    f.write(f'  * umi_treshold: {umi_treshold}\n')
    f.write(f'  * p_treshold: {p_treshold}\n')
    f.write(f'  * max_ratio_treshold: {max_ratio_treshold}\n')
    f.write(f'  * normalized_abundance_treshold: {normalized_abundance_treshold}\n')
    f.write(f'\n')
    
    # Read and count
    #sc_df, bulk_df = read_data(path_bulk, path_sample_map, path_sc, sample=sample)
    sc_df= read_data(path_sc)

    # Count
    COUNTS = {}
    COUNTS['raw'] = count_UMIs(sc_df)
    print('Correction')
    # Correction
    sc_map = map_GBCs(sc_df, sc_correction_treshold=sc_correction_treshold)
    #bulk_map = map_GBCs(sc_df, bulk_df, bulk_correction_treshold=bulk_correction_treshold, ncores=ncores)
    sc_df['GBC_reference-free'] = sc_df['GBC'].map(sc_map)
    #sc_df['GBC_reference'] = sc_df['GBC'].map(bulk_map)
    COUNTS['reference-free'] = count_UMIs(sc_df, gbc_col='GBC_reference-free')
    #COUNTS['reference'] = count_UMIs(sc_df, gbc_col='GBC_reference')
    

    # % of raw counts lost with bulk
    assert COUNTS['raw']['count'].sum() == COUNTS['reference-free']['count'].sum()
    #perc_read_retained_bulk = COUNTS['reference']['count'].sum() / COUNTS['raw']['count'].sum()
    
    # % GBCs retained with sc and bulk corrections
    perc_gbc_retained_sc = COUNTS['reference-free']['GBC_reference-free'].value_counts().size / \
                           COUNTS['raw']['GBC'].value_counts().size
    #perc_gbc_retained_bulk = COUNTS['reference']['GBC_reference'].value_counts().size / \
                           #COUNTS['raw']['GBC'].value_counts().size
    
    
    ##
    
    print('Viz correction effect')
    # Viz correction effect
    fig, axs = plt.subplots(1,3,figsize=(15,5))
    
    viz_UMIs(COUNTS['raw'], axs[0])
    axs[0].set(title='Raw')
    axs[0].text(.53, .95, f'Total reads: {COUNTS["raw"]["count"].sum():.2e}', transform=axs[0].transAxes)
    axs[0].text(.53, .91, f'Total GBCs: {COUNTS["raw"]["GBC"].unique().size:.2e}', transform=axs[0].transAxes)
    axs[0].text(.53, .87, f'n reads: {COUNTS["raw"]["count"].median():.2f} (+-{COUNTS["raw"]["count"].std():.2f})', 
                transform=axs[0].transAxes)
    
    viz_UMIs(COUNTS['reference-free'], axs[1])
    axs[1].set(title='Reference-free correction')
    axs[1].text(.53, .95, f'% reads retained: 100%', transform=axs[1].transAxes)
    axs[1].text(.53, .91, f'% GBCs retained: {perc_gbc_retained_sc*100:.2f}%', transform=axs[1].transAxes)
    axs[1].text(.53, .87, f'n reads: {COUNTS["reference-free"]["count"].median():.2f} (+-{COUNTS["reference-free"]["count"].std():.2f})', 
                transform=axs[1].transAxes)
    """"
    viz_UMIs(COUNTS['reference'], axs[2])
    axs[2].set(title='Bulk-DNA reference correction (PT)')
    axs[2].text(.53, .95, f'Reads retained: {perc_read_retained_bulk*100:.2f}%', transform=axs[2].transAxes)
    axs[2].text(.53, .91, f'GBCs retained: {perc_gbc_retained_bulk*100:.2f}%', transform=axs[2].transAxes)
    axs[2].text(.53, .87, f'n reads: {COUNTS["reference"]["count"].median():.2f} (+-{COUNTS["reference"]["count"].std():.2f})', 
                transform=axs[2].transAxes)
    """
    # Save
    fig.tight_layout()
    fig.savefig('CBC_GBC_UMI_read_distribution.png', dpi=300)

    # Save counts as pickle
    with open('counts.pickle_GC', 'wb') as p:
        pickle.dump(COUNTS, p)

    
    ##

    print('Mark noisy')
    # Mark noisy UMIs, from chosen correction method
    counts = mark_UMIs(COUNTS[correction_type], coverage_treshold=coverage_treshold, n_mixtures=3)

    # Viz UMI distribution and selected UMIs
    fig, ax = plt.subplots(figsize=(5,5))
    viz_UMIs(counts, ax, by='status')
    ax.set(title=f'Filtered UMIs, after {correction_type} correction')
    fig.tight_layout()
    fig.savefig('selected_UMIs.png', dpi=300)

    # Get CBC-GBC combos with selected UMIs
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')

    # Filter CBC-GBC combos
    M, _ = filter_and_pivot(
        df_combos, 
        umi_treshold=umi_treshold, 
        p_treshold=p_treshold, 
        max_ratio_treshold=max_ratio_treshold,
        normalized_abundance_treshold=normalized_abundance_treshold
    )


    ## 


    # General checks
    f.write(f'# General checks \n')
    f.write(f'- Unique, "good" GBCs, post correction (with bulk reference whitelist): {df_combos["GBC"].unique().size}\n')
    f.write(f'- n unsupported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="unsupported").sum()}\n')
    f.write(f'- n supported CBC-GBC combos (N.B. only good GBCs): {(df_combos["status"]=="supported").sum()}\n')
    f.write(f'- Observed MOI (from supported CBC-GBC combos only): {(M>0).sum(axis=1).median():.2f} median, (+-{(M>0).sum(axis=1).std():.2f})\n')
    f.write(f'\n')
    """"
    # Relationship with bulk (GBC sequences) checks
    pseudobulk_sc = M.sum(axis=0) / M.sum(axis=0).sum()
    common = list(set(pseudobulk_sc.index) & set(bulk_df.index))
    if len(common)>0:
        pseudobulk_sc = pseudobulk_sc.loc[common]
        bulk = bulk_df.loc[common]['read_count'] / bulk_df.loc[common]['read_count'].sum()
        corr = np.corrcoef(pseudobulk_sc, bulk)[0,1]
        f.write(f'# Individual "good" GBC sequences checks\n')
        f.write(f'- n "good" GBC sequences, bulk ({bulk_df.shape[0]}) vs sc ({df_combos["GBC"].unique().size})\n')
        f.write(f'- Fraction of bulk GBCs found in sc: {bulk_df.index.isin(df_combos["GBC"].unique()).sum()/bulk_df.shape[0]:.2f}\n')
        f.write(f'- Common GBCs abundance (normalized nUMIs from pseudobulk scRNA-seq vs normalized read counts from bulk DNA-seq) correlation: {corr:.2f}\n')
        f.write(f'\n')
    """
    # GBC sets checks
    sets = get_clones(M)
    GBC_set = list(chain.from_iterable(sets['GBC_set'].map(lambda x: x.split(';')).to_list()))
    redundancy = 1-np.unique(GBC_set).size/len(GBC_set)
    occurrences = pd.Series(GBC_set).value_counts().sort_values(ascending=False)
    f.write(f'# GBC sets (i.e. "good" GBCs combinations which are supported in a given cell subset) checks \n')
    f.write(f'- Unique GBCs sets: {sets.shape[0]}\n')
    f.write(f'- Unique GBCs in these sets: {np.unique(GBC_set).size}\n')
    f.write(f'- GBCs redundancy across sets: {redundancy:.2f}\n')
    f.write(f'- GBCs occurrences across sets: {occurrences.median():.2f} (+- {occurrences.std():.2f})\n')
    f.write(f'\n')
    
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
    f.write(f'# Final clones (i.e., distinct populations of uniquely barcoded cells only) checks \n')
    f.write(f'- n starting CBC (STARSolo): {df_combos["CBC"].unique().size}\n')
    f.write(f'- n uniquely barcoded cells: {cells_df.shape[0]}\n')
    f.write(f'- n clones: {clones_df.shape[0]}\n')
    f.write(f'- n clones>=10 cells: {clones_df["n cells"].loc[lambda x:x>=10].size}\n')

    # Save
    df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')
    clones_df.to_csv('clones_summary_table.csv')
    cells_df.to_csv('cells_summary_table.csv')

    # Final CBC-GBC combo support plot

    # Viz p_poisson vs nUMIs
    fig, ax = plt.subplots(figsize=(5,5))
    scatter(df_combos, 'umi', 'p', by='max_ratio', marker='o', s=10, vmin=.2, vmax=.8, ax=ax, c='Spectral_r')
    format_ax(
        ax, title='p Poisson vs nUMIs, all CBC-GBC combinations', 
        xlabel='nUMIs', ylabel='p', reduce_spines=True
    )
    ax.axhline(y=p_treshold, color='k', linestyle='--')
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