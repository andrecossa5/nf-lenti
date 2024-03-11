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
from helpers import*


##





def read_data_GC(path_sc):    

    """
    Create a table of CBC-UMI-GBC combinations from single-cell data, after correcting 
    GBCs with a bulk reference.
    """
    
    
    # Reverse complement and value_counts
    sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
    sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
    d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))
    #return sc_df, bulk
    return sc_df



##



 
##

def cell_assignment_workflow_GC(
    path_sc,
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
    # Read and count
    sc_df= read_data_GC(path_sc)

    # Count
    COUNTS = {}
    COUNTS['raw'] = count_UMIs(sc_df)
    print('Correction')
    # Correction
    sc_map = map_GBCs(sc_df, sc_correction_treshold=sc_correction_treshold)
    sc_df['GBC_reference-free'] = sc_df['GBC'].map(sc_map)
    COUNTS['reference-free'] = count_UMIs(sc_df, gbc_col='GBC_reference-free')
    

    # % of raw counts lost with bulk
    assert COUNTS['raw']['count'].sum() == COUNTS['reference-free']['count'].sum()
    
    # % GBCs retained with sc and bulk corrections
    perc_gbc_retained_sc = COUNTS['reference-free']['GBC_reference-free'].value_counts().size / \
                           COUNTS['raw']['GBC'].value_counts().size
    
    
    
    ##
    
    print('Viz correction effect')
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

    # Save
    fig.tight_layout()
    fig.savefig('CBC_GBC_UMI_read_distribution.png', dpi=300)

    # Save counts as pickle
    with open('counts.pickle_GC', 'wb') as p:
        pickle.dump(COUNTS, p)

    
    ##

    
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