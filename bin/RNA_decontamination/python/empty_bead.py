"""
Utils for custom cell_assignement. Brings together elements from LARRY_2020, CR_2023 and
other (pre-)pubblication works.
g
"""

import os
import sys
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
import umap
import scipy.sparse as sp
from sklearn.cluster import KMeans
from scipy.sparse import csr_matrix

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")

from helpers import *


def viz_UMIs(counts, ax, log=True, by=None, nbins='sturges', colors = 'k' ):
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
        hist(counts, value_type, c= colors, ax=ax, n=nbins, a=.7)
    format_ax(
        xticks=np.logspace(0,4,5), ax=ax, 
        xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
        title='CBC-GBC-UMI combination'
    )
    if log:
        ax.set_yscale('log')

def viz_UMIs_difference(data1, data2, ax, log=True, nbins='sturges', colors = 'k' ):
    """
    Plot UMI n reads distribution.
    """
    value_type = 'log' if log else 'count'
    nbins = nbins if nbins != 'sturges' else sturges(counts[value_type])
    counts1, bins1, _ = plt.hist(data1, value_type, c= colors, n=nbins, a=.7, label='Histogram 1')
    counts2, bins2, _ = plt.hist(data2, value_type, c= colors, n=nbins, a=.7, label='Histogram 2')
    
    difference = counts1 - counts2

    plt.bar(bins1[:-1], difference, width=(bins1[1]-bins1[0]), align='edge', alpha=0.5, color='red', label='Difference')

    format_ax(
        xticks=np.logspace(0,4,5), ax=ax, 
        xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
        title='CBC-GBC-UMI combination'
    )
        # Aggiungi legenda e titolo
    plt.legend()
    plt.title('Difference between Histograms')
    plt.xlabel('Bins')
    plt.ylabel('Difference')
    if log:
        ax.set_yscale('log')




# Parametri scratch reference free
#path_sc = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/GBC_read_elements.tsv.gz'
path_sc = '/Users/ieo6943/Documents/Guido/export_cluster/GBC_read_elements.tsv.gz'
path_count='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/counts.pickle'

sample_params=None
correction_type='reference-free' 
sc_correction_treshold=3
ncores=8
filtering_method="fixed"
coverage_treshold=0 
umi_treshold=5
p_treshold=1
max_ratio_treshold=.8
normalized_abundance_treshold=.8
   
# Reverse complement and value_counts
#sc_df, bulk_df = read_data(path_bulk, path_sample_map, path_sc, sample=sample)
sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))
#
## Count
COUNTS_empty = {}
COUNTS_empty['raw'] = count_UMIs(sc_df)

# COUNT reference

sc_map = map_GBCs(sc_df, sc_correction_treshold=sc_correction_treshold)
#bulk_map = map_GBCs(sc_df, bulk_df, bulk_correction_treshold=bulk_correction_treshold, ncores=ncores)
sc_df['GBC_reference-free'] = sc_df['GBC'].map(sc_map)
#sc_df['GBC_reference'] = sc_df['GBC'].map(bulk_map)

COUNTS_empty['reference-free'] = count_UMIs(sc_df, gbc_col='GBC_reference-free')
COUNTS_empty['reference'] = count_UMIs(sc_df, gbc_col='GBC_reference')

path = '/Users/ieo6943/Documents/Guido/empty_bead/'
data_empty = COUNTS_empty
coverage_treshold = 0
#for coverage_treshold in list(range(0, 101, 20)):
print(coverage_treshold)
# Mark noisy UMIs, from chosen correction method
counts_empty = mark_UMIs(data_empty[correction_type], coverage_treshold=coverage_treshold)
# Get CBC-GBC combos with selected UMIs
df_combos_empty = get_combos(counts_empty, gbc_col=f'GBC_{correction_type}')
df = df_combos_empty
M_empty = df_combos_empty.pivot_table(index='CBC', columns='GBC', values='umi')
M_empty[M_empty.isna()] = 0
M_empty_l = pd.DataFrame(M_empty.idxmax(axis=1), columns=['GBC_set'])
counts_empty_new = counts_empty[counts_empty['GBC_reference-free'].isin(M_empty_l['GBC_set'])]
#M_empty = coo_matrix((df['umi'], (df['CBC'].astype('category').cat.codes, df['GBC'].astype('category').cat.codes)), dtype=custom_dtype)

plt.figure(figsize=(12, 6))
sns.kdeplot(data=df_combos_empty, x=df_combos_empty['max_ratio'], y=df_combos_empty['normalized_abundance'], label='Before decontamination')
plt.title(f'Density Plot Empty bead {coverage_treshold}')
plt.savefig(path + f'density_plot_{coverage_treshold}')







# DAL PICKLE


with open(path_count, 'rb') as p:
    data_pickle= pickle.load(p)
data = data_pickle
counts = mark_UMIs(data[correction_type], coverage_treshold=coverage_treshold)
df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
M =df_combos.pivot_table(index='CBC', columns='GBC', values='umi')
M[M.isna()] = 0
M_l = pd.DataFrame(M.idxmax(axis=1), columns=['GBC_set'])
counts_new = counts[counts['CBC'].isin(M_l.index) & counts['GBC_reference-free'].isin(M_l['GBC_set'])]



fig, ax = plt.subplots(figsize=(8, 6))

viz_UMIs(counts_new, ax, colors= 'g')
ax.set(title='good-bead:reference-free')
ax.text(.53, .85, '----------------------------------------------------------')
ax.text(.53, .83, f'Total reads: {counts_new["count"].sum():.2e}', transform=ax.transAxes)
ax.text(.53, .79, f'Total GBCs: {counts_new["GBC_reference-free"].unique().size:.2e}', transform=ax.transAxes)
ax.text(.53, .75, f'n reads: {counts_new["count"].median():.2f} (+-{counts_new["count"].std():.2f})', 
            transform=ax.transAxes)

viz_UMIs(counts_empty_new, ax, colors = 'r' )
ax.set(title='Empty bead reference-free')
ax.text(.53, .95, f'Total reads empty: {counts_empty_new["count"].sum():.2e}', transform=ax.transAxes)
ax.text(.53, .91, f'Total GBCs empty: {counts_empty_new["GBC_reference-free"].unique().size:.2e}', transform=ax.transAxes)
ax.text(.53, .87, f'n reads empty: {counts_empty_new["count"].median():.2f} (+-{counts_empty_new["count"].std():.2f})', 
        transform=ax.transAxes)



# Save
fig.tight_layout()
fig.savefig(path + f'Raw_CBC_GBC_UMI_read_distribution_empty_vs_not.png', dpi=300)


fig, ax = plt.subplots(figsize=(8, 6))

viz_UMIs_difference(counts_new,counts_empty_new, ax, colors= 'g')
# Save
fig.tight_layout()
fig.savefig(path + f'Raw_CBC_GBC_UMI_read_distribution_difference.png', dpi=300)
