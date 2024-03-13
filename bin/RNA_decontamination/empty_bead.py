"""
Visualization of the counts supporting CBC-GBC-UMI comparing with empty bead
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting_utils._plotting_base import *

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
from decontamination_utils import *
from helpers import *


##

# Path
#path_sc = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/GBC_read_elements.tsv.gz'
path_main = '/Users/ieo6943/Documents/data/8_clones/'
path_decontX = os.path.join(path_main,'count_matrix_&_DecontX')
path_count = os.path.join(path_main,'counts.pickle')
path_results = '/Users/ieo6943/Documents/results/8_clones/'
path_sc = '/Users/ieo6943/Documents/Guido/export_cluster/GBC_read_elements.tsv.gz'
path_bulk = '/Users/ieo6943/Documents/data/8_clones/'
path_sample_map = '/Users/ieo6943/Documents/data/8_clones/sample_map.csv'


##


#Args
sample_params=None
correction_type='reference' 
sc_correction_treshold=3
bulk_correction_treshold = 3
ncores=8
filtering_method="fixed"
coverage_treshold=0 
umi_treshold=5
p_treshold=1
max_ratio_treshold=.8
normalized_abundance_treshold=.8
sample = 'MDA_clones'
   

##


# Create count matrix for the empty beads
sc_df, bulk_df = read_data(path_bulk, path_sample_map, path_sc, sample=sample)
#sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))

## Count
COUNTS_empty = {}
COUNTS_empty['raw'] = count_UMIs(sc_df)

# COUNT reference

sc_map = map_GBCs(sc_df, sc_correction_treshold=sc_correction_treshold)
bulk_map = map_GBCs(sc_df, bulk_df, bulk_correction_treshold=bulk_correction_treshold, ncores=ncores)
sc_df['GBC_reference-free'] = sc_df['GBC'].map(sc_map)
sc_df['GBC_reference'] = sc_df['GBC'].map(bulk_map)
COUNTS_empty['reference-free'] = count_UMIs(sc_df, gbc_col='GBC_reference-free')
COUNTS_empty['reference-free'].to_csv(os.path.join(path_main,'count_matrix_&_DecontX', f'COUNTS_empty_ref-free.csv'))
COUNTS_empty['reference'] = count_UMIs(sc_df, gbc_col='GBC_reference')

# Mark noisy UMIs, from chosen correction method
counts_empty = mark_UMIs(COUNTS_empty[correction_type], coverage_treshold=coverage_treshold)
# Get CBC-GBC combos with selected UMIs
df_combos_empty = get_combos(counts_empty, gbc_col=f'GBC_{correction_type}')
M_empty = df_combos_empty.pivot_table(index='CBC', columns='GBC', values='umi')
M_empty[M_empty.isna()] = 0
#M_empty.to_csv(os.path.join(path_main,'count_matrix_&_DecontX', f'M_empty_{coverage_treshold}.csv'))
#M_empty = pd.read_csv(os.path.join(path_decontX, f'M_empty_{coverage_treshold}.csv'), index_col=0) 
M_empty_l = pd.DataFrame(M_empty.idxmax(axis=1), columns=['GBC_set'])
counts_empty_major = counts_empty[counts_empty['GBC_reference'].isin(M_empty_l['GBC_set'])]
counts_empty_other = counts_empty[~counts_empty['GBC_reference'].isin(M_empty_l['GBC_set'])]


##


# Goo bead from the PICKLE
with open(path_count, 'rb') as p:
    data_pickle= pickle.load(p)

data = data_pickle
counts = mark_UMIs(data[correction_type], coverage_treshold=coverage_treshold)
df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
M =df_combos.pivot_table(index='CBC', columns='GBC', values='umi')
M[M.isna()] = 0
M_l = pd.DataFrame(M.idxmax(axis=1), columns=['GBC_set'])
counts_major = counts[counts['CBC'].isin(M_l.index) & counts['GBC_reference'].isin(M_l['GBC_set'])]
counts_other = counts[~counts['GBC_reference'].isin(M_l['GBC_set'])]


##


fig, axs = plt.subplots(1, 3, figsize=(15, 5))

##

value_type = 'count'
nbins = sturges(counts[value_type])
sns.kdeplot(counts[value_type], color='r', ax=axs[0], label='good bead')
sns.kdeplot(counts_empty[value_type], color='g', ax=axs[0], label='empty bead')
format_ax(
    ax=axs[0], 
    xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
    title='CBC-GBC-UMI combination'
)
axs[0].set_xlim(right=counts[value_type].max())  
axs[0].set_ylim(bottom=0)  

##


nbins = sturges(counts_empty[value_type])
sns.kdeplot(counts_major[value_type], color='r', ax=axs[1], label='good bead')
sns.kdeplot(counts_empty_major[value_type], color='g', ax=axs[1], label='empty bead')
format_ax(
    ax=axs[1], 
    xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
    title='CBC-GBC-UMI combination'
)
axs[1].set_xlim(right=counts_empty[value_type].max())  
axs[1].set_ylim(bottom=0)  

##


nbins = sturges(counts_empty[value_type])
sns.kdeplot(counts_other[value_type], color='r', ax=axs[2], label='good bead')
sns.kdeplot(counts_empty_other[value_type], color='g', ax=axs[2], label='empty bead')
format_ax(
    ax=axs[2], 
    xlabel='n reads', ylabel='n CBC-GBC-UMI combination',
    title='CBC-GBC-UMI combination'
)
axs[2].set_xlim(right=counts_empty[value_type].max())  
axs[2].set_ylim(bottom=0)  

axs[0].set(title='All count')
axs[1].set(title='Counts supporting major GBC')
axs[2].set(title='Counts supporting other GBC')


axs[0].legend()
axs[1].legend()
axs[2].legend()


fig.tight_layout()
fig.savefig(os.path.join(path_results, 'Prova_3.png'), dpi=300)