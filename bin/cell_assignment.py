#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


##

# Paths
path_CBCs = sys.argv[1]
path_UMIs = sys.argv[2]
path_GBCs = sys.argv[3]
path_aligned_reads = sys.argv[4]

# Read components and aligned reads
CBCs = pd.read_csv(path_CBCs, index_col=0, sep='\t', header=None)
UMIs = pd.read_csv(path_UMIs, index_col=0, sep='\t', header=None)
GBCs = pd.read_csv(path_GBCs, index_col=0, sep='\t', header=None)
aligned_names = pd.read_csv(path_aligned_reads, index_col=0, sep='\t', header=None)

# Chek read names are identical
if not (CBCs.index == UMIs.index).all() and (UMIs.index == GBCs.index).all():
    raise ValueError('Read names of filtered CBCs, UMIs and GBCs tables are not identical. Something is wrong...')
# Merge in a single table
df = pd.concat([CBCs, UMIs, GBCs], axis=1)
df.columns = ['CBC', 'UMI', 'GBC']

# Filter for only reads aligned to the bulk reference, write all CBC-UMI-GBC table
df = df.loc[aligned_names.index,:]
df.to_csv('CBC_UMI_GBC_by_read.tsv.gz', sep='\t')

# df = pd.read_csv('/Users/IEO5505/Desktop/CBC_UMI_GBC_by_read.tsv.gz', sep='\t')

# Compute unique CBC-GBC combo table
df_read_counts = df.groupby(['CBC', 'GBC']).size().to_frame('read_counts')
df_umi_counts = df.reset_index().loc[:, ['CBC', 'GBC', 'UMI']].drop_duplicates().groupby(['CBC', 'GBC']).size().to_frame('umi_counts')
df_combos = df_read_counts.join(df_umi_counts).reset_index()
df_combos['coverage'] = df_combos['read_counts'] / df_combos['umi_counts']
df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')

# Cell classification: from Adamson et al. 2016

# Assigned and not-assigned CBC-GBC combos
test = (df_combos['read_counts'] > 30) & \
        (df_combos['umi_counts'] > 3) & \
        (df_combos['coverage'] > 5) # TOTUNE
df_combos['status'] = np.where(test, 1, 0)
# Cell MOI tally 
cell_MOI = df_combos.loc[:, ['CBC', 'GBC', 'status']
    ].query('status == 1').groupby('CBC').sum(
    ).rename(columns={'status': 'multiplicity'})
# Filter uniquely assigned cells and get their unique GBC
uniquely_assigned_cells = cell_MOI.query('multiplicity == 1').index

# Cell assignment plot
fig, ax = plt.subplots()
x = np.log10(df_combos['read_counts'])
y = np.log10(df_combos['umi_counts'])
ax.plot(x[df_combos['status'] == 1], y[df_combos['status'] == 1], '.', label='assigned', color='blue')
ax.plot(x[df_combos['status'] == 0], y[df_combos['status'] == 0], '.', label='not-assigned', color='grey')
ax.set(title='CBC-GBC combination status', xlabel='log10_read_counts', ylabel='log10_umi_counts')
ax.legend()
fig.savefig('CBC_GBC_combo_status.png')

# Compute summary tables: cells 
df_cells = df_combos.query('CBC in @uniquely_assigned_cells and status == 1').drop(
    columns='status').set_index('CBC')
# Assert we have taken the correct ones and save
assert df_cells.index.size == uniquely_assigned_cells.size
df_cells.to_csv('cells_summary_table.csv')

# Compute summary tables: clones
df_clones = df_cells.reset_index().groupby('GBC').size().to_frame('n_cells')
df_clones['prevalence'] = df_clones['n_cells'] / df_clones['n_cells'].sum()
df_clones = df_clones.sort_values('prevalence', ascending=False)
df_clones.to_csv('clones_summary_table.csv')