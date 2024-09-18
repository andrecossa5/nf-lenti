"""
Nanopore SCMseq exploration.
"""

import os
import numpy as np
import pandas as pd
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Explore
sample = 'sAML1'
path_counts = os.path.join(f'/Users/IEO5505/Desktop/example_mito/results/test_nanopore_{sample}', 'counts_table.csv') 
path_consensus = os.path.join(f'/Users/IEO5505/Desktop/example_mito/results/test_nanopore_{sample}', 'cons_stats.csv') 

allele_counts = pd.read_csv(path_counts)
consensus_df = pd.read_csv(path_consensus).pivot(index='cell', values='value', columns='metric')

consensus_df.median()
consensus_df['n_chrom'].describe()

df = allele_counts.query('MUT>0 or WT>0')
df['cell'].nunique()
df.query('gene=="FLT3"')

df.query('cell=="AAAGGATTCCCAGACA"')

# df['bin'] = np.where(df['MUT']>=.05*df['WT'],1,0).astype(int)
# df = df.pivot(index='cell', columns='gene', values='bin')
# df.apply(lambda x: np.sum(~x.isna()), axis=0)
# (df==1).sum(axis=0)
# (df==0).sum(axis=0)

MUT = df.pivot(index='cell', columns='gene', values='MUT').fillna(0)
WT = df.pivot(index='cell', columns='gene', values='WT').fillna(0)

coverage = MUT+WT

l = []
for n in range(1,7):
    l.append((coverage>0).sum(axis=1).loc[lambda x: x>=n].index.size)

fig, ax = plt.subplots(figsize=(5,5))
ax.plot(range(1,7), l, 'ko', markersize=10)
format_ax(title=sample, xlabel='# nuclear SNVs covered by SCM-seq', ylabel='# cells', ax=ax)
fig.tight_layout()
plt.show()











# MUT['RUNX1'].sum()
# WT['RUNX1'].sum()