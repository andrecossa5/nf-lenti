"""
Nanopore SCMseq exploration.
"""

import os
import numpy as np
import pandas as pd


##


# Explore
path_counts = os.path.join('/Users/IEO5505/Desktop/example_mito/results/test_nanopore', 'counts_table.csv') 
path_consensus = os.path.join('/Users/IEO5505/Desktop/example_mito/results/test_nanopore', 'cons_stats.csv') 

allele_counts = pd.read_csv(path_counts)
consensus_df = pd.read_csv(path_consensus).pivot(index='cell', values='value', columns='metric')

df = allele_counts.query('MUT>=2 or WT>=2')
df['bin'] = np.where(df['MUT']>=2*df['WT'],1,0).astype(int)

df = df.pivot(index='cell', columns='gene', values='bin')
df.apply(lambda x: np.sum(~x.isna()), axis=0)
(df==1).sum(axis=0)
(df==0).sum(axis=0)

top_three = ((~df.isna()).sum(axis=1)>=3).loc[lambda x: x].index
df.loc[top_three]

