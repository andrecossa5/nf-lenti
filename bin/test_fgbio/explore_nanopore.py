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
path_counts = os.path.join('/Users/IEO5505/Desktop/example_mito/test_nanopore', 'counts_table.csv') 
path_consensus = os.path.join('/Users/IEO5505/Desktop/example_mito/test_nanopore', 'cons_stats.csv') 

allele_counts = pd.read_csv(path_counts)
consensus_df = pd.read_csv(path_consensus).pivot(index='cell', values='value', columns='metric')

consensus_df.median()
df = allele_counts.query('MUT>0 or WT>0')
df['cell'].nunique()

df.query('gene=="FLT3"')

# df['bin'] = np.where(df['MUT']>=.05*df['WT'],1,0).astype(int)
# df = df.pivot(index='cell', columns='gene', values='bin')
# df.apply(lambda x: np.sum(~x.isna()), axis=0)
# (df==1).sum(axis=0)
# (df==0).sum(axis=0)

MUT = df.pivot(index='cell', columns='gene', values='MUT').fillna(0)
WT = df.pivot(index='cell', columns='gene', values='WT').fillna(0)

coverage = MUT+WT

(coverage>0).sum(axis=1).describe()

# top_three = ((~df.isna()).sum(axis=1)>=3).loc[lambda x: x].index
# df.loc[top_three]

MUT['RUNX1'].sum()
WT['RUNX1'].sum()

covered_transcripts = coverage.apply(lambda x: (x>0).sum(), axis=1)
SH = coverage.apply(lambda x: -np.sum( (x/x.sum()) * np.log2(x/x.sum()) ), axis=1)

#

sns.kdeplot(np.log10(coverage.sum(axis=1)))
sns.kdeplot(SH)
sns.histplot(covered_transcripts)
plt.show()

distribute_cov_cells = SH.loc[lambda x: x>0.5].index
np.sum(coverage.loc[distribute_cov_cells]>0, axis=1).median()

#

np.sum(MUT["RUNX1"]>0)
np.sum(WT["RUNX1"]>0)

##

coverage.loc[WT["RUNX1"].sort_values(ascending=False).index[:1]]