"""
New filtering procedure. With parallel pandas.
"""

import sys
import pickle
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import pairwise_distances
from plotting_utils._utils import *
from plotting_utils._plotting_base import *
sys.path.append("/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/RNA_decontamination")
sys.path.append("/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/sc_gbc")
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from decontamination_utils import *
from helpers import *


##


# Read counts
path_counts = '/Users/IEO5505/Desktop/mito_bench/data/AML_clones/counts.pickle'
with open(path_counts, 'rb') as p:
    counts = pickle.load(p)['raw']


# Example
counts.groupby(['CBC', 'UMI']).size()


# counts.query('GBC=="TATGGCGGTTGTGCTGTC"').groupby('CBC')['UMI'].nunique().describe()#apply(lambda x: x)#.sum().size().describe()
# counts.query('GBC=="CGGCGCTTCTGTGCCCGG"').groupby('CBC')['UMI'].nunique().describe()
# df = counts.query('CBC=="AAACCCAAGCCTCCAG" and UMI=="AACAGAGACTTC"').sort_values('count', ascending=False)


##


# Loop for each CBC/UMI combination
def process_consensus_UMI(df, correction_threshold=3):     # SPOSTARE in helpers, o simili, quando e' finita!
    """
    Process a CBC-UMI GBCs, and their read counts.
    """
    df = df.sort_values('count', ascending=False).set_index('GBC')
    GBCs = df.index

    if GBCs.size == 1:
        l = [ GBCs[0], df['count'].sum(), 1, np.nan ] 
    else:
        numeric_GBCs = np.array([ [ ord(char) for char in s ] for s in GBCs ])
        D = 18 * pdist(numeric_GBCs, metric='hamming')
        # D = 18 * pairwise_distances(numeric_GBCs, metric='hamming')                     # pdist
        if (D<=correction_threshold).mean():
            # Simplest alternative: majority barcode as consensus
            top_GBC = GBCs[0]                                   
            l = [ top_GBC, df['count'].sum(), df.loc[top_GBC, 'count']/df['count'].sum(), D.mean() ]
            # More advanced (but more computationally intensive...)
            """
            AAAC 40
            ACGC 35
            ACGT 30
            ACGA 25
            top_GBC = ACGC
            """
        else:
            l =  [np.nan] * 3 + [ D.mean() ]

    res = pd.DataFrame(l).T
    res.columns = ['GBC', 'total_reads', 'max_ratio', 'mean_hamming']

    return res


##


# # Process each cell UMI into a consensus GBC sequence.
# def f(ncores=1, n=10000, correction_threshold=5):
#     """
#     This is piece of code to optimize. Now embedded in a function f() just for time
#     monitoring convenience.
#     """
#     if ncores>1:
#         pandarallel.initialize(nb_workers=ncores)
#         processed = (
#             counts.sample(n)
#             .groupby(['CBC', 'UMI'])
#             .parallel_apply(lambda x: process_consensus_UMI(x, correction_threshold=correction_threshold))
#             .droplevel(2).reset_index()
#             .dropna()  
#         )
#     else:
#         processed = (
#             counts.sample(n)
#             .groupby(['CBC', 'UMI'])
#             .apply(lambda x: process_consensus_UMI(x, correction_threshold=correction_threshold))
#             .droplevel(2).reset_index()
#             .dropna()  
#         )
#     
#     print(f'Retained cells: {processed["CBC"].nunique()/counts.head(n)["CBC"].nunique()*100:.2f}%')     
#     print(f'Retained triplets: {processed.shape[0]/n*100:.2f}%') 
#     
#     return processed
# 
# 
# # Run code and monitor time
# run_command(f, n=10000, ncores=8, correction_threshold=7, verbose=True)


##


# Process
t = Timer()
t.start()
correction_threshold = 6
pandarallel.initialize(nb_workers=8)
processed = (
    counts
    .groupby(['CBC', 'UMI'])
    .parallel_apply(lambda x: process_consensus_UMI(x, correction_threshold=correction_threshold))
    .droplevel(2).reset_index()
    .dropna()  
)
print(f'Elapsed {t.stop()}')
print(f'Retained cells: {processed["CBC"].nunique()/counts["CBC"].nunique()*100:.2f}%')     
print(f'Retained triplets: {processed.shape[0]/counts.shape[0]*100:.2f}%') 
print(f'Retained GBCs: {processed["GBC"].nunique()/counts["GBC"].nunique()*100:.2f}%')   
processed["GBC"].nunique()
processed["GBC"].value_counts().head(10)#.describe()



##



# Before-after checks: CBCs and CBC-UMI-GBC retained
starting_cells = counts["CBC"].nunique()
final_cells = processed["CBC"].nunique()
print(f'Retained {final_cells} ({final_cells/starting_cells*100:.2f}%) CBCs')
starting_CBC_GBC_UMIs = counts.shape[0]
final_CBC_GBC_UMIs = processed.shape[0]
print(f'Retained {final_CBC_GBC_UMIs} ({final_CBC_GBC_UMIs/starting_CBC_GBC_UMIs*100:.2f}%) CBC-GBC-UMIs')

# Get CBC-GBC combos
df_combos = (
    processed
    .groupby(['CBC', 'GBC'])['UMI'].nunique()
    .to_frame('umi').reset_index()
    .assign(
        max_ratio=lambda x: \
        x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
        normalized_abundance=lambda x: \
        x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
    )
)

# Filter CBC-GBC combos and pivot
M, df_combos = filter_and_pivot(
    df_combos, 
    umi_treshold=5, 
    p_treshold=1, 
    max_ratio_treshold=.5, 
    normalized_abundance_treshold=.5
)
single_cbc = (M>0).sum(axis=1).loc[lambda x: x==1].index
M = M.loc[single_cbc]

# Final cell assignment
print(f'Assigned {single_cbc.size} ({single_cbc.size/starting_cells*100:.2f}%) cells')

# Correlation with known GBC reference
path_bulk = '/Users/IEO5505/Desktop/mito_bench/data/bulk_GBC_reference.csv'
bulk_df = pd.read_csv(path_bulk, index_col=0).query('sample=="AML_clones"')
pseudobulk_sc = M.sum(axis=0) / M.sum(axis=0).sum()
common = list(set(pseudobulk_sc.index) & set(bulk_df.index))

pseudobulk_sc = pseudobulk_sc.loc[common]
bulk = bulk_df.loc[common]['read_count'] / bulk_df.loc[common]['read_count'].sum()
corr = np.corrcoef(pseudobulk_sc, bulk)[0,1]


plt.scatter(sets.set_index('GBC_set').loc[common]['prevalence'], pseudobulk_sc)
plt.show()







sets = get_clones(M)
GBC_set = list(chain.from_iterable(sets['GBC_set'].map(lambda x: x.split(';')).to_list()))
redundancy = 1-np.unique(GBC_set).size/len(GBC_set)
occurrences = pd.Series(GBC_set).value_counts().sort_values(ascending=False)


##


# Visualize CBC-GBC filtering
fig, ax = plt.subplots(figsize=(5,5))
sns.scatterplot(data=df_combos, x='normalized_abundance', y='max_ratio', ax=ax)
fig.tight_layout()
plt.show()


