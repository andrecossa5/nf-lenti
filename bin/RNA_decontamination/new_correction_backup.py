"""
New filtering procedure. 
"""

import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting_utils._utils import *
from plotting_utils._plotting_base import *
from sklearn.metrics.pairwise import pairwise_distances
#sys.path.append("/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/RNA_decontamination")
#sys.path.append("/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/sc_gbc")

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from decontamination_utils import *
from helpers import *




##


# Path
# path_sc = ...
#path_counts = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/data/CTCs/clonal_info/PT_1_late.pickle'
path_main = '/Users/ieo6943/Documents/data/8_clones'
path_sc = os.path.join(path_main, 'GBC_read_elements.tsv.gz')
path_counts = os.path.join(path_main, 'counts.pickle')


# Parameters
correction_threshold = 3
umi_treshold = 5
p_treshold = 1
max_ratio_treshold = .5
normalized_abundance_treshold = .5
coverage_treshold = 100


##


# Read CBC-UMI-GBC triplets, and reverse complement GBCs
# sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
# sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
# d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
# sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))

# Count raw CBC-UMI-GBC reads
# counts = count_UMIs(sc_df) 

# Read, for the moment, from pickled counts
with open(path_counts, 'rb') as p:
    counts = pickle.load(p)['raw']


##


# Loop for each CBC/UMI combination
def process_consensus_UMI(df, correction_threshold=3):     # SPOSTARE in helpers, o simili, quando e' finita!
    """
    Process a CBC-UMI GBCs, and their read counts.
    """
    df = df.sort_values('count', ascending=False).set_index('GBC')
    GBCs = df.index.unique() # Should we get rid of the .unique()?

    if GBCs.size == 1:
        l = [ GBCs[0], df['count'].sum(), 1 ] 
    else:
        numeric_GBCs = np.array([ [ ord(char) for char in s ] for s in GBCs ]) 
        D = 18 * pairwise_distances(numeric_GBCs, metric='hamming')                     # pdist
        if (D<=correction_threshold).all():
            # Simplest alternative: majority barcode as consensus
            top_GBC = GBCs[0]                                   
            l = [ top_GBC, df['count'].sum(), df.loc[top_GBC, 'count']/df['count'].sum() ]
            # More advanced (but more computationally intensive...)
            """
            AAAC 40
            ACGC 35
            ACGT 30
            ACGA 25
            top_GBC = ACGC
            """
        else:
            l =  [np.nan] * 3

    res = pd.DataFrame(l).T
    res.columns = ['GBC', 'total_reads', 'consensus']

    return res


##


# Process each cell UMI into a consensus GBC sequence.
def f():
    """
    This is piece of code to optimize. Now embedded in a function f() just for time
    monitoring convenience.
    """
    processed = (
        counts#.head(100000)
        .groupby(['CBC', 'UMI'])
        .apply(lambda x: process_consensus_UMI(x)) # Need to be speeded up in parallel, if possible...
        # pandasparallel --> soluzione pandas-specifica
        # simpler and more general data structure + joblib parallelization --> soluzione general purpose
        .droplevel(2).reset_index()
        .dropna()                                  # Filter out CBC-UMIs with too heterogeneous GBC sets
    )
    return processed

# Run code and monitor time
processed = run_command(f, verbose=True)           # ~4 minutes. # ~12 minutes 8 clones


##


def find_bad_UMI(df, ham_th=1):
    UMIs = df['UMI'].unique()
    numeric_UMIs = np.array([ [ ord(char) for char in s ] for s in UMIs ]) 
    D = 12 * pairwise_distances(numeric_UMIs, metric='hamming')
    bad_UMI_count = ( df['GBC'].iloc[np.where(D==ham_th)[0]].values == df['GBC'].iloc[np.where(D==ham_th)[1]].values ).sum()
    bad_UMI_ratio = bad_UMI_count / D.shape[0]
    res = pd.DataFrame([bad_UMI_ratio])
    res.columns = ['bad_UMI_ratio']
    return res






##


#Micheals
#UMI-deduplication

count_bad_UMI = (
    processed
    .groupby(['CBC'])
    .apply(lambda x: find_bad_UMI(x))
    .droplevel(1).reset_index()
    .dropna()
    
)

count_bad_UMI['bad_UMI_ratio'].values.mean()
count_bad_UMI['bad_UMI_ratio'].values.std()

##

#Micheals
#cell loss

cell_loss = 1 - processed['CBC'].nunique()/counts['CBC'].nunique()



#Weinreb et al 

CBC_GBCs = df_combos[df_combos['status']=='supported'][['CBC', 'GBC']]
filtered = counts.loc[lambda x: x['status']=='Retain']
processed_W = get_combos(counts_Weinreb, gbc_col='GB')

M, processed_W = filter_and_pivot(
    processed_W, 
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold, 
    max_ratio_treshold=max_ratio_treshold, 
    normalized_abundance_treshold=normalized_abundance_treshold
)
processed_W_UMI = filtered.merge(processed_W[processed_W['status']=='supported'], on=['CBC', 'GBC'], how='inner')

count_bad_UMI_W = (
    processed_W_UMI
    .groupby(['CBC'])
    .apply(lambda x: find_bad_UMI(x))
    .droplevel(1).reset_index()
    .dropna()
    
)

count_bad_UMI_W['bad_UMI_ratio'].values.mean()

cell_loss_W = 1 - processed['CBC'].nunique()/counts['CBC'].nunique()

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
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold, 
    max_ratio_treshold=max_ratio_treshold, 
    normalized_abundance_treshold=normalized_abundance_treshold
)




"""
Extend process_consensus_UMI con metodo piu' fancy.
Quantificare UMI-duplication partedo dagli UMI Weinreb et al., o dai Micheals et al. --> Distribuzione (per cellula) della frazione di UMI<=1H distance.
N cellule e numero di UMI (CBC-UMI-GBC combinations) che rimangono (filtering UMI Weinreb et al., o dai Micheals et al.)
Sul campione a 8 GBC dimmi quanti GBC rimangono dopo Micheals. 

Date:
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold, 
    max_ratio_treshold=max_ratio_treshold, 
    normalized_abundance_treshold=normalized_abundance_treshold
in tre i casi (filtering UMI Weinreb et al. con treshold sul minimo, o dai Micheals et al.)
Quantifica MOI dalla M nei due casi, e il numero di cellule che si porterebbero a casa prendendo solo le singole.
"""










##


# Visualize CBC-GBC filtering
# fig, ax = plt.subplots(figsize=(5,5))
# sns.kdeplot(data=df_combos, x='normalized_abundance', y='max_ratio', ax=ax)
# fig.tight_layout()
# plt.show()

# Cell assignment, final
# unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
# filtered_M = M.loc[unique_cells]
# clones_df = get_clones(filtered_M)
# cells_df = (
#     filtered_M
#     .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
#     .to_frame('GBC')
# )

# Write out and save
# ...


##




