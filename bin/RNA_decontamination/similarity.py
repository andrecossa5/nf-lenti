"""
Script to evaluate similarity between raw matrix (raw and decontaminated), with the filtered count matrix
"""
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting_utils._plotting_base import *
from mito_utils.kNN import *
from mito_utils.clustering import *

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from decontamination_utils import *
from helpers import *


##


# Paths
path_main = '/Users/ieo6943/Documents/data/8_clones/'
path_decontX = os.path.join(path_main,'count_matrix_&_DecontX')
path_count = os.path.join(path_main,'counts.pickle')
path_results = '/Users/ieo6943/Documents/results/8_clones/'


##


# Parameters
correction_type='reference' 
umi_treshold=5
p_treshold=1
max_ratio_treshold=.8
normalized_abundance_treshold=.8
coverage_treshold = 100
coverage_treshold_raw = 0


##

## uploading of the count matrixs
m_df = pd.read_csv(os.path.join(path_decontX, f'M_{coverage_treshold_raw}.csv'), index_col=0) 
m_df_dec = pd.read_csv(os.path.join(path_decontX, f'M_{coverage_treshold_raw}_dec.csv'), index_col=0) 
with open(path_count, 'rb') as p:
    COUNTS= pickle.load(p)

# Filtering and pivot
counts = mark_UMIs(COUNTS[correction_type], coverage_treshold=coverage_treshold)
df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
M, _ = filter_and_pivot(
    df_combos, 
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold, 
    max_ratio_treshold=max_ratio_treshold,
    normalized_abundance_treshold=normalized_abundance_treshold
)
M_l = pd.DataFrame(M.idxmax(axis=1), columns=['GBC_set'])
m_df_dec.columns.name = 'GBC'
keys = ['raw', 'decontaminated']
similarity = {keys[0]: [], keys[1]: []}
n_z_unique = {keys[0]: [], keys[1]: []}
resolution = {keys[0]: [], keys[1]: []}
df = {keys[0]: m_df, keys[1]: m_df_dec}

for type in keys:

    # Log normalization
    X = df[type].apply(lambda x: x/x.sum(), axis=1)
    X = np.log2(X+1)

    # Calculate similarity

    for r in range(0, 15):
        res = r*0.2 + 0.05
        resolution[type].append(res)
        knn_indices, distances, connectivities = kNN_graph(X, k=15)
        z = leiden_clustering(connectivities, res=res)
        n_z_unique[type].append(len(np.unique(z)))
        X_l = pd.DataFrame(z,index=X.index,columns=['GBC_set']) 
        X_l_n, M_l_n = merge_common_rows(X_l, M_l)
        similarity[type].append(custom_ARI(X_l_n['GBC_set'], M_l_n['GBC_set']))

    fig, ax = plt.subplots() 
    plt.figure()
    ax.scatter(n_z_unique[type], similarity[type], s=5 )
    ax.set_xlabel('number of clusters')
    ax.set_ylabel('Similarity')
    ax.set_title(f'Filtered data vs lognorm {type} data & leiden clustering')
    #ax.set_legend()
    fig.savefig(os.path.join(path_results, f'similarity_scatter_lognorm_{type}.png'), dpi=300)


##