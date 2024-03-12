"""
Dimensionality reduction script.
"""

import os
import sys
import numpy as np
import pandas as pd
from mito_utils.kNN import *
from mito_utils.clustering import *
from plotting_utils._plotting_base import *
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination")
from helpers import *
from decontamination_utils import *


##


# Paths
path = '/Users/ieo6943/Documents/Guido/sAML1'
path_count_pre = os.path.join(path, 'M.csv')
path_count_dec = os.path.join(path, 'M_dec.csv')
path_count = '/Users/ieo6943/Documents/Guido/counts.pickle'


##


# Args
coverage_treshold = 10
coverage_treshold_D = 0
sample_params = None
correction_type = 'reference' 
sc_correction_treshold = 3
ncores = 8
filtering_method = "fixed"
umi_treshold = 5 
p_treshold = 1
max_ratio_treshold = .8
normalized_abundance_treshold = .8


##


# Read counts matrices and DecontX clusters
m_df = pd.read_csv(path_count_pre, index_col=0)
m_df_dec = pd.read_csv(path_count_dec, index_col=0)
z = pd.read_csv(os.path.join(path, 'z.csv'), index_col=0)
fig = plot_double_UMAP(m_df, m_df_dec)
# fig.savefig(os.path.join(path, 'UMAP_pre_post_decontX.png'), dpi=300)

# Filtering and pivot
counts = mark_UMIs(data[correction_type], coverage_treshold=coverage_treshold)
df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
M, _ = filter_and_pivot(
    df_combos, 
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold, 
    max_ratio_treshold=max_ratio_treshold,
    normalized_abundance_treshold=normalized_abundance_treshold
)
M.columns = list(range(M.shape[1]))
M_l = pd.DataFrame(M.idxmax(axis=1), columns=['GBC_set'])
m_df_dec.columns.name = 'GBC'
m_df_dec.index = m_df.index

# Log normalization
X = m_df_dec.apply(lambda x: x/x.sum(), axis=1)
X = np.log2(X+1)

# Calculate similarity
similarity = []
n_z_unique = []
resolution = []
for r in range(0, 15):
    res = r*1.5 + 0.05
    resolution.append(res)
    knn_indices, distances, connectivities = kNN_graph(X, k=15)
    z = leiden_clustering(connectivities, res=res)
    n_z_unique.append(len(np.unique(z)))
    X_l = pd.DataFrame(z,index=X.index,columns=['GBC_set'])
    X_l_n, M_l_n = merge_common_rows(X_l, M_l)
    similarity.append(custom_ARI(X_l_n['GBC_set'], M_l_n['GBC_set']))

# Change with fig, ax = plt.subplots() synthax
plt.figure()
plt.scatter(n_z_unique, similarity, s=5 )
# Show the plot
plt.xlabel('number of clusters')
plt.ylabel('Similarity')
plt.title('Data with Cossa filtering vs lognorm raw DecontX data & leiden clustering')
plt.legend()
plt.savefig(path+'similarity_scatter_lognorm_raw_DEC.png')


##
