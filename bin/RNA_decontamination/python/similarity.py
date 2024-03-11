
import sys
import pickle
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances, r2_score, mean_squared_error
import matplotlib.pyplot as plt
from plotting_utils._plotting_base import *
from mito_utils.kNN import *
from mito_utils.clustering import *

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
from helpers import *

def merge_common_rows(df1, df2):
    """
    Returns two DataFrames containing only the rows that are common.

    Arguments:
    df1, df2: Input DataFrames.

    Returns:
    Tuple containing the two filtered DataFrames.
    """
    # Find common indices
    common_index = df1.index.intersection(df2.index)
    
    # Select only rows with common indices
    df1_filtered = df1.loc[common_index]
    df2_filtered = df2.loc[common_index]
    
    return df1_filtered, df2_filtered

################################################################################################################

#path_sc = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/GBC_read_elements.tsv.gz'
path_count='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/counts.pickle'
sample_params=None
correction_type='reference' 
sc_correction_treshold=3
ncores=8
filtering_method="fixed"
umi_treshold=5
p_treshold=1
max_ratio_treshold=.8
normalized_abundance_treshold=.8

#path = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'
path = '/Users/ieo6943/Documents/Guido/NTS/'
coverage_treshold = 100
coverage_treshold_D = 0
path_count_pre = path + 'M_'+str(coverage_treshold_D)+'.csv'
path_count_dec = path +'M_'+str(coverage_treshold_D)+'_dec.csv'
path_count= path + 'counts.pickle'


##################################################################################################################

with open(path_count, 'rb') as p:
    data= pickle.load(p)


# Mark noisy UMIs, from chosen correction method
counts = mark_UMIs(data[correction_type], coverage_treshold=coverage_treshold)

# Get CBC-GBC combos with selected UMIs
df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')

# Filter CBC-GBC combos
M, _ = filter_and_pivot(
    df_combos, 
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold, 
    max_ratio_treshold=max_ratio_treshold,
    normalized_abundance_treshold=normalized_abundance_treshold
)

M.columns=list(range(M.shape[1]))
M_l = pd.DataFrame(M.idxmax(axis=1), columns=['GBC_set'])

m_df = pd.read_csv(path_count_pre, sep=',', index_col=0)
m_df_dec = pd.read_csv(path_count_dec, sep=',')
m_df_dec.columns.name = 'GBC'
m_df_dec.index=m_df.index

X = m_df_dec.apply(lambda x: x/x.sum(), axis=1)
X= np.log2(X+1)


similarity=[]
n_z_unique=[]
resolution=[]
for r in range(0, 15):
    res = r*1.5 + 0.05
    resolution.append(res)
    #res_dct={0:0.12,20:0.1, 40: 0.085, 60:0.09, 80: 0.1, 100: 0.1}
    knn_indices, distances, connectivities = kNN_graph(X, k=15, from_distances=False, nn_kwargs={})
    z = leiden_clustering(connectivities, res=res)
    n_z_unique.append(len(np.unique(z)))
    #assert n_z_unique==8, "we want the labels to be equal to the number of clone, m_dec "

    X_l = pd.DataFrame(z,index=X.index,columns=['GBC_set'])
    X_l_n, M_l_n = merge_common_rows(X_l, M_l)
    similarity.append(custom_ARI(X_l_n['GBC_set'], M_l_n['GBC_set']))


#size = list(map(lambda x: 3^(x), n_z_unique))
plt.figure()
plt.scatter(n_z_unique, similarity, s=5 ) #size)

# Show the plot
plt.xlabel('number of clusters')
plt.ylabel('Similarity')
plt.title('Data with Cossa filtering vs lognorm raw DecontX data & leiden clustering')
plt.legend()
plt.savefig(path+'similarity_scatter_lognorm_raw_DEC.png')

