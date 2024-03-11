"""
Utils for custom cell_assignement. Brings together elements from LARRY_2020, CR_2023 and
other (pre-)pubblication works.
g
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances, r2_score, mean_squared_error
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from scipy.stats import poisson
from itertools import chain
from plotting_utils._plotting_base import *
from plotting_utils._utils import Timer
import umap
import scipy.sparse as sp
from sklearn.cluster import KMeans
import seaborn as sns
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")

from helpers import *

# Parametri scratch reference free
path='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/'
#path_sc = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/GBC_read_elements.tsv.gz'
#path_sc = '/Users/ieo6943/Documents/Guido/export_cluster/GBC_read_elements.tsv.gz'
path_count='/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/data/counts.pickle'
#path_empty = '/Users/ieo6943/Documents/Guido/empty_bead/'
#path_count='/Users/ieo6943/Documents/Guido/NTS/counts.pickle'
sample_params=None
correction_type='reference' 
sc_correction_treshold=3
ncores=8
filtering_method="fixed"
coverage_treshold=0 
umi_treshold=5
p_treshold=1
max_ratio_treshold=.8
normalized_abundance_treshold=.8
   
# Reverse complement and value_counts
#sc_df, bulk_df = read_data(path_bulk, path_sample_map, path_sc, sample=sample)
#sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
#sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
#d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
#sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))
#
## Count
#COUNTS = {}
#COUNTS['raw'] = count_UMIs(sc_df)
#print('Correction')
## Correction
#sc_map = map_GBCs(sc_df, sc_correction_treshold=sc_correction_treshold)
#sc_df['GBC_reference-free'] = sc_df['GBC'].map(sc_map)
#COUNTS['reference-free'] = count_UMIs(sc_df, gbc_col='GBC_reference-free')


# Read counts as pickle
with open(path_count, 'rb') as p:
    data_pickle= pickle.load(p)

median_data={}

COUNTS = data_pickle


data = data_pickle
for coverage_treshold in list(range(0, 101, 20)):
    print(coverage_treshold)
    # Mark noisy UMIs, from chosen correction method
    counts = mark_UMIs(data[correction_type], coverage_treshold=coverage_treshold)
    #count = mark_UMIs(data_pickle[correction_type], coverage_treshold=coverage_treshold)

    #counts['CBC'].nunique()
    #counts[counts['status']=='Retain']
    # Get CBC-GBC combos with selected UMIs

    # Get CBC-GBC combos with selected UMIs
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')

    #plt.figure(figsize=(12, 6))


    sns.kdeplot(data=df_combos, x=df_combos['max_ratio'], y=df_combos['normalized_abundance'], label='Before decontamination')

    plt.title(f'Density Plot th = {coverage_treshold}')


    # Salva il plot
    #plt.savefig(path_empty + f'density_plot_empty_bead{coverage_treshold}')
    #sparse_matrix = sp.coo_matrix((df['umi'], (df['CBC'].astype('category').cat.codes, df['GBC'].astype('category').cat.codes)))
    #M=sparse_matrix
    M =df_combos.pivot_table(index='CBC', columns='GBC', values='umi')
    #M = pd.crosstab(df_combos['CBC'], df_combos['GBC'], values=df_combos['umi'], aggfunc='sum', sparse=True)

    M[M.isna()] = 0
    
    # Filter CBC-GBC combos
    #M, _ = filter_and_pivot(
    #   df_combos, 
    #  umi_treshold=umi_treshold, 
    # p_treshold=p_treshold, 
    #max_ratio_treshold=max_ratio_treshold,
    #normalized_abundance_treshold=normalized_abundance_treshold
    #)

    # Save M
    M.to_csv(path + f'/M_{coverage_treshold}.csv')
    #np.savetxt('sparse_matrix.csv', M.toarray(), delimiter=',', fmt='%d')

########################################
type(M.iloc[:10,:5])
type(M)

dir(M)
M.values
column_names = M.columns.tolist()

X = M[column_names]  # Extract the feature columns
reducer = umap.UMAP()  # Initialize UMAP reducer
embedding = reducer.fit_transform(X)  # Perform UMAP dimensionality reduction


plt.scatter(embedding[:, 0], embedding[:, 1], s=10)  # Scatter plot of UMAP embeddings
plt.title('UMAP Plot')
plt.xlabel('UMAP Component 1')
plt.ylabel('UMAP Component 2')
plt.show()

# Step 2: Perform clustering on the UMAP embeddings
num_clusters = 8  # Number of clusters you want
kmeans = KMeans(n_clusters=num_clusters)
cluster_labels = kmeans.fit_predict(embedding)

# Step 3: Plot the UMAP embeddings with cluster labels
plt.figure(figsize=(10, 8))
for i in range(num_clusters):
    plt.scatter(embedding[cluster_labels == i, 0], embedding[cluster_labels == i, 1], label=f'Cluster {i+1}')
plt.title('UMAP Clustering')
plt.xlabel('UMAP Dimension 1')
plt.ylabel('UMAP Dimension 2')
plt.legend()
plt.savefig('umap_clusters.png')
plt.show()


