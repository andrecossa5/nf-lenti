"""
Decontamination utils functions.
"""

import numpy as np
import umap
import matplotlib.pyplot as plt


##


def merge_common_rows(df1, df2):
    """
    Returns two DataFrames containing only the rows that are common.
    """
    common_index = df1.index.intersection(df2.index)
    df1_filtered = df1.loc[common_index]
    df2_filtered = df2.loc[common_index]
    return df1_filtered, df2_filtered


##


def calculateUmap(m_df, n_neighbors=15):
    """
    Calculate 2D embeddings with UMAP.
    """
    X = m_df.apply(lambda x: x/x.sum(), axis=1)
    X = np.log2(X+1)
    reducer = umap.UMAP(min_dist=0.01, spread=1, n_neighbors=n_neighbors)
    embedding = reducer.fit_transform(X)
    return embedding


##


def plot_double_UMAP(m_df, m_df_dec,umap_dec,coverage_treshold, z=None):           
    """
    Plot the UMAP embeddings with cluster labels.
    """
    # Original UMAP
    m_embedding = calculateUmap(m_df)
    if z != None :
        m_clonal_labels = z
        m_dec_clonal_labels = z               
        m_n_clones = np.unique(z.values) 
        m_dec_n_clones = np.unique(z.values)       
    else :
        m_clonal_labels = calculateClonal_labels(m_df)
        m_n_clones = np.unique(m_clonal_labels)
        m_dec_clonal_labels = calculateClonal_labels(m_df_dec)
        m_dec_n_clones = np.unique(m_dec_clonal_labels)

    # Decontaminated UMAP
    m_dec_embedding = umap_dec.values
    
    # Visualization
    fig, axs = plt.subplots(1,2,figsize=(8,4))

    for i in m_n_clones:
       axs[0].scatter(m_embedding[m_clonal_labels == i, 0], m_embedding[m_clonal_labels == i, 1], label=i, s=3)
    axs[0].set_title('No decontamination')
    for i in m_dec_n_clones:
        axs[1].scatter(m_dec_embedding[m_dec_clonal_labels == i, 0], m_dec_embedding[m_dec_clonal_labels  == i, 1], label=i, s=3)
    axs[1].set_title('after DecontX')
    fig.suptitle(f'UMAP thr={coverage_treshold}')
    return fig


##
    

def calculateClonal_labels(m_df):
    clonal_labels = m_df.idxmax(axis=1).values
    return clonal_labels


##