"""
Decontamination utils functions.
"""
import sys
import numpy as np
import umap
import matplotlib.pyplot as plt
from plotting_utils._plotting_base import *
from scipy.stats import pearsonr
import seaborn as sns

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
from helpers import *

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


def double_hist(measure_before, measure_after, title, value=None):
    """
    Plot two histograms of frequencies from tw list of a measure of the cells before and after decontamination
    """
    fig, axs = plt.subplots(1,2)
    axs[0].hist(measure_before, bins=10, color='skyblue', edgecolor='black')
    # Aggiunta di titoli e label agli assi
    axs[0].set_title('No decontamination')
    axs[0].set_xlabel(value)
    axs[0].set_ylabel('')

    axs[1].hist(measure_after, bins=10, color='skyblue', edgecolor='black')
    # Aggiunta di titoli e label agli assi
    axs[1].set_title( 'After decontX')
    axs[1].set_xlabel(value)
    axs[1].set_ylabel('')
    fig.suptitle(title)
    return fig

##


def UMI_GBC_variation(m_df, m_df_dec):

    'It returns the percentage of the variation in the UMi counts and GBC counts after decontamination'

    m_discrete = m_df.copy()
    m_discrete[m_discrete>0] = 1
    m_discrete_dec=m_df_dec.copy()
    m_discrete_dec[m_discrete_dec>0] = 1
    m_discrete_dec[m_df_dec<0] = 0
    diff_dis = m_discrete-m_discrete_dec
    dist_dis = []
    diff = m_df-m_df_dec
    dist = []
    assert (diff >= 0).all().all(), "The decontamination should only remove elements"
    assert (diff_dis >= 0).all().all(), "The decontamination should only remove elements"
    for i in range(m_df.shape[0]):
            dist_dis.append(sum(diff_dis.iloc[i]))
            dist.append(sum(diff.iloc[i]))
    UMI_variation = sum(dist)/(sum(sum(m_df.values[:])))
    GBC_variation = sum(dist_dis)/(sum(sum(m_discrete.values[:])))
    data=[GBC_variation, UMI_variation]
    return data


##


def scatter_correlation(other_feature, contamination,coverage_treshold, x_label= 'Other barcodes ratio',y_label='DecontX contamination ratio'):
    """
    Correlationa scatter plot of two set of measures changing the coverage filtering
    """
    fig = plt.figure()
    plt.scatter(other_feature, contamination, s=3)
    plt.gca().set_aspect('equal', adjustable='box')
    m, q = np.polyfit(other_feature, contamination, 1)
    plt.plot(np.array(other_feature), m*np.array(other_feature) + q, color='red')
    plt.xlim(0, 1)  
    plt.ylim(0, 1)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(f'Correlation with th= {coverage_treshold}')
    correlation_coef, p_value = pearsonr(other_feature, contamination)
    # Annotate plot with correlation coefficient and p-value
    plt.annotate(f'Correlation coefficient: {correlation_coef:.2f}\nP-value: {p_value:.2f}',
                xy = (np.mean(correlation_coef), np.mean(contamination)),  # Position of the annotation
                xytext = (0.4, 0.05),             # Offset of the text from the annotation position
                textcoords = 'axes fraction',  # Coordinate system of xytext
                fontsize = 8,     )         # Font size of the annotation text
                #arrowprops=dict(facecolor='black', arrowstyle='->'))  # Arrow properties
        
    return fig


##


def complete_corr_relevance(m_df, m_df_dec):
    """
    It returns the number of all the entrancies that have been changed by DecontX, 
    and the UMI normalized  cell abundancy of each GBC that has been completly removed by DecontX
    """
    total_correction = m_df_dec[(m_df!=m_df_dec)].count().sum()
    X = m_df.apply(lambda x: x/x.sum(), axis=1)
    level_of_correction = X[(m_df!=m_df_dec) & (m_df_dec==0)].values.flatten()
    nan_mask = np.isnan(level_of_correction)
    level_of_correction = level_of_correction[~nan_mask]

    return total_correction, level_of_correction


##



def create_boxplots(data_list, plot_title, y_label, positions = [0, 20, 40, 60, 80, 100]):
    """
    Create a grid of 6 boxplots from a list of data arrays.

    Parameters:
        data_list (list of arrays): List containing the data arrays for boxplots.

    Returns:
        fig
    """
    fig = plt.figure(figsize=(10, 6))
    sns.boxplot(data=data_list, width=0.5)
    plt.xticks(range(len(positions)), positions)
    plt.title(plot_title)
    plt.xlabel('Coverage threshold')
    plt.ylabel(y_label)
    plt.grid(True)

    return fig


##


def get_df_combos(m_df):
    """
    From a count matrix it get the combos dataframe with normalized abundance a max_ratio of each couple CBC-GBC
    """
    df_combos = m_df.stack().reset_index(name='umi').rename(columns={'level_0': 'CBC', 'level_1': 'GBC'})
    df_combos = df_combos[df_combos['umi']!=0]
    df_combos = df_combos.assign(
            max_ratio=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.max()),
            normalized_abundance=lambda x: \
            x.groupby('CBC')['umi'].transform(lambda x: x/x.sum())
        )
    return df_combos


##

