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
#path_main = '/Users/ieo6943/Documents/data/8_clones/'
path_main = '/Users/ieo6943/Documents/data/complex_experiment'
path_decontX = os.path.join(path_main,'count_matrix_&_DecontX')
path_count = os.path.join(path_main,'counts.pickle')
#path_results = '/Users/ieo6943/Documents/results/8_clones/'
path_results = '/Users/ieo6943/Documents/results/complex_experiment/'
th_list = list(range(0, 61, 12))


##


# Read counts matrices and DecontX clusters
for coverage_treshold in th_list:
    m_df = pd.read_csv(os.path.join(path_decontX, f'M_{coverage_treshold}.csv'), index_col=0)
    m_df_dec = pd.read_csv(os.path.join(path_decontX, f'M_{coverage_treshold}_dec.csv'), index_col=0)
    z = pd.read_csv(os.path.join(path_decontX, f'z{coverage_treshold}.csv'))
    umap_dec = pd.read_csv(os.path.join(path_decontX, f'umap_{coverage_treshold}.csv'))
    fig = plot_double_UMAP(m_df, m_df_dec, umap_dec,coverage_treshold=coverage_treshold, z=None )
    fig.savefig(os.path.join(path_results, f'UMAP_pre_post_decontX_{coverage_treshold}.png'), dpi=300)

    ##