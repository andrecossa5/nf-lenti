"""
Creation of the counts matrix with differente coverage thresholds
"""

import sys
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import umap
from sklearn.cluster import KMeans
import seaborn as sns
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
from helpers import *


##


# Path
path = '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/'
path_sc = path +'data/GBC_read_elements.tsv.gz'
path_count = path + 'data/counts.pickle'
path_bulk = None
path_sample_map = None


##


#parameters
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
bulk_correction_treshold = 3


##


#Read count from the path_sc
if path_bulk != None:
    # Reverse complement and value_counts
    sc_df, bulk_df = read_data(path_bulk, path_sample_map, path_sc, sample=sample)
    sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
    sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
    d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))

    ## Count
    COUNTS = {}
    COUNTS['raw'] = count_UMIs(sc_df)
    ## Correction
    sc_map = map_GBCs(sc_df, sc_correction_treshold=sc_correction_treshold)
    bulk_map = map_GBCs(sc_df, bulk_df, bulk_correction_treshold=bulk_correction_treshold, ncores=ncores)
    sc_df['GBC_reference-free'] = sc_df['GBC'].map(sc_map)
    sc_df['GBC_reference'] = sc_df['GBC'].map(bulk_map)
    COUNTS['reference-free'] = count_UMIs(sc_df, gbc_col='GBC_reference-free')
    COUNTS['reference'] = count_UMIs(sc_df, gbc_col='GBC_reference')

else:
    # Read counts from the pickle data pickle
    with open(path_count, 'rb') as p:
        data_pickle= pickle.load(p)
    COUNTS = data_pickle


##
    

for coverage_treshold in list(range(0, 101, 20)):
    # Mark noisy UMIs, from chosen correction method
    counts = mark_UMIs(COUNTS[correction_type], coverage_treshold=coverage_treshold)
    # Get CBC-GBC combos with selected UMIs
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
    #pivot table
    M =df_combos.pivot_table(index='CBC', columns='GBC', values='umi')
    M[M.isna()] = 0
    # Save M
    M.to_csv(path + f'/M_{coverage_treshold}.csv')

