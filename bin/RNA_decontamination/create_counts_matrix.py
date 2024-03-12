"""
Creation of the counts matrix with differente coverage thresholds.
"""

import os
import sys
import pickle
import pandas as pd
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_contamination")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from decontamination_utils import *
from helpers import *


##


# Path
path_main = '/Users/ieo6943/Documents/data/8_clones'
path_sc = os.path.join(path_main, 'GBC_read_elements.tsv.gz')
path_count = os.path.join(path_main, 'counts.pickle')
path_bulk = None
path_sample_map = None


##


# Args
sample_params = None
correction_type = 'reference' 
sc_correction_treshold = 3
ncores = 8
filtering_method = "fixed"
coverage_treshold = 0 
umi_treshold = 5
p_treshold = 1
max_ratio_treshold = .8
normalized_abundance_treshold = .8
bulk_correction_treshold = 3
sample = None

##


# Read count from the path_sc
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
    # Read counts from the path_count, if available
    with open(path_count, 'rb') as p:
        data_pickle = pickle.load(p)
    COUNTS = data_pickle


##
    

# Create counts matrices with different coverage filters
for coverage_treshold in list(range(0, 101, 20)):
    counts = mark_UMIs(COUNTS[correction_type], coverage_treshold=coverage_treshold)
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
    M = df_combos.pivot_table(index='CBC', columns='GBC', values='umi')
    M[M.isna()] = 0
    M.to_csv(os.path.join(path_main,'count_matrix_&_DecontX', f'M_{coverage_treshold}.csv'))


##