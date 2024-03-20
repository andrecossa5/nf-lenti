"""
New filtering procedure.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting_utils._plotting_base import *
import umap
from mito_utils.kNN import *
from mito_utils.clustering import *
from sklearn.metrics.pairwise import pairwise_distances

sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from decontamination_utils import *
from helpers import *


##


# Path
path_main = '/Users/ieo6943/Documents/data/8_clones'
#path_main = '/Users/ieo6943/Documents/data/complex_experiment'
path_sc = os.path.join(path_main, 'GBC_read_elements.tsv.gz')
path_count = os.path.join(path_main, 'counts.pickle')
#path_results = '/Users/ieo6943/Documents/results/complex_experiment/'
path_bulk = None
path_sample_map = None


##


# Parameters
correction_threshold = 1
correction_type = 'reference' 
umi_treshold = 5
p_treshold = 1
max_ratio_treshold = .5
normalized_abundance_treshold = .5
th_list = list(range(0, 101, 20))


##


# Read 
# with open(path_count, 'rb') as p:
#     COUNTS= pickle.load(p)

# Rev complement GBC
sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))

# Count CBC


# Split CBC/UMI reads
D = {}
for i in range(sc_df.shape[0]):
    _, CBC, UMI, GBC = sc_df.iloc[i,:].values
    key = (CBC, UMI)
    if key not in D:
        D[key] = [GBC]
    else:
        D[key].extend(GBC)


##


# Loop for each CBC/UMI combination
def process_CBC_UMI(GBC_list, correction_threshold=2):

    if np.unique(GBC_list).size == 1:
        l = [ GBC_list[0], len(GBC_list), 1 ] 
    else:
        str_to_int = np.array([ [ ord(char) for char in s ] for s in GBC_list ])
        distance = 18 * pairwise_distances(str_to_int, metric='hamming')
        if np.any(distance<=correction_threshold):
            counts = pd.Series(GBC_list).value_counts()
            l = [ counts.index[0], counts.sum(), counts.values[0]/counts.sum() ] 

    return l


# Here we go
CBC_UMI_processed = {}
count = 0
for key in D:    
    print(f'{count+1}/{len(D)}')
    GBC_list = D[key]
    CBC_UMI_processed[key] = process_CBC_UMI(GBC_list, correction_threshold=2)
    if count == 10000:
        break


##
    