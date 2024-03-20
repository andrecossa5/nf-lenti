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
correction_type='reference' 
umi_treshold=5
p_treshold=1
max_ratio_treshold=.5
normalized_abundance_treshold=.5
th_list = list(range(0, 101, 20))


##


with open(path_count, 'rb') as p:
    COUNTS= pickle.load(p)


sc_df = pd.read_csv(path_sc, sep='\t', header=None, dtype='str')
sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))

## Count
COUNTS = {}
COUNTS['raw'] = count_UMIs(sc_df)

counts = mark_UMIs(COUNTS['raw'], coverage_treshold=0)
unique_GBC_df = counts.groupby(['CBC', 'UMI'])['GBC'].unique().to_frame('GBC').reset_index()


df_ = counts.query('CBC=="AAACCCATCAACACGT" and UMI=="AAAACCACACAC"').sort_values('count', ascending=False)

unique_GBC_df.iloc[0]
str_to_int = np.array([[ord(char) for char in s] for s in df_['GBC']])
distance = 18*pairwise_distances(str_to_int, metric='hamming')
distance.sort(axis=1)
distance[:,1]
distance.shape[0]

CBC_UMI_dct = {}
CBC_UMI_good = {}
CBC_UMI_to_correct = {}
a = 0
b = 0
for index, GBC_list in enumerate(unique_GBC_df['GBC']):
    print(index+1, ' su ', unique_GBC_df.shape[0])
    str_to_int = np.array([[ord(char) for char in s] for s in GBC_list])
    distance = 18*pairwise_distances(str_to_int, metric='hamming')
    if np.any(distance>1):
        CBC_UMI_dct[(unique_GBC_df['CBC'][index], unique_GBC_df['UMI'][index])] = None
    else:
        
        a +=1
        #print('-----a-----',a)
        CBC = unique_GBC_df['CBC'][index]
        UMI = unique_GBC_df['UMI'][index]

        if len(GBC_list)==1:

            
            CBC_UMI_good[(CBC, UMI)] = GBC_list, 1

        else:
            #element = counts[counts[['CBC','UMI']]].apply(tuple, axis=1).isin([(CBC, UMI)])['GBC']
            element= counts[(counts['CBC'] == CBC ) & (counts['UMI']== UMI)]
            GBC_major = element['GBC'][element['count'].idxmax()]
            count_major_density = element['count'].max()/element['count'].sum()
            CBC_UMI_to_correct[(CBC, UMI)] = GBC_major, count_major_density
            b += 1
            #print('-----b----------',b)
    #print(index)
    if distance.shape[0]==10:
        break



CBC_UMI_good[('AACAAGATCGCCAATA', 'AACGTAAGCCGT')][0]


a = np.array(CBC_UMI_good.items)

stringhe = CBC_UMI_dct[('AAACCCACATACTGAC', 'AGACTGACCATT')][0]
str_to_int = np.array([[ord(char) for char in s] for s in stringhe])
str_to_int.shape
distance = 18*pairwise_distances(str_to_int, metric='hamming')
bad_CBC_GBC_UMI = ((distance>1).sum(axis=1))

CBC_UMI_good.keys()

counts_filtered = counts[counts[['CBC','UMI']]].apply(tuple, axis=1).isin(CBC_UMI_good.keys())['GBC']


# DataFrame di esempio
data = {'A': [1, 2, 3, 4, 5],
        'B': [10, 20, 30, 40, 50],
        'C': [100, 200, 300, 400, 500],
        'D': [1000, 2000, 3000, 4000, 5000]}
df = pd.DataFrame(data)

# Lista di liste di valori per filtrare la colonna 'D'
valori_da_cercare = [[1, 10, 100], [3, 30, 300]]

# Converte la lista di liste in una lista di tuple per poterla passare a isin()
valori_da_cercare = [tuple(x) for x in valori_da_cercare]

# Filtra il DataFrame sulla base della colonna 'A', 'B' e 'C' usando isin() e quindi seleziona solo la colonna 'D'
df_filtrato = df[df[['A', 'B', 'C']].apply(tuple, axis=1).isin(valori_da_cercare)]['D']

print(df_filtrato)