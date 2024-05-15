#!/usr/bin/python

"""
New filtering procedure. With parallel pandas.
"""

import sys
import pickle
import numpy as np
import pandas as pd
from pandarallel import pandarallel
import argparse
from plotting_utils._utils import *
from plotting_utils._plotting_base import *
import csv
#sys.path.append("/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/RNA_decontamination")
#sys.path.append("/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/sc_gbc")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination")
sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from decontamination_utils import *
from helpers import *



##


# Parameters

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='preprocessing filtering ',
    description=
    """
    Script for preprocessing stats.
    """
)

my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Name of the sample. Default: None.'
)

my_parser.add_argument(
    '--read_elements_path', 
    type=str,
    default=None,
    help='Path to input tsv. Default: None.'
)

my_parser.add_argument(
    '--barcodes', 
    type=str,
    default=None,
    help='Path to the CBC barcodes associated to th input. Default: None.'
)

my_parser.add_argument(
    '--path_counts',
    type=str,
    default=None,
    help='pickle count path. Default: None.'
)

my_parser.add_argument(
    '--path_bulk', 
    type=str,
    default=None,
    help='Path to input bulk reference. Default: None.'
)

my_parser.add_argument(
    '--method', 
    type=str,
    default=None,
    help='Weinreb or Michael preprocessing method Default: None.'
)

my_parser.add_argument(
    '--modality', 
    type=str,
    default=None,
    help='Type of feature of the dataset. Default: None.'
)

my_parser.add_argument(
    '--correction_threshold', 
    type=int,
    default=6,
    help='threshold on the hamming distance for the feture filtering in Michael method. Default: 6.'
)

my_parser.add_argument(
    '--coverage_treshold',
    type=int,
    default=10,
    help='''
        Min number of reads to consider a CBC-GBC-UMI combination "true". Default: 10. This
        is used only if the --filtering method is not GMM (default) or medstd.
        '''
)

# Spikeins
my_parser.add_argument(
    '--umi_treshold',
    type=int,
    default=5,
    help='Min number of UMIs to consider a CBC-GBC combination supported. Default: 5.'
)

# p_treshold
my_parser.add_argument(
    '--p_treshold',
    type=float,
    default=1,
    help='Max p_poisson treshold to consider a CBC-GBC combination supported. Default: 1.'
)

# ratio_to_most_abundant_treshold
my_parser.add_argument(
    '--max_ratio_treshold',
    type=float,
    default=.5,
    help='Min ratio between a GBC nUMIs and the most abundant (nUMIs) GBC found for a given CBC. Default: .5.'
)

# Normalized abundance
my_parser.add_argument(
    '--normalized_abundance_treshold',
    type=float,
    default=.5,
    help='Min abundance (nUMIs fraction within a cell) of a CBC-GBC combination. Default: .5.'
)

# number of core for paralelization
my_parser.add_argument(
    '--ncore',
    type=int,
    default= 8 ,
    help='Number of cores'
)


##


# Parse arguments
args = my_parser.parse_args()
sample =args.sample
read_elements_path = args.read_elements_path
path_counts = args.path_counts
#method = args.method
barcodes_path = args.barcodes
modality = args.modality
correction_threshold = args.correction_threshold
coverage_treshold = args.coverage_treshold
umi_treshold = args.umi_treshold
p_treshold = args.p_treshold
max_ratio_treshold = args.max_ratio_treshold
normalized_abundance_treshold = args.normalized_abundance_treshold
path_bulk = args.path_bulk
n_cores = args.ncore


##


#sample = 'AML'
read_elements_path = '/Users/ieo6943/Documents/data/AML_clones/maester/grepped_head.txt'
#path_counts = '/Users/ieo6943/Documents/data/AML_clones/counts.pickle'
#method = "Micheals"
#modality = 'feature'
barcodes_path = '/Users/ieo6943/Documents/data/AML_clones/barcodes.tsv.gz'
#correction_threshold = 6
#coverage_treshold = 10
#umi_treshold = 5
#p_treshold = 1
#max_ratio_treshold = 0.5
#normalized_abundance_treshold = 0.5
#path_bulk = '/Users/ieo6943/Documents/data/AML_clones/bulk_reference.csv'
#n_cores = 10


##

if read_elements_path != None:
# Counts n reads per CBC-UMI-feature
    sc_df = pd.read_csv(read_elements_path, sep='\t', dtype='str')
    sc_df.columns = ['name', 'CBC', 'UMI', modality]
    if modality == 'GBC':
        d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
        sc_df[modality] = (
            sc_df[modality]
            .map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))
        )
    counts = count_UMIs(sc_df, gbc_col=modality)
    barcodes = pd.read_csv(barcodes_path)
    barcodes = barcodes.iloc[:,0].tolist()
    counts = counts[counts['CBC'].isin(barcodes)]
else:
    with open(path_counts, 'rb') as p:
        counts = pickle.load(p)['raw']

##
counts[counts['count']>4]
len(counts['feature'].iloc[0])
sc_df

# Stats
STATS = []
n_reads_total = counts['count'].sum()
median_n_reads_per_UMI = counts.groupby(['CBC', 'UMI'])['count'].sum().median()
print(f'n_reads_total: {n_reads_total}')   
print(f'median_n_reds per UMI: {median_n_reads_per_UMI}') 

# UMI redundancy
pandarallel.initialize(nb_workers=n_cores)
count_bad_UMI = (
    counts
    .groupby(['CBC'])
    .parallel_apply(lambda x: find_bad_UMI(x, modality))
    .droplevel(1).reset_index()
    .dropna()
    
)
print(f"mean of UMI with the same {modality} information {count_bad_UMI['bad_UMI_ratio'].values.mean()*100}%")
print(f"mean of UMI with the same {modality} information {count_bad_UMI['bad_UMI_ratio'].values.std()*100}%")
STATS.append(count_bad_UMI['bad_UMI_ratio'].values.mean()*100)
STATS.append(count_bad_UMI['bad_UMI_ratio'].values.std()*100)

##




df = counts[counts['CBC']=='AAACCCACACCAGTTA']



methods = ["Weinreb", "Micheals"]
for method in methods:
    # UMI filtering/processing
    t = Timer()
    t.start()

    if method == "Weinreb":
        print('Im using Weinreb')
        counts = mark_UMIs(counts, coverage_treshold=coverage_treshold)
        processed = counts.loc[lambda x: x['status']=='Retain']
    elif method == "Micheals":
        print('Im using Micheals')
        processed = (
            counts
            .groupby(['CBC', 'UMI'])
            .parallel_apply(lambda x: process_consensus_UMI(x,modality=modality, correction_threshold=correction_threshold))
            .droplevel(2).reset_index()
            .dropna()  
        )
    t_stop = t.stop()
    print(f'Elapsed {t_stop}')
    STATS.append(t_stop)


    ##


    # Before-after checks: CBCs and CBC-UMI-GBC retained

    # More stats to collect
    starting_cells = counts["CBC"].nunique()
    final_cells = processed["CBC"].nunique()
    print(f'Retained {final_cells} ({final_cells/starting_cells*100:.2f}%) CBCs')
    STATS.append(final_cells/starting_cells*100)
    starting_CBC_feature_UMIs = counts.shape[0]
    final_CBC_feature_UMIs = processed.shape[0]
    print(f'Retained {final_CBC_feature_UMIs} ({final_CBC_feature_UMIs/starting_CBC_feature_UMIs*100:.2f}%) CBC-GBC-UMIs')
    STATS.append(final_CBC_feature_UMIs/starting_CBC_feature_UMIs*100)

    ##


    if modality == 'GBC':

        # Get CBC-GBC combos
        df_combos = get_combos_not_filtered(processed)

        # Filter CBC-GBC combos and pivot
        M, df_combos = filter_and_pivot(                    # Argomenti fissi, per ora, da mettere ad inizio script
            df_combos, 
            umi_treshold=5,     
            p_treshold=1, 
            max_ratio_treshold=.5, 
            normalized_abundance_treshold=.5
        )

        # More stats
        print(f'Median MOI: {(M>0).sum(axis=1).median()}')
        STATS.append((M>0).sum(axis=1).median())
        # Final cell assignment
        single_cbc = (M>0).sum(axis=1).loc[lambda x: x==1].index 
        M = M.loc[single_cbc]
        print(f'Assigned {single_cbc.size} ({single_cbc.size/starting_cells*100:.2f}%) cells')
        STATS.append(single_cbc.size/starting_cells*100)

        ##


        # Correlation with known GBC reference if available (optional)           
        if os.path.exists(path_bulk):
            bulk_df = pd.read_csv(path_bulk, index_col=0).query(f'sample=="{sample}"')
            pseudobulk_sc = M.sum(axis=0) / M.sum(axis=0).sum()
            common = list(set(pseudobulk_sc.index) & set(bulk_df.index))
            pseudobulk_sc = pseudobulk_sc.loc[common]
            bulk = bulk_df.loc[common]['read_count'] / bulk_df.loc[common]['read_count'].sum()
            corr = np.corrcoef(pseudobulk_sc, bulk)[0,1]
            print(f'Correlation bulk read counts vs pseudobulk umi counts {corr:.2f}%) cells')
            STATS.append(corr)
        else:
            STATS.append(np.nan)

    ##

    else:
        STATS.append(np.nan, np.nan)

# Percorso del file CSV
STATS_path = f'/hpcnfs/home/ieo6943/mito_preprocessing/bin/RNA_decontamination/STATS_{sample}.csv'


# Scrivere la lista nel file CSV
with open(STATS_path, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(STATS)

# Create and save .csv
# Save STATS...
