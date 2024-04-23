"""
Creation of the counts matrix with differente coverage thresholds.
"""

import os
import sys
import pickle
import pandas as pd
import argparse
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
#sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination")
from decontamination_utils import *

#sys.path.append("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/sc_gbc")

from helpers import *


##

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
    '--path_counts',
    type=str,
    default=None,
    help='pickle count path. Default: None.'
)

my_parser.add_argument(
    '--path_results',
    type=str,
    default=None,
    help='results Default: None.'
)



##


# Parse arguments
args = my_parser.parse_args()
#sample =args.sample
path_sc = args.read_elements_path
path_counts = args.path_counts
path_results = args.path_results
#method = args.method



##


# Args
#sample_params = None
correction_type = 'raw' 
#sc_correction_treshold = 3
#ncores = 8
#coverage_treshold = 0 
#umi_treshold = 5
#p_treshold = 1
#max_ratio_treshold = .8
#normalized_abundance_treshold = .8
#bulk_correction_treshold = 3
#sample = '8_clones'

##

# Path
#sample = 'AML_clones'
#path_main = '/Users/ieo6943/Documents/data/'
#path_count = os.path.join(path_main,f'{sample}/counts.pickle')
#path_main = '/Users/ieo6943/Documents/data/complex_experiment'
#path_sc = '/Users/ieo6943/Documents/data/AML_clones/GBC_read_elements_scratch.tsv'
#path_count = path_counts
#path_bulk = None
#path_sample_map = None

##


# Read count from the path_sc
#if path_bulk != None:
if path_sc != None:
    print('sto aprendo il grepped')
    # Reverse complement and value_counts
    #sc_df, bulk_df = read_data(path_bulk, path_sample_map, path_sc, sample=sample)
    sc_df = pd.read_csv(path_sc, sep='\t', dtype='str')
    sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
    d_rev = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    sc_df['GBC'] = sc_df['GBC'].map(lambda x: ''.join([ d_rev[x] for x in reversed(x) ]))
    #sc_df.to_csv(os.path.join(path_main,'count_matrix_&_DecontX', 'sc_df.csv'))
    ## Count
    COUNTS = {}
    COUNTS['raw'] = count_UMIs(sc_df)
    print('creato la matrice di count')

    ## Correction
    #sc_map = map_GBCs(sc_df, sc_correction_treshold=sc_correction_treshold)
    #bulk_map = map_GBCs(sc_df, bulk_df, bulk_correction_treshold=bulk_correction_treshold, ncores=ncores)
    #sc_df['GBC_reference-free'] = sc_df['GBC'].map(sc_map)
    #sc_df['GBC_reference'] = sc_df['GBC'].map(bulk_map)
    #COUNTS['reference-free'] = count_UMIs(sc_df, gbc_col='GBC_reference-free')
    #COUNTS['reference'] = count_UMIs(sc_df, gbc_col='GBC_reference')

else:
    # Read counts from the path_count, if available
    with open(path_counts, 'rb') as p:
        data_pickle = pickle.load(p)
    COUNTS = data_pickle


##

# % of raw counts lost with bulk
#assert COUNTS['raw']['count'].sum() == COUNTS['reference-free']['count'].sum()
#perc_read_retained_bulk = COUNTS['reference']['count'].sum() / COUNTS['raw']['count'].sum()
#
## % GBCs retained with sc and bulk corrections
#perc_gbc_retained_sc = COUNTS['reference-free']['GBC_reference-free'].value_counts().size / \
#                        COUNTS['raw']['GBC'].value_counts().size
#perc_gbc_retained_bulk = COUNTS['reference']['GBC_reference'].value_counts().size / \
#                        COUNTS['raw']['GBC'].value_counts().size
    

# Viz correction effect
fig, axs = plt.subplots(1,1,figsize=(5,5))

viz_UMIs(COUNTS['raw'], axs[0])
axs[0].set(title='Raw')
axs[0].text(.53, .95, f'Total reads: {COUNTS["raw"]["count"].sum():.2e}', transform=axs[0].transAxes)
axs[0].text(.53, .91, f'Total GBCs: {COUNTS["raw"]["GBC"].unique().size:.2e}', transform=axs[0].transAxes)
axs[0].text(.53, .87, f'n reads: {COUNTS["raw"]["count"].median():.2f} (+-{COUNTS["raw"]["count"].std():.2f})', 
            transform=axs[0].transAxes)

#viz_UMIs(COUNTS['reference-free'], axs[1])
#axs[1].set(title='Reference-free correction')
#axs[1].text(.53, .95, f'% reads retained: 100%', transform=axs[1].transAxes)
#axs[1].text(.53, .91, f'% GBCs retained: {perc_gbc_retained_sc*100:.2f}%', transform=axs[1].transAxes)
#axs[1].text(.53, .87, f'n reads: {COUNTS["reference-free"]["count"].median():.2f} (+-{COUNTS["reference-free"]["count"].std():.2f})', 
#            transform=axs[1].transAxes)
#
#viz_UMIs(COUNTS['reference'], axs[2])
#axs[2].set(title='Bulk-DNA reference correction (PT)')
#axs[2].text(.53, .95, f'Reads retained: {perc_read_retained_bulk*100:.2f}%', transform=axs[2].transAxes)
#axs[2].text(.53, .91, f'GBCs retained: {perc_gbc_retained_bulk*100:.2f}%', transform=axs[2].transAxes)
#axs[2].text(.53, .87, f'n reads: {COUNTS["reference"]["count"].median():.2f} (+-{COUNTS["reference"]["count"].std():.2f})', 
#           transform=axs[2].transAxes)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'CBC_GBC_UMI_read_distribution.png'), dpi=300)


    ##



# Create counts matrices with different coverage filters
#for coverage_treshold in list(range(0, 61, 12)):
#    counts = mark_UMIs(COUNTS[correction_type], coverage_treshold=coverage_treshold)
#    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')
#    M = df_combos.pivot_table(index='CBC', columns='GBC', values='umi')
#    M[M.isna()] = 0
#    M.to_csv(os.path.join(path_results, f'M_{coverage_treshold}.csv'))
    

##