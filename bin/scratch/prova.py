# 


import os
import sys


#


# Args
# sample = 'AC_AC_mets_3'
# path_bulk = '/Users/IEO5505/Desktop/example_mito/scratch_data/'
# path_sample_map = '/Users/IEO5505/Desktop/example_mito/scratch_data/sample_map.csv'
# path_sc = '/Users/IEO5505/Desktop/example_mito/scratch_data/GBC_read_elements.tsv.gz'
# ncores = 8
# bulk_correction_treshold = 3
# sc_correction_treshold = 3
# correction_type = 'reference-free'
# filtering_method = 'fixed'
# coverage_treshold = 10
# umi_treshold = 5
# p_treshold = 0.5
# max_ratio_treshold = .5
# max_ratio_treshold = .5


# Import code
sys.path.append('/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/sc_gbc')
os.chdir('/Users/IEO5505/Desktop/example_mito/scratch_results')
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from helpers import *

from mito_utils.phylo import build_tree


cell_assignment_workflow(
    path_bulk, 
    path_sample_map, 
    path_sc, 
    sample, 
    correction_type=correction_type, 
    sc_correction_treshold=sc_correction_treshold,
    bulk_correction_treshold=bulk_correction_treshold,
    ncores=ncores,
    filtering_method=filtering_method,
    coverage_treshold=coverage_treshold,
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold,
    max_ratio_treshold=max_ratio_treshold,
    normalized_abundance_treshold=normalized_abundance_treshold,
    sample_params=params
)
