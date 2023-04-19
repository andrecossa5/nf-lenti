#!/usr/bin/python

import os
import sys
import re
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from anndata import AnnData


##


def check_files(path_folder):
    """
    Check if all files are present.
    """
    patterns = ['_refAllele', 'C.txt', 'G.txt', 'A.txt', 'T.txt', 'depthTable', 'coverage']
    checklist = []
    for x in os.listdir(path_folder):
        checklist.append(any([ bool(re.search(p, x)) for p in patterns ]))
    assert all(checklist)


##


def sparse_from_long(df, cells_map, covariate, nrow, ncol):
    """
    Make a long df a sparse matrix. 
    """
    df['cell_id'] = df['cell'].map(lambda x: cells_map[x])
    df = df.loc[:, ['cell_id', 'pos', covariate]]
    df['pos'] = df['pos']-1 # index is 0-based. Positions are 1.
    rows, cols, vals = list(zip(*df.values.tolist()))
    matrix = csr_matrix((vals, (rows, cols)), shape=(nrow, ncol))

    return matrix


##


# Path
path_folder = sys.argv[1]

# Check all files
check_files(path_folder)

# Load coverage and reference alleles
cov = pd.read_csv(
    os.path.join(path_folder, 'maegatk.coverage.txt.gz'),
    header=None, 
    names=['pos', 'cell', 'cov']
)

cells = cov['cell'].unique()
cells_map = { c:i for i,c in enumerate(cells) }
ref_alleles = pd.read_csv(
    os.path.join(path_folder, 'chrM_refAllele.txt'),
    sep='\t', 
    header=None, 
    names=['pos', 'wt_allele']
)

# Get dimensions, for all matrices
nrow = cells.size
ncol = ref_alleles.shape[0]

# Make all sparse matrices
d = {}

# Coverage
d['cov'] = sparse_from_long(cov, cells_map, 'cov', nrow, ncol)

# Base counts and qualities
covariates = ['counts_fw', 'qual_fw', 'counts_rev', 'qual_rev']
base_files = { 
    x.split('.')[1] : x for x in os.listdir(path_folder) \
    if any(bool(re.search(p, x)) for p in ['A.txt', 'C.txt', 'T.txt', 'G.txt']) 
}
# Here we go
for k in base_files:
    df_base = pd.read_csv(
        os.path.join(path_folder, base_files[k]), 
        header=None, 
        names=['pos', 'cell'] + covariates
    )
    for covariate in covariates:
        d[f'{k}_{covariate}'] = sparse_from_long(df_base, cells_map, covariate, nrow, ncol)

# Build the AnnData
cells_meta = pd.Series(cells_map).to_frame().reset_index().iloc[:,[0]].set_index('index')
pos_meta = ref_alleles.set_index('pos')
afm = AnnData(X=d['cov'], obs=cells_meta, var=pos_meta)

# Fill afm
for k in d:
    afm.layers[k] = d[k]

# Save
afm.write('AFM.h5ad')