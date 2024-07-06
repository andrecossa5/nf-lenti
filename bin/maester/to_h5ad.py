#!/usr/bin/python

import os
import sys
import re
import pandas as pd
from scipy.sparse import csr_matrix
from anndata import AnnData
import pysam


##


def check_files(path_folder='tables'):
    """
    Check if all files are present.
    """
    patterns = ['C.txt.gz', 'G.txt.gz', 'A.txt.gz', 'T.txt.gz', 'coverage.txt.gz']
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


# Run
def main():

    # Read ref_alleles ad format as df
    pattern = sys.argv[1]
    sample = sys.argv[2]
    fasta = pysam.FastaFile(f'{pattern}.fa')
    mt_genome = fasta.fetch(pattern)
    ref_alleles = pd.DataFrame(
        [ (i+1,x) for i, x in enumerate(list(mt_genome)) ], 
        columns=['pos', 'ref']
    ).set_index('pos')

    # Check all tables in place
    check_files()

    # Load coverage and reference alleles
    cov = pd.read_csv(
        os.path.join('tables', 'coverage.txt.gz'),
        header=None, 
        names=['pos', 'cell', 'cov']
    )
    cells = cov['cell'].unique()
    cells_map = { c:i for i,c in enumerate(cells) }

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
        x.split('.')[0] : x for x in os.listdir('tables') \
        if any(bool(re.search(p, x)) for p in ['A.txt', 'C.txt', 'T.txt', 'G.txt']) 
    }

    # Here we go
    for k in base_files:
        df_base = pd.read_csv(
            os.path.join('tables', base_files[k]), 
            header=None, 
            names=['pos', 'cell'] + covariates
        )
        for covariate in covariates:
            d[f'{k}_{covariate}'] = sparse_from_long(df_base, cells_map, covariate, nrow, ncol)

    # Build the AnnData
    cells_meta = pd.Series(cells_map).to_frame().reset_index().iloc[:,[0]].set_index('index')
    afm = AnnData(X=d['cov'], obs=cells_meta, var=ref_alleles)
    afm.obs_names = afm.obs_names.map(lambda x: f'{x}_{sample}')

    # Fill afm
    for k in d:
        afm.layers[k] = d[k]

    # Save
    afm.write('AFM.h5ad')


##


# Run 
if __name__ == '__main__':
    main()


##