#!/usr/bin/python

import os
from scipy.sparse import csr_matrix, save_npz
import pandas as pd


##


def read_call(path):
    cell = path.split('_')[0]
    df = pd.read_csv(path, sep='\t', header=None)
    df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'AD', 'DP']
    df['AD'] = df['AD'].map(lambda x: int(x.split(',')[1]))
    df = df.assign(cell=cell)
    return df  


##


def main():

    calls = [ read_call(x) for x in os.listdir() if x.endswith('filtered.tsv') ]
    calls = pd.concat(calls)
    cell_map = { x:i for i,x in enumerate(calls['cell'].unique()) }
    calls['cell'] = calls['cell'].map(cell_map).astype(int)

    pd.Series(cell_map).to_csv('cells.txt', header=False)
    save_npz(csr_matrix(calls[['POS','cell','AD']].values), 'AD.npz')
    save_npz(csr_matrix(calls[['POS','cell','DP']].values), 'AD.npz')


##


if __name__ == '__main__':
    main()


