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
    calls.to_csv('samtools_allele_table.csv.gz')


##


if __name__ == '__main__':
    main()


