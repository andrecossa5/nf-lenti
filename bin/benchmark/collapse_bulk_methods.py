#!/usr/bin/python

import os
import sys
import pandas as pd


##


def read_samtools_call(path):
    cell = path.split('_')[0]
    df = pd.read_csv(path, sep='\t', header=None)
    df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'AD', 'DP']
    df['AD'] = df['AD'].map(lambda x: int(x.split(',')[1]))
    df = df.assign(cell=cell)
    return df  


##


# Args 
method = sys.argv[1] 


##


def main():

    if method == 'samtools':
        calls = [ read_samtools_call(x) for x in os.listdir() if x.endswith('filtered.tsv') ]

    calls = pd.concat(calls)
    calls.to_csv(f'{method}_allele_table.csv.gz')


##


if __name__ == '__main__':
    main()


