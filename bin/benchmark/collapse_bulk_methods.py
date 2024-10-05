#!/usr/bin/python

import os
import pandas as pd


##


def read_call(path):
    cell = path.split('_')[0]
    df = pd.read_csv(path, sep='\t', header=None)
    df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'AD', 'DP']
    df = df.loc[df['REF'].map(lambda x: len(x)==1)]                  # Only SNVs here
    df['AD'] = df['AD'].map(lambda x: int(x.split(',')[1]))
    df = df.assign(cell=cell)
    return df  


##


def main():

    calls = []
    for x in os.listdir():
        try:
            if x.endswith('filtered.tsv'):
                calls.append(read_call(x))
        except:
            pass
        
    calls = pd.concat(calls).query('DP>=10') # At least 10 deduplicated reads
    calls.to_csv(f'allele_table.csv.gz')


##


if __name__ == '__main__':
    main()


