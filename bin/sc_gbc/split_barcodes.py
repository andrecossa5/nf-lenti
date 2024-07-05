#!/usr/bin/python

"""
Split a .tsv.gz file of barcodes into chuncks of 3000 barcodes.
"""

import sys
import pandas as pd


##


# Args
path_cbs = sys.argv[1]
chunk_size = int(sys.argv[2])

# path_cbs = '/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/toy/allowed_cbs.csv'
# chunk_size = 3000


##


def main():
    """
    Split CBs into multiple .csv files.
    """

    cbs = pd.read_csv(path_cbs, sep='\t', header=None).iloc[:,0]
    n_cbs = cbs.size

    for file_index, i in enumerate(range(0,n_cbs,chunk_size)):
        cbs[i:i+chunk_size].to_csv(f'barcodes_{file_index}.csv', index=False, header=False)


# Run
if __name__ == "__main__":
    main()
