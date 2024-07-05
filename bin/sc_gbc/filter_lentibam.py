#!/usr/bin/python

"""
Filter a bam file for records having the needed CBs.
"""

import os
import sys
import pandas as pd


##


# Args
i = sys.argv[1]
o = sys.argv[2]
path_cbs = sys.argv[3]

# i = 'toy_dataset.bam'
# o = 'output.bam'
# path_cbs = 'allowed_cbs.csv'


##


def main():
      """
      Filter all records from a set of CBs.
      """

      cbs = pd.read_csv(path_cbs, sep='\t', header=None).iloc[:,0].to_list()
      picard_call = f'picard FilterSamReads -I {i} -O {o} -TAG CB ' + \
                    ' '.join([ f'-TV {x}' for x in cbs ]) + \
                    ' -FILTER includeTagValues' \
                    ' --QUIET true' \
                    ' --COMPRESSION_LEVEL 1' \
                    ' --MAX_RECORDS_IN_RAM 10000000'
      os.system(picard_call)


# Run
if __name__ == "__main__":
    main()
