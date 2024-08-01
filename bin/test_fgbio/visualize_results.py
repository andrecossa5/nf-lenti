"""
Viz fgbio test results
"""

import os
import pandas as pd
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Read
path = '/Users/IEO5505/Desktop/example_mito/results/test_fgbio'
L = []
for x in os.listdir(path):
    L.append(pd.read_csv(os.path.join(path, x, 'stats.csv'), index_col=0).assign(run=x))

pd.concat(L).pivot(columns='run', values='value', index='metric')