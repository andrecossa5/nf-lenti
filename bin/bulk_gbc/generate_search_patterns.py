#!/usr/bin/python

import os
import sys
import pandas as pd

# Anchor sequence
sequence = sys.argv[1]

# Patterns
L = []
for i in range(len(sequence)):
    s = list(sequence)
    s[i] = '.'
    L.append(''.join(s+['.*']))

# Save
(
    pd.DataFrame(L)
    .to_csv(
        os.path.join(os.getcwd(), 'search_patterns.tsv'),
        header=False, index=False, sep='\t'
    )
)


##

