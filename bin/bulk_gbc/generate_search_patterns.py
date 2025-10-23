"""
generate_search_patterns.py

Generate search patterns for a given anchor sequence, replacing each position with a wildcard ('.') one at a time.
- Takes the anchor sequence as a command-line argument.
- Outputs all single-position wildcard patterns as a TSV file.
"""

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

