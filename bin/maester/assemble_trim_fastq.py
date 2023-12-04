#!/usr/bin/python

import sys
import dnaio
import math

# Args
R1 = sys.argv[1]
R2 = sys.argv[2]
cpus = sys.argv[3]

# Set cores
cpus = math.floor(float(cpus))

# Read-format write
with dnaio.open(R1, R2, mode='r', open_threads=cpus) as reader, \
    dnaio.open(f'assembled_fastq.gz', mode='w', open_threads=cpus) as writer:

    for r1, r2 in reader: 
        cbc, umi = r1.sequence[:16], r1.sequence[16:] 
        r2.name = '_'.join( r2.name.replace('/', '.').split(' ') + [cbc] + [umi] )
        writer.write(r2)