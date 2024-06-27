#!/usr/bin/python

import sys
import dnaio
import math

# Args
R1 = sys.argv[1]
R2 = sys.argv[2]
cpus = sys.argv[3]

def main(R1, R2, cpus):

    cpus = math.floor(float(cpus))
    with dnaio.open(R1, R2, mode='r', open_threads=cpus) as reader, \
        dnaio.open('assembled.fastq.gz', mode='w', open_threads=cpus) as writer:
        for r1, r2 in reader: 
            cbc, umi = r1.sequence[:16], r1.sequence[16:] 
            r2.name = '_'.join( r2.name.replace('/', '.').split(' ') + [cbc] + [umi] )
            writer.write(r2)


##


if __name__ == '__main__':
    main(R1, R2, cpus)

    