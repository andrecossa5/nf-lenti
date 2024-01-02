#!/usr/bin/python

import sys
import pysam

# Args
bam_path = sys.argv[1]

# Read-format write
input_bam = pysam.AlignmentFile(bam_path, 'rb')
output_bam = pysam.AlignmentFile('mitobam_fixed.bam', 'wb', template=input_bam)

for r in input_bam:
    n1, n2, cb, ub = r.query_name.split('_')
    r.query_name = ' '.join([n1, n2])
    r.set_tag('CR', cb)
    r.set_tag('UR', ub)
    output_bam.write(r)

# Close
input_bam.close()
output_bam.close()