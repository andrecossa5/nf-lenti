#!/usr/bin/python

import sys
import pysam

# Args
bam_path = sys.argv[1]


def main():

    input_bam = pysam.AlignmentFile(bam_path, 'rb')
    output_bam = pysam.AlignmentFile('mitobam_fixed_tags.bam', 'wb', template=input_bam)

    for r in input_bam:
        n1, n2, cb, ub = r.query_name.split('_')
        r.query_name = ' '.join([n1, n2])
        r.set_tag('CB', cb)
        r.set_tag('UB', ub)
        output_bam.write(r)

    # Close
    input_bam.close()
    output_bam.close()

if __name__ == '__main__':
    main()