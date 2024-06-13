#!/usr/bin/python

"""
Filter a bam file for records having the needed CBs.
"""

import sys
import pysam
import pandas as pd


##


# Args
allowed_cbs_file = sys.argv[1]

##

def main():
    allowed_cbs = pd.read_csv(allowed_cbs_file, sep='\t').iloc[:,0].values.tolist()
    with pysam.AlignmentFile('merged_mitobam.bam', "rb") as bam_in:
        with pysam.AlignmentFile('filtered_mitobam.bam', "wb", header=bam_in.header) as bam_out:
            for read in bam_in:
                if read.get_tag("CB") in allowed_cbs:
                    bam_out.write(read)

# Run
if __name__ == "__main__":
    main()
