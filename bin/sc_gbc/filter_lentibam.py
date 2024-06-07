#!/usr/bin/python


import os
import sys
import pysam
import gzip

#input_bam = '/Users/ieo6943/prova_pipeline/100/lentibam.bam'
#allowed_cbs_file = '/Users/ieo6943/prova_pipeline/100/barcodes.tsv.gz'
#output_bam = '/Users/ieo6943/prova_pipeline/100/filtered_lentibam.bam'
def filter_bam_by_cb(input_bam, output_bam, allowed_cbs_file):
    # Read allowed CB tags into a set from a gzipped TSV file
    allowed_cbs = set()
    with gzip.open(allowed_cbs_file, 'rt') as f: #rifare con pandas
        for line in f:
            cb_tag = line.strip().split()[0]  # Assuming the CB tags are in the first column
            allowed_cbs.add(cb_tag)

    # Open input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:
            for read in bam_in:# .fetch(until_eof=True): #togliere fetch vedere di ottimizzare
                if read.get_tag("CB") in allowed_cbs:
                    bam_out.write(read)

if __name__ == "__main__":
    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    allowed_cbs_file = sys.argv[3]
    
    filter_bam_by_cb(input_bam, output_bam, allowed_cbs_file)
