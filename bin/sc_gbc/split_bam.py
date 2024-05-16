import os
import pandas as pd
import argparse
import pysam


##


def cbc_split_bam(bam_file, output_folder):
    '''
    It creates "n" folders for all the different CBCs in the fastq (that are also presents in the barcodes list).
    In each folders write down the part of the fastq (read1 and read2 in separeted files) associated only to the folder CBC
    '''
    os.makedirs(output_folder, exist_ok=True)
    cbc_tag = "CR:Z:"
    cbc_writers = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam_input:
        for alignment in bam_input:
            cbc = alignment.get_tag(cbc_tag)
            bam_cbc_path = os.path.join(output_folder,f'{cbc}/')
            os.makedirs(bam_cbc_path, exist_ok=True)
            if cbc not in cbc_writers:
                cbc_writers[cbc] = pysam.AlignmentFile(os.path.join(bam_cbc_path, f"{cbc}.bam"),'wb',header=bam_input.header )
            cbc_writers[cbc].write(alignment)
    for writers in cbc_writers.values():
        writers.close()


##


#ARGpars
# Create the parser
my_parser = argparse.ArgumentParser(
    prog='split_fastq',
    description=
    """
    It splits fastq in the different cbc fastq files 
    """
)

my_parser.add_argument(
    '--output_folder', 
    type=str,
    default=None,
    help='Where to save the files'
)

my_parser.add_argument(
    '--bam_input', 
    type=str,
    default=None,
    help='Read1'
)


# Parse arguments
args = my_parser.parse_args()
bam_input = args.bam_input
output_folder = args.output_folder

##



#
cbc_split_bam(bam_input, output_folder)

