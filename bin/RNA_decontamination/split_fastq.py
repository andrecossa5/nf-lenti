import os
from dnaio import open as open_fastq
import dnaio
import pandas as pd
import argparse


##


def cbc_split_fastq(input_fastq1_path, input_fastq2_path, output_folder, barcodes):
    '''
    It creates "n" folders for all the different CBCs in the fastq (that are also presents in the barcodes list).
    In each folders write down the part of the fastq (read1 and read2 in separeted files) associated only to the folder CBC
    '''
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Dictionary to store output file handles based on CBC
    cbc_writers = {}
    
    with dnaio.open(input_fastq1_path, input_fastq2_path) as reader:
        #i=0
        for record1, record2 in reader:
            # Extract CBC from the read1 ID (assuming the CBC is at the end of the read ID)
            cbc = record1.sequence[:16]

            if cbc in barcodes:
                fastq_cbc = os.path.join(output_folder,f'{cbc}/')
                os.makedirs(fastq_cbc, exist_ok=True)
                
                # Check if output file writer for CBC exists, create if not
                if cbc not in cbc_writers:
                    #cbc_writers[cbc] = {
                    #    "read1": dnaio.open(os.path.join(fastq_cbc, f"reads_{cbc}_R1.fastq"), mode='w', fileformat='fastq'),
                    #    "read2": dnaio.open(os.path.join(fastq_cbc, f"reads_{cbc}_R2.fastq"), mode='w', fileformat='fastq')
                    #}

                    cbc_writers[cbc] = dnaio.open(os.path.join(fastq_cbc, f"reads_{cbc}_R1.fastq"),os.path.join(fastq_cbc, f"reads_{cbc}_R2.fastq"), mode='w', fileformat='fastq')

                # Write record1 and record2 to corresponding output files
                #cbc_writers[cbc]["read1"].write(record1)
                #cbc_writers[cbc]["read2"].write(record2)
                cbc_writers[cbc].write(record1, record2)
                #i+=1
                #if i==1000000: break
    
    # Close all output files
    for writers in cbc_writers.values():
        #for writer in writers.values():
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
    '--output_path', 
    type=str,
    default=None,
    help='Where to save the files'
)

my_parser.add_argument(
    '--read1_path', 
    type=str,
    default=None,
    help='Read1'
)

my_parser.add_argument(
    '--read2_path', 
    type=str,
    default=None,
    help='read2'
)

my_parser.add_argument(
    '--barcodes_path', 
    type=str,
    default=None,
    help='path of the barcodes'
)




#Path
#main = '/Users/ieo6943/Documents/data/AML_clones/'
#input_fastq1_path = os.path.join(main, 'head_S45545_854_enrichment_PCR_S2_L001_R1_001.fastq.gz')
#input_fastq2_path = os.path.join(main, 'head_S45545_854_enrichment_PCR_S2_L001_R2_001.fastq.gz')
#output_path = os.path.join(main,'cbc/')
#barcodes_path = os.path.join(main, 'barcodes.tsv.gz')

# Parse arguments
args = my_parser.parse_args()
input_fastq1_path = args.read1_path
input_fastq2_path = args.read2_path
barcodes_path = args.barcodes_path
output_path = args.output_path

##


solo_CBCs = pd.read_csv(
    barcodes_path, 
    header=None, index_col=0
)
solo_CBCs = list(solo_CBCs.index.unique())
barcodes = solo_CBCs
##


#
cbc_split_fastq(input_fastq1_path, input_fastq2_path, output_path, barcodes)

