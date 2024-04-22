#!/usr/bin/python

import dnaio
import os
import argparse
import pandas as pd


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='get_grepped_file',
    description=
    """
    Get grepped file.
    """
)

my_parser.add_argument(
    '--type', 
    type=str,
    default=None,
    help='type of experiment: lentiviral, GEX, mitocondrial'
)

my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Datasample: AML_clones, MDA_clones, ...'
)

my_parser.add_argument(
    '--read1', 
    type=str,
    default=None,
    help='Read1'
)

my_parser.add_argument(
    '--read2', 
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


# Parse arguments
args = my_parser.parse_args()
type =args.type
sample = args.sample
read1 = args.read1
read2 = args.read2

##


#Arg
main_path = f'/hpcnfs/scratch/PGP/MI_TO_benchmark/preprocessing/pp/data/{type}/{sample}'
input_read_1 = os.path.join(main_path, read1)
input_read_2 = os.path.join(main_path,read2)
path_barcodes = args.barcodes_path
path_barcodes = os.path.join(path_barcodes, 'barcodes.tsv.gz')
output_read = f"/hpcnfs/home/ieo6943/results/{type}/{sample}/grepped.txt"


#main_path = '/Users/ieo6943/Documents/data/AML_clones/'
#input_read_1 = os.path.join(main_path, 'head_S45545_854_enrichment_PCR_S2_L001_R1_001.fastq.gz')
#input_read_2 = os.path.join(main_path,'head_S45545_854_enrichment_PCR_S2_L001_R2_001.fastq.gz')
#path_barcodes = os.path.join(main_path, 'barcodes.tsv.gz')
#output_read = os.path.join(main_path,'grepped_scratch.txt.gz')


##


solo_CBCs = pd.read_csv(
    path_barcodes, 
    header=None, index_col=0
)
solo_CBCs = list(solo_CBCs.index.unique())
        
with dnaio.open(input_read_1, input_read_2) as reader, open(output_read,'w') as writer:
    writer.write(f"Name\tcbc\tumi\tfeature\n")
    for record_1, record_2 in reader:
        cbc = record_1.sequence[:16]
        umi = record_1.sequence[16:]
        feature = record_2.sequence
        if  cbc in solo_CBCs:
            writer.write(f"{record_1.name}\t{cbc}\t{umi}\t{feature}\n")