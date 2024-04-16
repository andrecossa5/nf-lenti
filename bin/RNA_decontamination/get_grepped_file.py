#!/usr/bin/python

import dnaio
import os
import argparse
import gzip


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
output_read = f"/hpcnfs/home/ieo6943/results/{type}/{sample}/grepped.txt"

print('path_read1=', input_read_1 )
print('path_read2=', input_read_2)

print('read2')
#main_path = f'/Users/ieo6943/Documents/data/AML_clones/'
#input_read_1 = os.path.join(main_path, 'head_S45545_854_enrichment_PCR_S2_L001_R1_001.fastq.gz')
#input_read_2 = os.path.join(main_path,'head_S45545_854_enrichment_PCR_S2_L001_R2_001.fastq.gz')
#output_read = os.path.join(main_path,'enrich_grepped.txt')
#    
##
#
#print('input_read=', input_read_1)
#print(output_read)
#
with gzip.open(input_read_2, 'rb') as file:
    # Leggi le prime 8 righe
    for _ in range(8):
        riga = file.readline()
        print(riga.decode().strip())
#
#with open(output_read, 'r') as file:
#    # Leggi le prime dieci righe
#    for i in range(10):
#        riga = file.readline()
#        print(riga.strip())
        
with dnaio.open(input_read_1, input_read_2) as reader, open(output_read,'w') as writer:
    writer.write(f"Name\tcbc\tumi\tfeature\n")
    for record_1, record_2 in reader:
        writer.write(f"{record_1.name}\t{record_1.sequence[:16]}\t{record_1.sequence[16:]}\t{record_2.sequence}\n")