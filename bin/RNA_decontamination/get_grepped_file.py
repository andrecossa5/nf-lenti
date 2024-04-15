#!/usr/bin/python

import dnaio
import os
import argparse


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
output_read = f"/hpcnfs/home/ieo6943/results/{type}/grepped.txt"
    
#

print('input_read=', input_read_1)
        
with dnaio.open(input_read_1, input_read_2) as reader, open(output_read,'w') as writer:
    writer.write(f"Name\tcbc\tumi\tfeature\n")
    for record_1, record_2 in reader:
        writer.write(f"{record_1.name}\t{record_1.sequence[:16]}\t{record_1.sequence[16:]}\t{record_2.sequence}\n")