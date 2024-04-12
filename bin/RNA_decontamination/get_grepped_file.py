#!/usr/bin/python

import dnaio
import os

# File path
sample = 'a'
type = 'b'
name = 'c'

#main_path = f'/hpcnfs/scratch/PGP/MI_TO_benchmark/preprocessing/pp/data/maester/AML_clones'
main_path = f'/hpcnfs/scratch/PGP/MI_TO_benchmark/preprocessing/pp/data/sc_tenx/AML_clones'
input_read_1 = os.path.join(main_path, "S45555_854_S11_L001_R1_001.fastq.gz")
input_read_2 = os.path.join(main_path,"S45555_854_S11_L001_R2_001.fastq.gz")
#output_read = "/hpcnfs/home/ieo6943/scracth/grepped.txt"

#main_path = f'/hpcnfs/scratch/PGP/MI_TO_benchmark/preprocessing/pp/data/{type}/{sample}'
#input_read_1 = os.path.join(main_path, f"{name}_R1_001.fastq.gz")
#input_read_2 = os.path.join(main_path,f"{name}_R2_001.fastq.gz")
output_read = "/hpcnfs/home/ieo6943/scracth/grepped_sc_tenex.txt"
    
#
        
with dnaio.open(input_read_1, input_read_2) as reader, open(output_read,'w') as writer:
    writer.write(f"Name\tcbc\tumi\tfeature\n")
    for record_1, record_2 in reader:
        writer.write(f"{record_1.name}\t{record_1.sequence[:16]}\t{record_1.sequence[16:]}\t{record_2.sequence}\n")