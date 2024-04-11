#!/usr/bin/python

import dnaio
import os

# File path
main_path = '/hpcnfs/scratch/PGP/MI_TO_benchmark/preprocessing/pp/data/maester/AML_clones'
input_read_1 = os.path.join(main_path, "S46393_AML_854_MAESTER_S2_L001_R1_001.fastq.gz")
input_read_2 = os.path.join(main_path,"S46393_AML_854_MAESTER_S2_L001_R2_001.fastq.gz")
output_read = "/hpcnfs/home/ieo6943/scracth/grepped.txt"
    
#
        
with dnaio.open(input_read_1, input_read_2) as reader, open(output_read,'w') as writer:
    writer.write(f"Name\tcbc\tumi\tfeature\n")
    for record_1, record_2 in reader:
        writer.write(f"{record_1.name}\t{record_1.sequence[:16]}\t{record_1.sequence[16:]}\t{record_2.sequence}\n")