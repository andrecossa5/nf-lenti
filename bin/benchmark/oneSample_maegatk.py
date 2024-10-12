#!/usr/bin/python

import os
import sys
import pysam


##


# Parse input bam, cell name and reference .fa file
inputbam = sys.argv[1]
sample = sys.argv[2]
fasta_file = sys.argv[3]

# Define default maegatk parameters within the script

# mito_genome = config["mito_chr"]
# mito_length = str(config["mito_length"])
# fasta_file = config["fasta_file"]

# Get this info from fasta_file
mito_genome = 'chrM' 
mito_length = str(16569)

# Hard-code all the other variables according to https://github.com/caleblareau/maegatk and the MAESTER paper (Miller et el., 2022)
# umi_barcode = config["umi_barcode"]
umi_barcode = 'UB'
# base_qual = str(config["base_qual"])
base_qual = str(0)
# alignment_quality = config["alignment_quality"]
alignment_quality = str(0)
# NHmax = config["NHmax"]
NHmax = str(2) 
# NMmax = config["NMmax"]
NMmax = str(15)
# min_reads = str(config["min_reads"])
min_reads = str(3)
# max_javamem  = config["max_javamem"]

# Software paths
python = "python"
fgbio = 'fgbio'

# Script locations: oneSample_maegatk.py and other script are placed in the same folder
script_dir = os.path.dirname(os.path.abspath(__file__))
filtclip_py = os.path.join(script_dir, "filterClipBam.py")
sumstatsBP_py = os.path.join(script_dir, "sumstatsBP.py")

# Prepare intermediate filenames: skip of the original script .replace synthax needed by Snakemake 
filtlog = 'rmlog.log'
filtlog = 'filter.log'
temp_bam0 = 'temp0.bam'
temp_bam1 = 'temp1.bam' 
temp_sam15 = 'sam15.sam'
temp_bam15 = 'bam15.bam'
temp_bam2 = 'bam2.bam'
temp_fastq = 'temp0.fastq'
outputbam = 'consensus.bam' # The last one, actually used for pileup by "sumstatsBP.py" 
prefixSM = sample
outputdepth = sample + ".depth.txt"


##


# From here on, the code is exactly as it is in https://github.com/caleblareau/maegatk

# 1) Filter bam files
proper_paired = "False"
pycall = " ".join([python, filtclip_py, inputbam, filtlog, mito_genome, proper_paired, NHmax, NMmax]) + " > " + temp_bam0
os.system(pycall)

# 2) Sort the filtered bam file
fgcallone =  fgbio + " GroupReadsByUmi -s Identity -e 0 -i " + temp_bam0 + " -o " + temp_bam1 + " -t " + umi_barcode
os.system('echo "'+fgcallone+'"')
os.system(fgcallone)

# 2.5) Modify the UB tag
samtoolscall = 'samtools view -H ' + temp_bam1 + '> ' + temp_sam15 + '; samtools view ' + temp_bam1 + '| awk \'OFS="\t" {$13=$13""$4; print $0}\' >> ' + temp_sam15 + '; samtools view -b ' + temp_sam15 + '> ' + temp_bam15
os.system('echo "'+samtoolscall+'"')
os.system(samtoolscall)

# 3) Call consensus reads
fgcalltwo = fgbio + " CallMolecularConsensusReads -t "+umi_barcode+" -i "+temp_bam15+" -o " + temp_bam2 +" -M " + min_reads
os.system(fgcalltwo)
print(fgcalltwo)

# 4) Convert consensus bam to fastq
# bedtools_call = "bedtools bamtofastq -i "+ temp_bam2 +" -fq " + temp_fastq # Bedtools stopped working for some reason, replacing it with samtools fastq
samtoolscall2 = 'samtools fastq -T cM ' + temp_bam2 + " | sed 's/\tcM:i:/_/g' > " + temp_fastq
os.system(samtoolscall2)

# 5) Remap + sort bam files
bwa_call = "bwa mem " + fasta_file + " " + temp_fastq + " | samtools sort -o "+ outputbam +" -"
os.system(bwa_call)
pysam.index(outputbam)

# 6) Get allele counts per sample / base pair and per-base quality scores
alleleCountcall = " ".join([python, sumstatsBP_py, outputbam, prefixSM, mito_genome, mito_length, base_qual, sample, fasta_file, alignment_quality])
os.system(alleleCountcall)

# 7) Get depth from the coverage sparse matrix
with open(prefixSM + ".coverage.txt", 'r') as coverage:
	depth = 0
	for row in coverage:
		s = row.split(",")
		depth += int(s[2].strip())

with open(outputdepth, 'w') as d:
	d.write(sample + "\t" + str(round(float(depth)/float(mito_length),2)) + "\n")


## 