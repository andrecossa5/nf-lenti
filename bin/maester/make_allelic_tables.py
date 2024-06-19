#!/usr/bin/python

"""
Make allelic tables.
"""
 

##


# Libraries
import os
import sys
import argparse
import pysam


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='make_allelic_tables',
    description=
    """
    Create cell-specific allelic tables.
    """
)

# Input
# my_parser.add_argument(
#     '--sample', 
#     type=str,
#     default=None,
#     help='Sample name. Default: None.'
# )




##


# Parse arguments
args = my_parser.parse_args()
bamfile = args.input_bam
cell = args.cell
base_qual = args.min_base_qual
alignment_quality = args.min_alignment_quality


##
outpre = cell
outputdepth = f'{cell}.depth.txt'
sample = cell
maxBP = 16569
base_qual = 30
#cell = 'TCCATCGAGAATAACC'
#bamfile = '/Users/ieo6943/Desktop/bam_file_maester/TCCATCGAGAATAACC_consensus_filtered_mapped.bam'
#outpre = f'/Users/ieo6943/Documents/Guido/scratch/Lareu/{cell}'

#outputdepth = f'/Users/ieo6943/Documents/Guido/scratch/Lareu/{cell}.depth.txt'

#mito_genome = '/Users/ieo6943/Documents/Guido/scratch/ref/new_genome_masked.fa'
#alignment_quality = float(30)


##


# Export Functions
def writeSparseMatrix(mid, vec):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec[i])+"\n")
def writeSparseMatrix2(mid, vec1, vec2):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0 or vec2[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+"\n")
def writeSparseMatrix4(mid, vec1, vec2, vec3, vec4):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0 or vec3[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+","+str(vec3[i])+","+str(vec4[i])+"\n")

def main():
    n = int(maxBP)

    # initialize with a pseudo count to avoid dividing by zero
    countsA_fw = [0.00000001] * n 
    countsC_fw = [0.00000001] * n 
    countsG_fw = [0.00000001] * n 
    countsT_fw = [0.00000001] * n 

    qualA_fw = [0.0] * n
    qualC_fw = [0.0] * n
    qualG_fw = [0.0] * n
    qualT_fw = [0.0] * n

    countsA_rev = [0.00000001] * n 
    countsC_rev = [0.00000001] * n 
    countsG_rev = [0.00000001] * n 
    countsT_rev = [0.00000001] * n 

    qualA_rev = [0.0] * n
    qualC_rev = [0.0] * n
    qualG_rev = [0.0] * n
    qualT_rev = [0.0] * n


    bam2 = pysam.AlignmentFile(bamfile, "rb")
    for read in bam2:
        seq = read.seq
        reverse = read.is_reverse
        quality = read.query_qualities
        align_qual_read = read.mapping_quality
        for qpos, refpos in read.get_aligned_pairs(True):
            if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
                if(seq[qpos] == "A" and quality[qpos] > base_qual):
                    if(reverse):
                        qualA_fw[refpos] += quality[qpos]
                        countsA_fw[refpos] += 1
                    else:
                        qualA_rev[refpos] += quality[qpos]
                        countsA_rev[refpos] += 1
                elif(seq[qpos] == "C" and quality[qpos] > base_qual):
                    if(reverse):
                        qualC_fw[refpos] += quality[qpos]
                        countsC_fw[refpos] += 1
                    else:
                        qualC_rev[refpos] += quality[qpos]
                        countsC_rev[refpos] += 1
                elif(seq[qpos] == "G" and quality[qpos] > base_qual):
					
                    if(reverse):
                        qualG_fw[refpos] += quality[qpos]
                        countsG_fw[refpos] += 1
                    else:
                        qualG_rev[refpos] += quality[qpos]
                        countsG_rev[refpos] += 1
                elif(seq[qpos] == "T" and quality[qpos] > base_qual):
                    if(reverse):
                        qualT_fw[refpos] += quality[qpos]
                        countsT_fw[refpos] += 1
                    else:
                        qualT_rev[refpos] += quality[qpos]
                        countsT_rev[refpos] += 1
    
    meanQualA_fw = [round(x/y,1) for x, y in zip(qualA_fw, countsA_fw)]
    meanQualC_fw = [round(x/y,1) for x, y in zip(qualC_fw, countsC_fw)]
    meanQualG_fw = [round(x/y,1) for x, y in zip(qualG_fw, countsG_fw)]
    meanQualT_fw = [round(x/y,1) for x, y in zip(qualT_fw, countsT_fw)]

    countsA_fw = [ int(round(elem)) for elem in countsA_fw ]
    countsC_fw = [ int(round(elem)) for elem in countsC_fw ]
    countsG_fw = [ int(round(elem)) for elem in countsG_fw ]
    countsT_fw = [ int(round(elem)) for elem in countsT_fw ]

    meanQualA_rev = [round(x/y,1) for x, y in zip(qualA_rev, countsA_rev)]
    meanQualC_rev = [round(x/y,1) for x, y in zip(qualC_rev, countsC_rev)]
    meanQualG_rev = [round(x/y,1) for x, y in zip(qualG_rev, countsG_rev)]
    meanQualT_rev = [round(x/y,1) for x, y in zip(qualT_rev, countsT_rev)]

    countsA_rev = [ int(round(elem)) for elem in countsA_rev ]
    countsC_rev = [ int(round(elem)) for elem in countsC_rev ]
    countsG_rev = [ int(round(elem)) for elem in countsG_rev ]
    countsT_rev = [ int(round(elem)) for elem in countsT_rev ]

    # Allele Counts

    writeSparseMatrix4("A", countsA_fw, meanQualA_fw, countsA_rev, meanQualA_rev)
    writeSparseMatrix4("C", countsC_fw, meanQualC_fw, countsC_rev, meanQualC_rev)
    writeSparseMatrix4("G", countsG_fw, meanQualG_fw, countsG_rev, meanQualG_rev)
    writeSparseMatrix4("T", countsT_fw, meanQualT_fw, countsT_rev, meanQualT_rev)

    zipped_list = zip(list(countsA_fw),list(countsC_fw),list(countsG_fw),list(countsT_fw), list(countsA_rev),list(countsC_rev),list(countsG_rev),list(countsT_rev))
    sums = [sum(item) for item in zipped_list]
    writeSparseMatrix("coverage", sums)

#Get depth from the coverage sparse matrix
    with open(outpre + ".coverage.txt", 'r') as coverage:
        depth = 0
        for row in coverage:
            s = row.split(",")
            depth += int(s[2].strip())
    with open(outputdepth, 'w') as d:
    	d.write(sample + "\t" + str(round(float(depth)/float(maxBP),2)) + "\n")


##


# Run
if __name__ == '__main__':
    main()


