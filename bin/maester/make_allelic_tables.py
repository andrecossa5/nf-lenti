#!/usr/bin/python

"""
Make allelic tables.
"""

import argparse
import numpy as np
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

#Input
my_parser.add_argument(
    '--consensus_bam', 
    type=str,
    default=None,
    help='path of the consensus bam'
)

my_parser.add_argument(
    '--cell', 
    type=str,
    default=None,
    help='Name of the cell'
)

my_parser.add_argument(
    '--min_base_qual', 
    type=int,
    default=30,
    help='min base quality of the cell'
)

my_parser.add_argument(
    '--min_alignment_quality', 
    type=int,
    default=60,
    help='min mapping quality of the consensus read to the reference.'
)

my_parser.add_argument(
    '--min_base_depth', 
    type=int,
    default=3,
    help='min base depth (i.e., n reads at that position in the UMI group).'
)

my_parser.add_argument(
    '--max_base_consensus_error', 
    type=float,
    default=.25,
    help='Max base consensus error (i.e., n reads not supporting the consesus call at that position in the UMI group).'
)


##


# Parse arguments
args = my_parser.parse_args()
cell = args.cell
consensus_bam = args.consensus_bam
base_quality_thr = args.min_base_qual
mapping_quality_thr = args.min_alignment_quality
base_depth_thr = args.min_base_depth
base_consensus_error_thr = args.max_base_consensus_error

# Default
maxBP = 16569                                                       # rRCS MT-genome length


##


def writeSparseMatrix(cell, mid, vec):
	with open(f'{cell}.{mid}.txt', 'w') as V:
		for i in range(0,int(maxBP)):
			if vec[i] > 0:
				V.write(f'{(i+1)},{cell},{str(vec[i])}\n')
                        
def writeSparseMatrix2(cell, mid, vec1, vec2):
    with open(f'{cell}.{mid}.txt', 'w') as V:
	    for i in range(0,int(maxBP)):
		    if vec1[i] > 0 or vec2[i] > 0:
			    V.write(f'{str(i+1)},{cell},{str(vec1[i])},{str(vec2[i])}\n')
                        
def writeSparseMatrix8(cell, mid, vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8):
    with open(f'{cell}.{mid}.txt', 'w') as V:
        for i in range(0,int(maxBP)):
            if(vec1[i] > 0 or vec5[i] > 0):
                V.write(f'{str(i+1)},{cell},{str(vec1[i])},{str(vec2[i])},{str(vec3[i])},{str(vec4[i])},{str(vec5[i])},{str(vec6[i])},{str(vec7[i])},{str(vec8[i])}\n')

##


def main():

    # initialize with a pseudo count to avoid dividing by zero
    n = int(maxBP)

    countsA_fw = [0.00000001] * n 
    countsC_fw = [0.00000001] * n 
    countsG_fw = [0.00000001] * n 
    countsT_fw = [0.00000001] * n 

    qualA_fw = [0.0] * n
    qualC_fw = [0.0] * n
    qualG_fw = [0.0] * n
    qualT_fw = [0.0] * n

    consensusA_fw = [0.0] * n 
    consensusC_fw = [0.0] * n 
    consensusG_fw = [0.0] * n 
    consensusT_fw = [0.0] * n 

    groupsizeA_fw = [0.0] * n 
    groupsizeC_fw = [0.0] * n 
    groupsizeG_fw = [0.0] * n 
    groupsizeT_fw = [0.0] * n 

    countsA_rev = [0.00000001] * n 
    countsC_rev = [0.00000001] * n 
    countsG_rev = [0.00000001] * n
    countsT_rev = [0.00000001] * n 

    qualA_rev = [0.0] * n
    qualC_rev = [0.0] * n
    qualG_rev = [0.0] * n
    qualT_rev = [0.0] * n

    consensusA_rev = [0.0] * n 
    consensusC_rev = [0.0] * n 
    consensusG_rev = [0.0] * n 
    consensusT_rev = [0.0] * n 

    groupsizeA_rev = [0.0] * n 
    groupsizeC_rev = [0.0] * n 
    groupsizeG_rev = [0.0] * n 
    groupsizeT_rev = [0.0] * n 


    # Read bam
    bam = pysam.AlignmentFile(consensus_bam, "rb", require_index=False)

    # Loop for each read and base
    n_umis_unfiltered = 0
    n_umis_filtered = 0
    read_groups_size = []

    for read in bam:
        
        seq = read.seq
        reverse = read.is_reverse
        base_qualities = read.query_qualities
        mapping_quality = read.mapping_quality
        base_dephts = read.get_tag('cd')
        base_n_discordant = read.get_tag('ce')

        if mapping_quality >= mapping_quality_thr:                  # UMI sequence test 

            for qpos, refpos in read.get_aligned_pairs(True):

                n_umis_unfiltered += 1
                q = base_qualities[qpos]
                d = base_dephts[qpos]
                n_discordant = base_n_discordant[qpos]
                r = n_discordant / (d+0.0000001)                    # Avoid dividing by zero
                base_test = ( qpos is not None ) and \
                            ( refpos is not None ) and \
                            ( q >= base_quality_thr ) and \
                            ( r <= base_consensus_error_thr ) and \
                            ( d >= base_depth_thr )

                if base_test:                                       # Individual base test 
                    
                    n_umis_filtered += 1
                    read_groups_size.append(d)

                    if seq[qpos] == "A":
                        if reverse:
                            qualA_fw[refpos] += q
                            countsA_fw[refpos] += 1
                            consensusA_fw[refpos] += 1-r 
                            groupsizeA_fw[refpos] = d
                        else:
                            qualA_rev[refpos] += q
                            countsA_rev[refpos] += 1
                            consensusA_rev[refpos] += 1-r 
                            groupsizeA_rev[refpos] = d
                    elif seq[qpos] == "C":
                        if reverse:
                            qualC_fw[refpos] += q
                            countsC_fw[refpos] += 1
                            consensusC_fw[refpos] += 1-r 
                            groupsizeC_fw[refpos] = d 
                        else:
                            qualC_rev[refpos] += q
                            countsC_rev[refpos] += 1
                            consensusC_rev[refpos] += 1-r
                            groupsizeC_rev[refpos] = d 
                    elif seq[qpos] == "G":
                        if reverse:
                            qualG_fw[refpos] += q
                            countsG_fw[refpos] += 1
                            consensusG_fw[refpos] += 1-r 
                            groupsizeG_fw[refpos] = d 
                        else:
                            qualG_rev[refpos] += q
                            countsG_rev[refpos] += 1
                            consensusG_rev[refpos] += 1-r 
                            groupsizeG_rev[refpos] = d 
                    elif seq[qpos] == "T":
                        if reverse:
                            qualT_fw[refpos] += q
                            countsT_fw[refpos] += 1
                            consensusT_fw[refpos] += 1-r 
                            groupsizeT_fw[refpos] = d 
                        else:
                            qualT_rev[refpos] += q
                            countsT_rev[refpos] += 1
                            consensusT_rev[refpos] += 1-r 
                            groupsizeT_rev[refpos] = d 

    # Get mean quality, mean consensus, and base counts lists
    meanQualA_fw = [round(x/y,1) for x, y in zip(qualA_fw, countsA_fw)]
    meanQualC_fw = [round(x/y,1) for x, y in zip(qualC_fw, countsC_fw)]
    meanQualG_fw = [round(x/y,1) for x, y in zip(qualG_fw, countsG_fw)]
    meanQualT_fw = [round(x/y,1) for x, y in zip(qualT_fw, countsT_fw)]

    meanConsA_fw = [round(x/y,1) for x, y in zip(consensusA_fw, countsA_fw)]
    meanConsC_fw = [round(x/y,1) for x, y in zip(consensusC_fw, countsC_fw)]
    meanConsG_fw = [round(x/y,1) for x, y in zip(consensusG_fw, countsG_fw)]
    meanConsT_fw = [round(x/y,1) for x, y in zip(consensusT_fw, countsT_fw)]

    meanGSizeA_fw = [round(x/y,1) for x, y in zip(groupsizeA_fw, countsA_fw)]
    meanGSizeC_fw = [round(x/y,1) for x, y in zip(groupsizeC_fw, countsC_fw)]
    meanGSizeG_fw = [round(x/y,1) for x, y in zip(groupsizeG_fw, countsG_fw)]
    meanGSizeT_fw = [round(x/y,1) for x, y in zip(groupsizeT_fw, countsT_fw)]

    countsA_fw = [ int(round(elem)) for elem in countsA_fw ]
    countsC_fw = [ int(round(elem)) for elem in countsC_fw ]
    countsG_fw = [ int(round(elem)) for elem in countsG_fw ]
    countsT_fw = [ int(round(elem)) for elem in countsT_fw ]

    meanQualA_rev = [round(x/y,1) for x, y in zip(qualA_rev, countsA_rev)]
    meanQualC_rev = [round(x/y,1) for x, y in zip(qualC_rev, countsC_rev)]
    meanQualG_rev = [round(x/y,1) for x, y in zip(qualG_rev, countsG_rev)]
    meanQualT_rev = [round(x/y,1) for x, y in zip(qualT_rev, countsT_rev)]

    meanConsA_rev = [round(x/y,1) for x, y in zip(consensusA_rev, countsA_rev)]
    meanConsC_rev = [round(x/y,1) for x, y in zip(consensusC_rev, countsC_rev)]
    meanConsG_rev = [round(x/y,1) for x, y in zip(consensusG_rev, countsG_rev)]
    meanConsT_rev = [round(x/y,1) for x, y in zip(consensusT_rev, countsT_rev)]

    meanGSizeA_rev = [round(x/y,1) for x, y in zip(groupsizeA_rev, countsA_rev)]
    meanGSizeC_rev = [round(x/y,1) for x, y in zip(groupsizeC_rev, countsC_rev)]
    meanGSizeG_rev = [round(x/y,1) for x, y in zip(groupsizeG_rev, countsG_rev)]
    meanGSizeT_rev = [round(x/y,1) for x, y in zip(groupsizeT_rev, countsT_rev)]

    countsA_rev = [ int(round(elem)) for elem in countsA_rev ]
    countsC_rev = [ int(round(elem)) for elem in countsC_rev ]
    countsG_rev = [ int(round(elem)) for elem in countsG_rev ]
    countsT_rev = [ int(round(elem)) for elem in countsT_rev ]

    # Write as sparse
    writeSparseMatrix8(cell, "A", countsA_fw, meanQualA_fw, meanConsA_fw, meanGSizeA_fw, countsA_rev, meanQualA_rev, meanConsA_rev, meanGSizeA_rev)
    writeSparseMatrix8(cell, "C", countsC_fw, meanQualC_fw, meanConsC_fw, meanGSizeC_fw, countsC_rev, meanQualC_rev, meanConsC_rev, meanGSizeC_rev)
    writeSparseMatrix8(cell, "G", countsG_fw, meanQualG_fw, meanConsG_fw, meanGSizeG_fw, countsG_rev, meanQualG_rev, meanConsG_rev, meanGSizeG_rev)
    writeSparseMatrix8(cell, "T", countsT_fw, meanQualT_fw, meanConsT_fw, meanGSizeT_fw, countsT_rev, meanQualT_rev, meanConsT_rev, meanGSizeT_rev)

    zipped_list = zip(list(countsA_fw),list(countsC_fw),list(countsG_fw),list(countsT_fw), list(countsA_rev),list(countsC_rev),list(countsG_rev),list(countsT_rev))
    sums = [sum(item) for item in zipped_list]
    writeSparseMatrix(cell, 'coverage', sums)

    # Coverage
    with open(f'{cell}.coverage.txt', 'r') as coverage:
        depth = 0
        for row in coverage:
            s = row.split(',')
            depth += int(s[2].strip())

    # Depth
    with open(f'{cell}.depth.txt', 'w') as d:
    	d.write(f'{cell},{round(float(depth)/float(maxBP),2)}\n')
         
    # Added metrics
    median_base_umi_group_size = np.median(read_groups_size)
    with open(f'{cell}.median_base_umi_group_size.txt', 'w') as f:
        f.write(f'{cell},{median_base_umi_group_size}\n')
         
    with open(f'{cell}.n_umis_unfiltered.txt', 'w') as f:
        f.write(f'{cell},{n_umis_unfiltered}\n')
         
    with open(f'{cell}.n_umis_filtered.txt', 'w') as f:
        f.write(f'{cell},{n_umis_filtered}\n')
         

##


# Run 
if __name__ == '__main__':
    main()
