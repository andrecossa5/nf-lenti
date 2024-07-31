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
    '--min_base_consensus_error', 
    type=float,
    default=.25,
    help='min base consensus error (i.e., n reads not supporting the consesus call at that position in the UMI group).'
)

my_parser.add_argument(
    '--filtering', 
    type=str,
    default='custom',
    help='Type of filtering: default, custom.'
)


##


# Parse arguments
args = my_parser.parse_args()
cell = args.cell
consensus_bam = args.consensus_bam
base_quality_thr = args.min_base_qual
mapping_quality_thr = args.min_alignment_quality
base_depth_thr = args.min_base_depth
base_consensus_error_thr = args.min_base_consensus_error
filtering = args.filtering

# Test
# import os
# path_ = '/Users/IEO5505/Desktop/example_mito/scratch'
# os.chdir(path_)
# cell = 'AAA'
# consensus_bam = os.path.join(path_, 'cons_mapped.bam')
# base_quality_thr = 30
# mapping_quality_thr = 30
# base_depth_thr = 3
# base_consensus_error_thr = .2

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
                        
def writeSparseMatrix4(cell, mid, vec1, vec2, vec3, vec4):
    with open(f'{cell}.{mid}.txt', 'w') as V:
        for i in range(0,int(maxBP)):
            if(vec1[i] > 0 or vec3[i] > 0):
                V.write(f'{str(i+1)},{cell},{str(vec1[i])},{str(vec2[i])},{str(vec3[i])},{str(vec4[i])}\n')


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

    countsA_rev = [0.00000001] * n 
    countsC_rev = [0.00000001] * n 
    countsG_rev = [0.00000001] * n
    countsT_rev = [0.00000001] * n 

    qualA_rev = [0.0] * n
    qualC_rev = [0.0] * n
    qualG_rev = [0.0] * n
    qualT_rev = [0.0] * n


    # Read bam
    bam = pysam.AlignmentFile(consensus_bam, "rb", require_index=False)

    # Loop for each read and base
    n_umis_unfiltered = 0
    n_umis_filtered = 0
    n_reads_unfiltered = 0
    n_reads_filtered = 0
    read_groups_size = []
    base_consensus = []
    median_read_base_qualities = []

    for read in bam:
        
        n_reads_unfiltered += 1

        seq = read.seq
        reverse = read.is_reverse
        base_qualities = read.query_qualities
        mapping_quality = read.mapping_quality
        base_dephts = read.get_tag('cd')
        base_n_discordant = read.get_tag('ce')

        if mapping_quality >= mapping_quality_thr:                           # Entire sequence test 
            
            n_reads_filtered += 1
            median_read_base_qualities.append(np.median(base_qualities))

            for qpos, refpos in read.get_aligned_pairs(True):

                n_umis_unfiltered += 1
                q = base_qualities[qpos]
                d = base_dephts[qpos]
                n_discordant = base_n_discordant[qpos]
                r = n_discordant / (d+0.0000001)                             # Avoid dividing by zero

                if filtering == 'custom':
                    base_test = ( qpos is not None ) and \
                                ( refpos is not None ) and \
                                ( q >= base_quality_thr ) and \
                                ( r <= base_consensus_error_thr ) and \
                                ( d >= base_depth_thr )
                else:
                    base_test = ( qpos is not None ) and \
                                ( refpos is not None ) and \
                                ( q >= base_quality_thr )

                if base_test:                                                # Individual base test 
                    
                    n_umis_filtered += 1
                    read_groups_size.append(d)
                    base_consensus.append(r)

                    if seq[qpos] == "A":
                        if reverse:
                            qualA_fw[refpos] += q
                            countsA_fw[refpos] += 1
                        else:
                            qualA_rev[refpos] += q
                            countsA_rev[refpos] += 1
                    elif seq[qpos] == "C":
                        if reverse:
                            qualC_fw[refpos] += q
                            countsC_fw[refpos] += 1
                        else:
                            qualC_rev[refpos] += q
                            countsC_rev[refpos] += 1
                    elif seq[qpos] == "G":
                        if reverse:
                            qualG_fw[refpos] += q
                            countsG_fw[refpos] += 1
                        else:
                            qualG_rev[refpos] += q
                            countsG_rev[refpos] += 1
                    elif seq[qpos] == "T":
                        if reverse:
                            qualT_fw[refpos] += q
                            countsT_fw[refpos] += 1
                        else:
                            qualT_rev[refpos] += q
                            countsT_rev[refpos] += 1


    # Get Qual and base counts lists
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

    # Write as sparse
    writeSparseMatrix4(cell, "A", countsA_fw, meanQualA_fw, countsA_rev, meanQualA_rev)
    writeSparseMatrix4(cell, "C", countsC_fw, meanQualC_fw, countsC_rev, meanQualC_rev)
    writeSparseMatrix4(cell, "G", countsG_fw, meanQualG_fw, countsG_rev, meanQualG_rev)
    writeSparseMatrix4(cell, "T", countsT_fw, meanQualT_fw, countsT_rev, meanQualT_rev)

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
    with open(f'{cell}.median_filtered_base_umi_group_size.txt', 'w') as f:
        f.write(f'{cell},{median_base_umi_group_size}\n')

    median_base_consensus = np.median(base_consensus)
    with open(f'{cell}.median_filtered_base_consensus_error.txt', 'w') as f:
        f.write(f'{cell},{median_base_consensus}\n')
    
    median_read_quality = np.median(median_read_base_qualities)
    with open(f'{cell}.median_filtered_read_quality.txt', 'w') as f:
        f.write(f'{cell},{median_read_quality}\n')
         
    with open(f'{cell}.n_umis_unfiltered.txt', 'w') as f:
        f.write(f'{cell},{n_umis_unfiltered}\n')
         
    with open(f'{cell}.n_umis_filtered.txt', 'w') as f:
        f.write(f'{cell},{n_umis_filtered}\n')

    with open(f'{cell}.n_reads_unfiltered.txt', 'w') as f:
        f.write(f'{cell},{n_reads_unfiltered}\n')
    
    with open(f'{cell}.n_reads_filtered.txt', 'w') as f:
        f.write(f'{cell},{n_reads_filtered}\n')
         

##


# Run 
if __name__ == '__main__':
    main()


