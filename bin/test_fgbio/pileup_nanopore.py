#!/usr/bin/python

"""
Nanopore consensus.
"""

import os
import sys
import numpy as np
import pandas as pd
import pysam
# from sklearn.metrics import pairwise_distances


##


# Utils
def find_query_base_at_position(read, target_reference_name, target_position):

    if read.reference_name != target_reference_name:
        # print('Different chrom')
        return 'Different chrom'  # The read is not aligned to the desired chromosome

    ref_pos = read.reference_start
    read_pos = 0  # This will track the position in the read's query sequence

    # Iterate through the CIGAR operations
    for cigar_tuple in read.cigartuples:
        cigar_op = cigar_tuple[0]
        cigar_len = cigar_tuple[1]

        if cigar_op == 0:  # Match or mismatch (M)
            if ref_pos <= target_position < ref_pos + cigar_len:
                # The target position falls within this match/mismatch block
                read_index = read_pos + (target_position - ref_pos)
                # print(read.query_sequence[read_index])
                return read.query_sequence[read_index]
            
            ref_pos += cigar_len
            read_pos += cigar_len

        elif cigar_op == 1:  # Insertion to the reference (I)
            read_pos += cigar_len
        elif cigar_op == 2:  # Deletion from the reference (D)
            ref_pos += cigar_len
        elif cigar_op == 3:  # Skipped region from the reference (N)
            ref_pos += cigar_len
        elif cigar_op == 4:  # Soft clip (S)
            read_pos += cigar_len
        elif cigar_op == 5:  # Hard clip (H)
            continue  # Hard clipping doesn't consume query or reference
        elif cigar_op == 6:  # Padding (P)
            continue  # Padding doesn't consume query or reference
        else:
            raise ValueError(f"Unhandled CIGAR operation: {cigar_op}")
        
    # print('Pos uncovered')
    return 'Pos uncovered'  # The target position was not covered by this read


##


def filter_UMIs(path_bam, positions=None, mapping_quality_thr=60, median_base_quality_thr=20):
    """
    Get .bam stats.
    """

    # Read bam
    bam = pysam.AlignmentFile(path_bam, "rb", require_index=False)
    chromosomes = positions.iloc[:,0].unique()

    # Loop for each read and base
    n_umis_unfiltered = 0
    n_reads_unfiltered = 0
    n_reads_filtered = 0
    median_read_base_qualities = []
    median_mapping_qualities = []
    read_lengths = []
    
    # Here we go
    UMIs = {}
    for read in bam:

        if read.seq is None:    # BAH
            continue
 
        n_reads_unfiltered += 1
        n_umis_unfiltered += len(read.seq)
        base_qualities = read.query_qualities
        mapping_quality = read.mapping_quality
        BQ = np.median(base_qualities)
        median_read_base_qualities.append(BQ)
        median_mapping_qualities.append(mapping_quality)
        read_lengths.append(len(read.seq))
        # assert read.reference_name in chromosomes
        seq_test = mapping_quality >= mapping_quality_thr and BQ >= median_base_quality_thr and read.reference_name in chromosomes
    
        if seq_test:                                    # Entire sequence test 
            
            n_reads_filtered += 1
            UB = read.get_tag('UB')
            if UB in UMIs:
                UMIs[UB].append(read)
            else:
                UMIs[UB] = [read]
        
    # Record .bam stats
    UMI_group_sizes = pd.Series(np.array([ len(x) for x in UMIs.values() ]), index=UMIs.keys())
    d = {}
    d['n_reads_filtered'] = n_reads_filtered
    d['n_reads_unfiltered'] = n_reads_unfiltered
    d['read_lengths'] = np.median(read_lengths)
    d['median_mapping_qualities'] = np.median(median_mapping_qualities)
    d['median_read_base_qualities'] = np.median(median_read_base_qualities)
    d['n_UMIs'] = len(UMIs)
    d['median_nUMI_group_size'] = np.median(UMI_group_sizes)
    
    # TO FIX AND CONTROL
    # mapping = { 'A' : 0, 'C' : 1, 'G' : 2, 'T': 3, 'N': 4} 
    # X = np.array([ [ mapping[x] for x in umi ] for umi in UMIs.keys() ])
    # D = pairwise_distances(X, metric='hamming') * 12
    # (D==1).sum(axis=1)

    return UMIs, pd.Series(d), pd.Series(UMI_group_sizes)


##


# pos = 208248389
# chrom = 'chr2'
# wt_allele = 'G'
# mut_allele = 'A'


def get_allele_counts(UMIs, UMI_group_sizes=None, chrom=None, pos=None, wt_allele=None, mut_allele=None, consensus_thr=.75, min_reads=3):
    """
    Consensus pileup at a target genomic position.
    """

    filtered_umis = UMI_group_sizes.loc[lambda x: x>=min_reads].index

    wt = 0
    mut = 0
    mis = 0

    for umi in filtered_umis:

        L = UMIs[umi]  
        bases = [ find_query_base_at_position(x, chrom, pos) for x in L ]
        bases = pd.Series(bases)
        freqs = bases.value_counts(normalize=True)

        if freqs.size == 1:
            consensus_base = freqs.index[0]
        elif freqs.values[0] >= consensus_thr:
            consensus_base = freqs.index[0]
        else:
            consensus_base = np.nan
        
        mismatched_alleles = set(['A', 'C', 'T', 'T'])-set([wt_allele, mut_allele])
        if consensus_base == wt_allele:
            wt += 1
        elif consensus_base == mut_allele:
            mut += 1
        elif consensus_base in mismatched_alleles:
            mis += 1

    return wt, mut, mis
        

##


# gene = 'IDH1'
# chrom = 'chr2'
# pos = 208248389
# wt_allele = 'G'
# mut_allele = 'A'
# min_reads =  3
# consensus_thr = .75
# positions = pd.read_csv(path_bed, sep='\t', header=None)

def make_allelic_table(UMIs, positions, UMI_group_sizes=None):

    MUT = []
    WT  = []
    MIS = []
    GENES = []

    for i in range(positions.shape[0]):
        chrom, pos, _, wt_allele, mut_allele, gene = positions.iloc[i,:].to_list()
        wt, mut, mis = get_allele_counts(UMIs, UMI_group_sizes=UMI_group_sizes, chrom=chrom, pos=pos, wt_allele=wt_allele, mut_allele=mut_allele)
        WT.append(wt)
        MUT.append(mut)
        MIS.append(mis)
        GENES.append(gene)

    allele_table = pd.DataFrame({'MUT':MUT, 'WT':WT, 'MIS':MIS}, index=GENES)

    return allele_table


##


# Filter UMIs

# Test
# path_ = '/Users/IEO5505/Desktop/example_mito/scratch'
# os.chdir(path_)
# cell = 'AAA'
# path_bam = os.path.join(path_, 'pre_filtered_nanopore.bam')
# path_bed = os.path.join(path_, 'sAML1.bed')
# base_quality_thr = 25                                   # 20 --> 1%, 30, 0.1 % errors --> NOT LESS THEN 25-30
# mapping_quality_thr = 30                                # 20 --> 1%, 30, 0.1 % errors --> NOT LESS THEN 30
# base_depth_thr = 3
# base_consensus_error_thr = .2


##


def main():

    cell = sys.argv[1]
    path_bam = sys.argv[2]
    path_bed = sys.argv[3]

    positions = pd.read_csv(path_bed, sep='\t', header=None)
    UMIs, d, UMI_group_sizes = filter_UMIs(path_bam, positions)
    allelic_table = make_allelic_table(UMIs, positions=positions, UMI_group_sizes=UMI_group_sizes)

    d.to_frame('value').reset_index().rename(columns={'index':'metric'}).assign(cell=cell).to_csv(f'{cell}_consensus_stats.csv', index=False, header=None)
    allelic_table.reset_index().rename(columns={'index':'gene'}).assign(cell=cell).to_csv(f'{cell}_allelic_table.csv', index=False, header=None)


##


# Run
if __name__ == '__main__' :
    main()





# selected_umis = UMI_group_sizes[UMI_group_sizes>=3].index
# 
# bam_in = pysam.AlignmentFile(path_bam, "rb", require_index=False)
# bam_out = pysam.AlignmentFile('pre_filtered_nanopore.bam', "wb", template=bam_in)

# for read in bam_in:
#     if read.get_tag('UB') in selected_umis:
#         bam_out.write(read)
# 
# bam_in.close()
# bam_out.close()

# UMI_group_sizes.size


##



