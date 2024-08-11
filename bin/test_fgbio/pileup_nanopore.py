#!/usr/bin/python

"""
Nanopore consensus.
"""

import os
import sys
import numpy as np
import pandas as pd
import pysam


##


# Utils

# target_position = 47120267
# target_reference_name = 'chr3'
# read = L[0]

def find_query_base_at_position(read, target_reference_name, target_position):

    if read.reference_name != target_reference_name:
        return (np.nan, np.nan, np.nan)              # The read is not aligned to the desired chromosome        
    
    seq = read.seq
    base_qualities = read.query_qualities 
    strand = 'fw' if read.mate_is_forward else 'rev'
    
    found_position = False
    for qpos, refpos in read.get_aligned_pairs(True):
        if refpos == (target_position-1):           # PD
            a = seq[qpos]
            q = base_qualities[qpos]
            found_position = True
            break
    
    record = ( a, q, strand ) if found_position else (np.nan, np.nan, np.nan)  

    return record


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

        seq_test = mapping_quality >= mapping_quality_thr and BQ >= median_base_quality_thr and read.reference_name in chromosomes
    
        if seq_test:                                    # Entire sequence test 
            
            n_reads_filtered += 1
            UB = read.get_tag('UB')
            if UB in UMIs:
                UMIs[UB].append((read.reference_name, read))
            else:
                UMIs[UB] = [(read.reference_name, read)]
        
    # Record .bam stats
    UMI_group_sizes = [ len(UMIs[k]) for k in UMIs ]
    chrom = [ [ x[0] for x in UMIs[k] ] for k in UMIs ]
    chrom = [ ch[0] if all(x == ch[0] for x in ch) else np.nan for ch in chrom ]
    UMI_stats = pd.DataFrame({'n_reads':UMI_group_sizes, 'chrom':chrom}, index=UMIs.keys())
    d = {}
    d['n_reads_filtered'] = n_reads_filtered
    d['n_reads_unfiltered'] = n_reads_unfiltered
    d['read_lengths'] = np.median(read_lengths)
    d['median_mapping_qualities'] = np.median(median_mapping_qualities)
    d['median_read_base_qualities'] = np.median(median_read_base_qualities)
    d['n_UMIs'] = len(UMIs)
    d['median_nUMI_group_size'] = np.median(UMI_stats['n_reads'])
    d['n_chrom'] = UMI_stats['chrom'].dropna().unique().size
    p = (UMI_stats['n_reads'] / UMI_stats['n_reads'].sum()).sort_values(ascending=False)
    n_UMIs_cumsum_75 = p.size-(p.cumsum()>=.75).sum()
    d['n_UMIs_cumsum_75'] = n_UMIs_cumsum_75 if n_UMIs_cumsum_75>0 else 1

    return UMIs, pd.Series(d), UMI_stats


##


# pos = 28018505
# chrom = 'chr13'
# wt_allele = 'C'
# mut_allele = 'A'

def get_allele_counts(UMIs, UMI_stats=None, chrom=None, pos=None, wt_allele=None,
                     mut_allele=None, consensus_thr=.75, min_reads=3, min_base_quality=25):
    """
    Consensus pileup at a target genomic position.
    """

    filtered_umis = UMI_stats.dropna().loc[lambda x: x['n_reads']>=min_reads].index

    wt = 0
    mut = 0
    mis = 0

    for umi in filtered_umis:

        L = [ x[1] for x in UMIs[umi] ]  
        bases = [ find_query_base_at_position(x, chrom, pos) for x in L ]
        bases = pd.DataFrame(bases, columns=['allele', 'q', 'strand'])

        consensus_base = np.nan
        if bases['allele'].isna().all():                                    # Temporary...
            if mut_allele[0] == '-':
                bases_minus_one = pd.DataFrame([ find_query_base_at_position(x, chrom, pos-1) for x in L ], columns=['allele', 'q', 'strand'])
                bases_plus_one = pd.DataFrame([ find_query_base_at_position(x, chrom, pos+1) for x in L ], columns=['allele', 'q', 'strand'])
                if bases_minus_one.dropna().shape[0] > 0 and bases_plus_one.dropna().shape[0] > 0:
                    consensus_base = mut_allele
            else:
                consensus_base = np.nan
        else:
            bases = bases.dropna().query('q>=@min_base_quality')
            if bases.shape[0]>0:
                freqs = bases['allele'].value_counts(normalize=False)
                freqs = freqs.to_frame('n').assign(f=lambda x: x['n']/x['n'].sum())
                test = freqs['n'][0]>=min_reads and freqs['f'][0]>=consensus_thr
                consensus_base = freqs.index[0] if test else np.nan
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

def make_allelic_table(UMIs, positions, UMI_stats=None, consensus_thr=.75, min_reads=3, min_base_quality=25):

    MUT = []
    WT  = []
    MIS = []
    GENES = []

    for i in range(positions.shape[0]):
        chrom, pos, _, wt_allele, mut_allele, gene = positions.iloc[i,:].to_list()
        wt, mut, mis = get_allele_counts(
            UMIs, UMI_stats=UMI_stats, chrom=chrom, pos=pos, wt_allele=wt_allele, mut_allele=mut_allele, 
            consensus_thr=consensus_thr, min_reads=min_reads, min_base_quality=min_base_quality
        )
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
# path_bam = os.path.join(path_, 'nanopore.bam')
# path_bed = os.path.join(path_, 'sAML1.bed')
# base_quality_thr = 25                                   # 20 --> 1%, 30, 0.1 % errors --> NOT LESS THEN 25-30
# mapping_quality_thr = 60                                # 20 --> 1%, 30, 0.1 % errors --> NOT LESS THEN 30
# min_reads = 1
# base_consensus_error_thr = .25


##


def main():

    cell = sys.argv[1]
    path_bam = sys.argv[2]
    path_bed = sys.argv[3]
    min_reads = int(sys.argv[4])
    base_consensus_error_thr = float(sys.argv[5])
    base_quality_thr = int(sys.argv[6])

    positions = pd.read_csv(path_bed, sep='\t', header=None)
    UMIs, d, UMI_stats = filter_UMIs(path_bam, positions)
    allelic_table = make_allelic_table(
        UMIs, positions=positions, UMI_stats=UMI_stats, min_reads=min_reads, 
        consensus_thr=1-base_consensus_error_thr, min_base_quality=base_quality_thr
    )
    d.to_frame('value').reset_index().rename(columns={'index':'metric'}).assign(cell=cell).to_csv(f'{cell}_consensus_stats.csv', index=False, header=None)
    allelic_table.reset_index().rename(columns={'index':'gene'}).assign(cell=cell).to_csv(f'{cell}_allelic_table.csv', index=False, header=None)


##


# Run
if __name__ == '__main__' :
    main()


##


# Controls

# BLAT
# chrom = 'chr21'
# pos = 34880691
# 
# UMI_stats.query('chrom==@chrom and n_reads>=3')
# 
# umi = 'TATACCCGCAAT'
# 
# seqs = [ x[1].seq for x in UMIs[umi] ]
# [ find_query_base_at_position(x[1], chrom, pos+1) for x in UMIs[umi] ]
# 
# n = 5
# i = 0
# 
# seqs[i]
# read = UMIs[umi][i][1]
# len(read.seq)
# 
# read.mate_is_forward
# read.query_alignment_start
# read.query_alignment_end
# read.reference_start
# 
# [ x for x in read.get_aligned_pairs(True) if x[1] == pos-2 ][0]
# [ x for x in read.get_aligned_pairs(True) if x[1] == pos ][0]
# 
# 
# qpos, refpos = [ x for x in read.get_aligned_pairs(True) if x[1] == pos-1 ][0]
# flanks = [ x for x in read.get_aligned_pairs(True) if x[0] > (qpos-n) and x[0] < (qpos+n) ]
# 
# fasta = pysam.FastaFile('GRCh38.d1.vd1.fa')
# 
# query_seq = [ read.seq[x[0]] for x in flanks ]
# ref_seq = [ fasta.fetch(chrom, start=x[1], end=x[1]+1) for x in flanks ] 
# 
# len(query_seq)
# len(ref_seq)
# 
# query_seq
# ref_seq
# read.seq[qpos]
# fasta.fetch(chrom, start=refpos, end=refpos+1)
# list(read.query_qualities[flanks[0][0]:flanks[-1][0]+1])


##