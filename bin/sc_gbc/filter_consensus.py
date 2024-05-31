import pysam
import numpy as np
import pandas as pd
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='cell_assignment',
    description=
    """
    Script for clone calling and cell assignment.
    """
)

my_parser.add_argument(
    '--input', 
    type=str,
    default=None,
    help='Path input bam. should come from fgbio CallMolecularConsensus'
)

my_parser.add_argument(
    '--output', 
    type=str,
    default=None,
    help='Path output bam'
)

my_parser.add_argument(
    '--min_quality', 
    type=30,
    default=None,
    help='Mean quality of the whole reads expressed with Phred score, reads with lower mean quality are filtered out '
)

my_parser.add_argument(
    '--read_max_consensus_error', 
    type=float,
    default=0.1,
    help='The maximum raw-read error rate across the entire consensus read. Reads with higher erreor are filtered out '
)

my_parser.add_argument(
    '--base_consensus_error', 
    type=float,
    default=0.1,
    help='The maximum error rate for a single consensus base. bases with higer error are masked with N'
)

my_parser.add_argument(
    '--read_max_N', 
    type=float,
    default=0.2,
    help='The maximum density of N in a sequence after masking.'
)

my_parser.add_argument(
    '--consensus_filter_mode', 
    type=str,
    default=None,
    help='The maximum density of N in a sequence after masking.'
)

my_parser.add_argument(
    '--GBC_max_N', 
    type=float,
    default=0,
    help='The maximum density of N in a sequence after masking.'
)
##

# Parse arguments
args = my_parser.parse_args()
bam_file = args.input
output_bam = args.output
min_quality = args.min_quality
read_max_consensus_error = args.read_max_consensus_error
base_consensus_error = args.base_consensus_error
read_max_N = args.read_max_N
consensus_filter_mode = args.consensus_filter_mode
GBC_max_N =args.GBC_max_N



def filter_consensus(bam_file, output_bam, min_quality,read_max_consensus_error,base_consensus_error,read_max_N,consensus_filter_mode,GBC_max_N):



    bam_in = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)
    df = pd.DataFrame(columns=['read', 'ce', 'quality','UMI','n_read_consensus', 'SUPPORT'])
    GBC_bases = True

    j=0
    for read in bam_in:

        n_consensus_read = read.get_tag('cM:')
        mean_consensus_error = read.get_tag('cE:')
        base_consensus_error = np.array(read.get_tag('ce:'))/n_consensus_read
        high_consensus_err_index = np.where(base_consensus_error>base_consensus_error)[0]
        umi = read.get_tag('MI:')

        sequence = list(read.query_sequence)
        qualities = np.array(read.query_qualities)
        mean_qualities = qualities.sum()/qualities.shape[0]

        #MASKs for the bad bases #mettere np.where per la maschera
        low_quality_indices = np.where(qualities<min_quality)[0]
        # Replace low quality bases with 'N'
        for i in low_quality_indices:
            sequence[i] = 'N'
        for i in high_consensus_err_index:
            sequence[i] = 'N'
        bad_reads_density=sequence.count('N')/len(sequence)
        sequence = np.array(sequence)
        np.where(sequence=='N',1,0).sum()

        GBC_n_N = np.where(sequence[33:33+18]=='N',1,0).sum()
        #GBC_n_N = 0

        supported = 'not supported'

        if mean_consensus_error<read_max_consensus_error and mean_qualities>=min_quality and bad_reads_density<read_max_N :
            if consensus_filter_mode=='GBC':
                if GBC_n_N==GBC_max_N:
                    supported = 'supported'

            supported = 'supported'

        if supported=='supported':
            read.query_sequence = ''.join(sequence)
            bam_out.write(read)

        final_seq = ''.join(sequence)
        #save metadata for filtering analysis    
        df.loc[j] = [final_seq, base_consensus_error, qualities,umi,n_consensus_read, supported]
        j+=1
    #if j==10: break
    bam_in.close()
    bam_out.close()

def main():
        filter_consensus(bam_file, output_bam, min_quality,read_max_consensus_error,base_consensus_error,read_max_N,consensus_filter_mode,GBC_max_N)




if __name__ == "__main__":
    main()


#PEr ogni UMI che passa nel filter bam srivere la size del UMI fare un dizionario? 