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
    '--base_quality_th', 
    type=int,
    default=None,
    help='Quality threshold of each base expressed with Phred score, bases with lower quality are masked with N '
)

my_parser.add_argument(
    '--mean_quality_th', 
    type=int,
    default=None,
    help='Mean quality of the whole reads expressed with Phred score, reads with lower mean quality are filtered out '
)

my_parser.add_argument(
    '--E', 
    type=float,
    default=None,
    help='The maximum raw-read error rate across the entire consensus read. Reads with higher erreor are filtered out '
)

my_parser.add_argument(
    '--e', 
    type=float,
    default=None,
    help='The maximum error rate for a single consensus base. bases with higer error are masked with N'
)

my_parser.add_argument(
    '--mask_th', 
    type=float,
    default=None,
    help='The maximum density of N in a sequence after masking.'
)

##

# Parse arguments
args = my_parser.parse_args()
bam_file = args.input
output_bam = args.output
quality_th = args.base_quality_th
mean_qual_th = args.mean_quality_th
max_Er = args.E
base_er = args.e
mask_th = args.mask_th

bam_file = '/Users/ieo6943/Downloads/cell_folder/consensus.bam'
output_bam = '/Users/ieo6943/Downloads/cell_folder/filter_local.bam'


def filter_consensus(bam_file, output_bam, quality_th, mean_qual_th,max_Er,base_er,mask_th):

    bam_in = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)
    df = pd.DataFrame(columns=['read', 'ce', 'quality', 'supported'])
    j=0
    for read in bam_in.fetch(until_eof=True):
        
        n_consensus_read = read.get_tag('cM:')
        mean_consensus_error = read.get_tag('cE:')
        base_consensus_error = np.array(read.get_tag('ce:'))/n_consensus_read
        high_consensus_err_index = np.where(base_consensus_error>base_er)[0]
        

        sequence = list(read.query_sequence)
        qualities = np.array(read.query_qualities)
        mean_qualities = qualities.sum()/qualities.shape[0]

        
        
        #MASKs for the bad bases #mettere np.where per la maschera
        low_quality_indices = np.where(qualities < quality_th)[0]
        # Replace low quality bases with 'N'
        for i in low_quality_indices:
            sequence[i] = 'N'
        for i in high_consensus_err_index:
            sequence[i] = 'N'
        bad_reads_density=sequence.count('N')/len(sequence)
        #filter if the mask is too high
        supported = 'not supported'

        if bad_reads_density<mask_th and mean_consensus_error<max_Er and mean_qualities>mean_qual_th:
            supported = 'supported'


        #save metadata for filtering analysis    
        df.loc[j] = [sequence, base_consensus_error, qualities , supported]


        if bad_reads_density>mask_th:
            continue

        #filter mean error
        if mean_consensus_error>max_Er:
            continue

        # filter mean quality
        if mean_qualities<mean_qual_th: 
            continue

        read.query_sequence = ''.join(sequence)
        bam_out.write(read)
        j+=1
        #if j==10: break

    df.to_csv('seq_qual_err.csv',sep='\t')
    bam_in.close()
    bam_out.close()

def main():
        filter_consensus(bam_file, output_bam, quality_th, mean_qual_th,max_Er,base_er,mask_th)




if __name__ == "__main__":
    main()


#PEr ogni UMI che passa nel filter bam srivere la size del UMI fare un dizionario? 