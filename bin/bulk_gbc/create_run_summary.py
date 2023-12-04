#!/usr/bin/python

# Create run summary

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='create_run_summary',
    description=
    """
    Create run summary.
    """
)

# indir
my_parser.add_argument(
    '--indir', 
    type=str,
    default=None,
    help='Path to input directory: params.bulk_indir.'
)

# outdir
my_parser.add_argument(
    '--outdir', 
    type=str,
    default=None,
    help='Path to input directory: params.bulk_outdir.'
)

# anchor_sequence
my_parser.add_argument(
    '--anchor_sequence', 
    type=str,
    default=None,
    help='Anchor sequence: params.anchor_sequence.'
)

# read_counts
my_parser.add_argument(
    '--raw_counts',
    type=str,
    default=None,
    help='Path to raw_counts, all observed GBC sequences. Default: . .'
)

# Sample
my_parser.add_argument(
    '--sample',
    type=str,
    default=None,
    help='Sample name.'
)

# correction_df
my_parser.add_argument(
    '--correction_df',
    type=str,
    default=None,
    help='Correction_df as generated in correct_and_count.py script.'
)

# stats_table
my_parser.add_argument(
    '--corrected_counts',
    type=str,
    default=None,
    help='Corrected counts for filtered, high quality GBC sequences.'
)

# 
my_parser.add_argument(
    '--output',
    type=str,
    default=os.getcwd(),
    help='Path to write output. Default: cwd.'
)

# Min_n_reads
my_parser.add_argument(
    '--min_n_reads',
    type=int,
    default=1000,
    help='Min number of reads for a filtered GBC. Default: 1000.'
)

# Min_n_reads
my_parser.add_argument(
    '--hamming_treshold',
    type=int,
    default=3,
    help='Hamming treshold to find GBCs communities via UMItools clusterer. Default: 3.'
)


##


# Parse arguments
args = my_parser.parse_args()
indir = args.indir
outdir = args.outdir
sample = args.sample
anchor_sequence = args.anchor_sequence
min_n_reads = args.min_n_reads
hamming_treshold = args.hamming_treshold
raw_counts = args.raw_counts
corrected_counts = args.corrected_counts
correction_df = args.correction_df
path_o = args.output


##


# Import code
import getpass
import datetime
import pandas as pd


##


########################################################################

def main():

    # Calculate stats

    df_ = pd.read_csv(raw_counts, index_col=0)
    total_GBC_reads = df_['read_counts'].sum()
    n_unique_GBCs = df_.shape[0]
    del df_

    df_ = pd.read_csv(corrected_counts, index_col=0)
    n_filtered_GBCs = df_.shape[0]
    del df_

    df_ = pd.read_csv(correction_df)
    n_corrected = df_['correct'].unique().size
    n_degenerated = df_['degenerated'].unique().size
    median_n_reads_added = df_['n_reads_degenerated'].median()
    del df_


    ##


    # Write summary
    f = open(os.path.join(path_o, 'run_summary.txt'), 'w')

    # Write stats
    f.write(f'Summary Step 1 (bulk GBC), sample {sample} \n')
    f.write('-------------------------------------')
    f.write('\n')
    f.write('Overview: \n')
    f.write(f'- Date of analysis:               {datetime.datetime.now().strftime("%d-%m-%Y")} \n')
    f.write(f'- User:                           {getpass.getuser()} \n')
    f.write(f'- Working directory:              {"/".join(os.getcwd().split("/")[:-2])} \n')
    f.write('\n')
    f.write(f'Parameters \n')
    f.write(f'--indir:                          {indir} \n')
    f.write(f'--outdir:                         {outdir} \n')
    f.write(f'--anchor_sequence:                {anchor_sequence} \n')
    f.write(f'--min_n_reads:                    {min_n_reads} \n')
    f.write(f'--hamming_treshold:               {hamming_treshold} \n')
    f.write('\n')
    f.write('Numbers: \n')
    f.write(f'- Reads with GBC:                 {int(total_GBC_reads)} \n')
    f.write(f'- Unique GBCs:                    {int(n_unique_GBCs)} \n')
    f.write(f'- Filtered GBCs:                  {int(n_filtered_GBCs)} \n')
    f.write(f'- Correction:                     \n')
    f.write(f'  * n corrected GBCs:             {n_corrected} \n')
    f.write(f'  * n degenerated GBCs:           {n_degenerated} \n')
    f.write(f'  * Median n added reads:         {median_n_reads_added:.2f} \n')
    f.write('\n')

    f.close()
 

##


########################################################################

# Run
if __name__ == '__main__':
    main()