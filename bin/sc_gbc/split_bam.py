#!/usr/bin/python

import os
import sys
import pandas as pd
import pysam


##


# Args
path_bam = sys.argv[1]
path_cells = sys.argv[2]
# path_cells = '/Users/IEO5505/Desktop/example_mito/scratch_data/barcodes.tsv.gz'
# path_bam = '/Users/IEO5505/Desktop/example_mito/scratch_data/Aligned.sortedByCoord.out.bam'
# os.chdir('/Users/IEO5505/Desktop/example_mito/scratch_data/')


##


def main():
    """
    Sript to split a .bam file into a cell_bams folder with <CB> subfolders each /<CB>.bam,
    with <CB> being the STAR-Solo error-corrected cellular barcode for an alignment record.
    """

    # Prep cell writers and folders
    output_folder = 'cell_bams'
    os.makedirs(output_folder, exist_ok=True)
    cbc_writers = {}
    
    # Read good quality, STAR-Solo corrected CBs
    cells = pd.read_csv(path_cells, header=None, index_col=0)

    # Parse .bam and write it
    with pysam.AlignmentFile(path_bam, "rb") as bam_input:
        for alignment in bam_input:
            cbc = alignment.get_tag('CB')
            # If CB in STAR-solo CBs
            if cbc in cells.index:
                bam_cbc_path = os.path.join(output_folder, f'{cbc}/')
                os.makedirs(bam_cbc_path, exist_ok=True)
                if cbc not in cbc_writers:
                    cbc_writers[cbc] = pysam.AlignmentFile(
                        os.path.join(bam_cbc_path, f"{cbc}.bam"), 
                        'wb', header=bam_input.header
                    )
                cbc_writers[cbc].write(alignment)

    # Close all streams
    for writers in cbc_writers.values():
        writers.close()


##


# Run 
if __name__ == '__main__':
    main()


