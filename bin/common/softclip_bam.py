



def filter_softclipping(softclipping_thr, input_bam, output_bam, drop_bam):
    # Import necessary modules
    import os
    import pysam

    # Open the input BAM file for reading
    bam_i = pysam.AlignmentFile(input_bam,  "rb")
    # Create an output BAM file for writing filtered reads, using the same header as the input BAM
    bam_o = pysam.AlignmentFile(output_bam, "wb", template=bam_i)

    # Dictionary to store reads that will be dropped due to excessive soft-clipping
    clipped_rid_dict = dict()

    # Iterate through all reads in the input BAM file
    for read in bam_i.fetch():
        # Initialize counters for matched (M) and soft-clipped (S) bases
        no_m, no_s = 0, 0

        # Loop through the CIGAR tuples of the read
        for cigar in read.cigartuples:
            # CIGAR operation 0 represents a match (M)
            if cigar[0] == 0:
                no_m += cigar[1]
            # CIGAR operations 4 and 5 represent soft-clipping (S) and hard-clipping (H)
            if cigar[0] == 4 or cigar[0] == 5:
                no_s += cigar[1]

        # Calculate the fraction of matched bases (M / (M + S))
        if (no_m / (no_m + no_s)) >= softclipping_thr:
            # If the fraction is greater than or equal to the threshold, write the read to the output BAM
            bam_o.write(read)
        else:
            # Otherwise, add the read to the dictionary of clipped reads
            clipped_rid_dict.update({read.query_name: 1})

    # Close the input and output BAM files
    bam_i.close()
    bam_o.close()

    # Re-open the input BAM file for reading
    bam_i = pysam.AlignmentFile(input_bam,  "rb")
    # Create a new BAM file for writing reads that were dropped
    bam_d = pysam.AlignmentFile(drop_bam,   "wb", template=bam_i)

    # Iterate through all reads again to find and write the dropped reads
    for read in bam_i.fetch():
        # If the read's query name is in the clipped reads dictionary, write it to the drop BAM file
        if read.query_name in clipped_rid_dict:
            bam_d.write(read)

    # Close the input and drop BAM files
    bam_i.close()
    bam_d.close()
