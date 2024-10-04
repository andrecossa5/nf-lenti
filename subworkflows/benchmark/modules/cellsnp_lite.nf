// freebayes module

nextflow.enable.dsl = 2

//

process CELLSNP {

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(barcodes)

    output:
    tuple val(sample_name),
    path('AD.mtx'), 
    path('DP.mtx'),
    path('cells.txt'), emit: matrices

    script:
    """ 
    samtools index ${bam}
    cellsnp-lite -s ${bam} -O out_cellsnp -p 10 --minMAF 0.1 --minCOUNT 100 --chrom M -b ${barcodes}
    """

    stub:
    """
    touch AD.mtx
    touch DP.mtx
    touch cells.txt
    """

}