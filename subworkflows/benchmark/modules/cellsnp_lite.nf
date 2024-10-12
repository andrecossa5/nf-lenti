// freebayes module

nextflow.enable.dsl = 2

//

process CELLSNP {

    tag "${sample_name}"
    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), path(bam), path(barcodes)

    output:
    tuple val(sample_name), path('out_cellsnp'), emit: ch_output

    script:
    """ 
    samtools index ${bam}
    cellsnp-lite -s ${bam} -O out_cellsnp -p ${task.cpus} --minMAF 0.1 --minCOUNT 100 --chrom M -b ${barcodes}
    """

    stub:
    """
    mkdir out_cellsnp
    """

}