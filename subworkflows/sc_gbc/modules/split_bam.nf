// SPLIT_BAM module

nextflow.enable.dsl = 2

//

process SPLIT_BAM {

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(index)

    output:
    //tuple val(sample_name), path('cell_bams/*'), emit: cell_bams
    tuple val(sample_name), path('output/*'), emit: cell_bams

    script:
    """
    mkdir output 
    cd output
    samtools split -M -1 -u unrecognized.bam -d CB -f '%!.bam' ../${bam}
    """

    stub:
    """
    mkdir cell_bams
    cd cell_bams
    mkdir AAAA BBBB
    """

}
