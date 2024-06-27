// SPLIT_BAM module

nextflow.enable.dsl = 2

//

process SPLIT_BAM {

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(index)

    output:
    tuple val(sample_name), path('*.bam'), emit: cell_bams

    script:
    """
    samtools split -M -1 -@ ${task.ncpus} -d CB -f '%!.bam' ${bam}
    rm -f *lentibam.* *mitobam.*
    """

    stub:
    """
    touch AAAA.bam BBBB.bam
    rm -f *lentibam.* *mitobam.*
    """

}
