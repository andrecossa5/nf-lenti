// SPLIT_BAM module

nextflow.enable.dsl = 2

//

process SPLIT_BAM {

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path('*.bam'), emit: cell_bams

    script:
    """ 
    samtools split -M -1 -@ ${task.cpus} -d CB -f '%!.bam' ${bam}
    rm -f *lentibam.* *mitobam.*
    rm -f -- -.bam
    """

    stub:
    """
    python ${baseDir}/bin/sc_gbc/random_mock.py \
    rm -f *filtered_bam_cb*
    """

}
