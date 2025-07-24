// FILTER_MITOBAM module

nextflow.enable.dsl = 2

//

process FILTER_BAM_CB {

    label 'scLT'
    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(bam), path(CBs)
    
    output:
    tuple val(sample_name), path("*.bam"), emit: bam

    script:
    """
    python ${baseDir}/bin/common/filter_bam_cb.py ${bam} ${CBs}
    """

    stub:
    """
    touch filtered.bam
    echo ${bam} > filtered.bam
    """

} 
