// FILTER_MITOBAM module

nextflow.enable.dsl = 2

//

process FILTER_BAM_CB {

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(bam), path(filtered)
    
    output:
    tuple val(sample_name), path("filtered_bam_cb.bam"), emit: bam

    script:
    """
    python ${baseDir}/bin/common/filter_bam_cb.py ${bam} filtered_bam_cb.bam ${filtered}/barcodes.tsv.gz
    """

    stub:
    """
    touch filtered_bam_cb.bam
    """

} 
