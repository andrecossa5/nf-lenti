// FILTER LENTIBAM module

nextflow.enable.dsl = 2

//

process FILTER_LENTIBAM {

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(bam), path(barcodes)

    output:
    tuple val(sample_name), path("filtered_lentibam.bam"), emit: filtered_lentibam

    script:
    """
    python ${baseDir}/bin/sc_gbc/filter_lentibam.py ${bam} filtered_lentibam.bam ${barcodes}
    """

    stub:
    """
    touch filtered_lentibam.bam
    """

} 
 