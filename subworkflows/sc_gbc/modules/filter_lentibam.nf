// FILTER LENTIBAM module

nextflow.enable.dsl = 2

//

process FILTER_LENTIBAM {

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(bam), path(filtered)

    output:
    tuple val(sample_name), path("filtered_lentibam.bam"), path("filtered_lentibam.bam.csi"), emit: filtered_lentibam

    script:
    """
    python ${baseDir}/bin/sc_gbc/filter_lentibam.py lentibam.bam filtered_lentibam.bam ${filtered}/barcodes.tsv.gz 
    """

    stub:
    """
    touch filtered_lentibam.bam
    touch filtered_lentibam.bam.csi 
    """

} 
 