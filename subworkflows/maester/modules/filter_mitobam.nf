// FILTER_MITOBAM module

nextflow.enable.dsl = 2

//

process FILTER_MITOBAM {

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(bam), path(filtered)
    
    output:
    tuple val(sample_name), path("filtered_mitobam.bam"), emit: filtered_mitobam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    python ${baseDir}/bin/maester/filter_mitobam.py ${filtered}/barcodes.tsv.gz
    """

    stub:
    """
    touch filtered_mitobam.bam
    """

} 
