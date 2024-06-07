// FILTER LENTIBAM module

nextflow.enable.dsl = 2

//

process FILTER_LENTIBAM {

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(bam), path(filtered)
    
    output:
    tuple val(sample_name), path("filtered_lentibam.bam"), path("filtered_lentibam.bam.bai"), emit: filtered_lentibam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    samtools view ${bam} -b -@ ${task.cpus} lentiCassette > lentibam.bam
    samtools index -@ ${task.cpus} lentibam.bam
    python ${baseDir}/bin/sc_gbc/filter_lentibam.py lentibam.bam filtered_lentibam.bam ${filtered}/barcodes.tsv.gz
    samtools index -@ ${task.cpus} filtered_lentibam.bam
    """

    stub:
    """
    touch filtered_lentibam.bam
    touch filtered_lentibam.bam.bai
    """

} 
