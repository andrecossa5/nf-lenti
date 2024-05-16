// FILTER LENTIBAM module

nextflow.enable.dsl = 2

//

process FILTER_LENTIBAM {

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(bam)
    
    output:
    tuple val(sample_name), path("lentibam.bam"), path("lentibam.bam.bai"), emit: lentibam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    samtools view ${bam} -b -@ ${task.cpus} lentiCassette > lentibam.bam
    samtools index -@ ${task.cpus} lentibam.bam
    """

    stub:
    """
    touch lentibam.bam
    touch lentibam.bam.bai
    """

} 
