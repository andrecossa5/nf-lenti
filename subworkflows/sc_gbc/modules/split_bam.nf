// SPLIT_BAM module

nextflow.enable.dsl = 2

//

process SPLIT_BAM {

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(index), path(filtered)

    output:
    //tuple val(sample_name), path('cell_bams/*'), emit: cell_bams
    tuple val(sample_name), path('output/*'), emit: cell_bams

    beforeScript 'source ~/.bashrc && mamba activate NextFlow_G'

    script:
    """
    #python ${baseDir}/bin/sc_gbc/split_bam.py ${bam} ${filtered}/barcodes.tsv.gz
    mkdir output 
    cd output
    python ${baseDir}/bin/sc_gbc/filter_lentibam ${bam} filtered_${bam} ${filtered}/barcodes.tsv.gz
    samtools split -u unrecognized.bam -d CB -f '%!.bam' ../filtered_${bam}
    """

    stub:
    """
    mkdir cell_bams
    cd cell_bams
    mkdir AAAA BBBB
    """

}
