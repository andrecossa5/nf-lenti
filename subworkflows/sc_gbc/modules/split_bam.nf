// SPLIT_BAM module

nextflow.enable.dsl = 2

//

process SPLIT_BAM {

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(index), path(filtered)

    output:
    tuple val(sample_name), path('cell_bams/*'), emit: cell_bams

    script:
    """
    python ${baseDir}/bin/sc_gbc/split_bam.py ${bam} ${filtered}/barcodes.tsv.gz
    """

    stub:
    """
    mkdir cell_bams
    cd cell_bams
    mkdir AAAA BBBB
    """

}
