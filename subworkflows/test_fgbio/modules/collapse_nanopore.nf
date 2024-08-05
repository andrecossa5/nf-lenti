// COLLAPSE_NANOPORE module

nextflow.enable.dsl = 2

process COLLAPSE_NANOPORE {

    tag "${sample_name}"
    publishDir "${params.test_outdir}", mode: 'copy'

    input:
    tuple val(sample_name), val(cells), path(consensus_stats), path(allelic_tables)

    output:
    tuple val(sample_name), path("counts_table.csv"), path("cons_stats.csv"), emit: tables

    script:
    """ 
    # Gather
    echo metric,value,cell > header
    cat header *consensus_stats.csv >> cons_stats.csv
    echo gene,MUT,WT,MIS,cell > header
    cat header *allelic_table.csv > counts_table.csv
    """

    stub:
    """
    mkdir counts_table.csv
    touch cons_stats.csv
    """

}