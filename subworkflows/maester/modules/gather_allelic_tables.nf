// GATHER_ALLELIC_TABLES module

nextflow.enable.dsl = 2

process GATHER_TABLES {

    tag "${sample}"

    input:
    tuple val(sample), path(files)

    output:
    tuple val(sample), path("T_allelic_tables_cell.tsv.gz"), emit: elements_T
    tuple val(sample), path("G_allelic_tables_cell.tsv.gz"), emit: elements_G
    tuple val(sample), path("A_allelic_tables_cell.tsv.gz"), emit: elements_A
    tuple val(sample), path("C_allelic_tables_cell.tsv.gz"), emit: elements_C
    tuple val(sample), path("coverage_allelic_tables_cell.tsv.gz"), emit: elements_coverage

    script:
    """
    #!/bin/bash
    for ext in T G A C coverage; do
    touch "\${ext}_cells.txt"
    done
    """

    stub:
    """
    for ext in "T" "G" "A" "C" "coverage"; do
        touch "\${ext}_allelic_tables_cell.tsv.gz"
    done
    """
}