// GATHER_ALLELIC_TABLES module

nextflow.enable.dsl = 2

process GATHER_TABLES {

    tag "${sample}"

    input:
    tuple val(sample), path(files)

    output:
    tuple val(sample), path("T_allelic_tables_cell.tsv.gz"), emit: elements_T

    script:
    """
    """

    stub:
    """
    for ext in "T" "G" "A" "C" "coverage"; do
        touch "\${ext}_allelic_tables_cell.tsv.gz"
    done
    """
}