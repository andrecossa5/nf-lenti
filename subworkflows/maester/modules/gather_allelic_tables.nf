// GATHER_ALLELIC_TABLES module

nextflow.enable.dsl = 2

process GATHER_ALLELIC_TABLES {

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
    files=(${files[@]})

    for ext in "T" "G" "A" "C" "coverage"; do
        touch "\${ext}_cells.txt"
    done

    #for f in "\${files[@]}"; do
    #    if [[ "\$f" == *.T.txt ]]; then
    #        cat "\$f" >> T_cells.txt
    #    elif [[ "\$f" == *.G.txt ]]; then
    #        cat "\$f" >> G_cells.txt
    #    elif [[ "\$f" == *.A.txt ]]; then
    #        cat "\$f" >> A_cells.txt
    #    elif [[ "\$f" == *.C.txt ]]; then
    #        cat "\$f" >> C_cells.txt
    #    elif [[ "\$f" == *.coverage.txt ]]; then
    #        cat "\$f" >> coverage_cells.txt
    #    fi
    #done

    #for ext in "T" "G" "A" "C" "coverage"; do
    #    sed 's/,/\t/g' \${ext}_cells.txt > \${ext}_allelic_tables_cell.tsv
    #    gzip --fast \${ext}_allelic_tables_cell.tsv
    #done
    """

    stub:
    """
    for ext in "T" "G" "A" "C" "coverage"; do
        touch "\${ext}_allelic_tables_cell.tsv.gz"
    done
    """
}
