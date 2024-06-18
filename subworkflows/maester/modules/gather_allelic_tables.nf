// GATHER_ALLELIC_TABLES  module
 
nextflow.enable.dsl = 2
 
//


process GATHER_ALLELIC_TABLES {
 
    tag "${sample}"
 
    input:
        tuple val(sample), path(files)
 
    output:
        tuple val(sample), path("allelic_tables_cell.tsv.gz"), emit: elements
 
    script:
    """
    files=(${files})
    touch "allelic_tables_cell.txt"
    `
    for f in "${files[@]}"; do
        cat "$f" >> cells.txt
    done
    sed 's/,/\t/g' cells.txt > allelic_tables_cell.tsv
    gzip --fast allelic_tables_cell.tsv
    """

    stub:
    """
    touch allelic_tables_cell.tsv.gz
    """
}