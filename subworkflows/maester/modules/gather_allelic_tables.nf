// GATHER_ALLELIC_TABLES module

nextflow.enable.dsl = 2

process GATHER_TABLES {

    tag "${sample}"

    input:
    tuple val(sample), val(cells), path(A), path(C), path(T), path(G), path(coverage)

    output:
    tuple val(sample), path(tables), emit: tables

    script:
    """ 
    cat *.A.txt | gzip --fast > A.txt.gz
    cat *.C.txt | gzip --fast > C.txt.gz
    cat *.T.txt | gzip --fast > T.txt.gz
    cat *.G.txt | gzip --fast > G.txt.gz
    cat *.coverage.txt | gzip --fast > coverage.txt.gz
    mkdir tables 
    mv *.gz tables
    """

    stub:
    """
    mkdir tables
    """

}