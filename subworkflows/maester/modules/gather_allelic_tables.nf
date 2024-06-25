// GATHER_ALLELIC_TABLES module

nextflow.enable.dsl = 2

process GATHER_TABLES {

    tag "${sample}"

    input:
    tuple val(sample), path(A), path(C), path(T), path(G), path(coverage)

    // .map { it -> tuple(it[0], it[2], it[3], it[4], it[5], it[6]) }
    // .groupTuple(by: 0)

    output:
    tuple val(sample), 
        path("T.txt.gz"), 
        path("G.txt.gz"),
        path("A.txt.gz"), 
        path("C.txt.gz"),
        path("coverage.txt.gz"), emit: tables

    script:
    """ 
    cat *.A.txt | gzip --fast > A.txt.gz
    cat *.C.txt | gzip --fast > C.txt.gz
    cat *.T.txt | gzip --fast > T.txt.gz
    cat *.G.txt | gzip --fast > G.txt.gz
    cat *.coverage.txt | gzip --fast > coverage.txt.gz
    """

    stub:
    """
    touch T.txt.gz 
    touch G.txt.gz
    touch A.txt.gz 
    touch C.txt.gz
    touch coverage.txt.gz
    """

}