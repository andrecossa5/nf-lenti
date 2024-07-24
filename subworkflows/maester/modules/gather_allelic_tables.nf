// GATHER_ALLELIC_TABLES module

nextflow.enable.dsl = 2

process GATHER_TABLES {

    tag "${sample}"

    input:
    tuple val(sample), 
        val(cells), 
        path(A), 
        path(C), 
        path(T), 
        path(G), 
        path(median_base_umi_group_size)
        path(n_umis_unfiltered)
        path(n_umis_filtered)
        path(depth)
        path(coverage)

    output:
    tuple val(sample), path(tables), emit: tables

    script:
    """ 
    cat *.A.txt | gzip --fast > A.txt.gz
    cat *.C.txt | gzip --fast > C.txt.gz
    cat *.T.txt | gzip --fast > T.txt.gz
    cat *.G.txt | gzip --fast > G.txt.gz
    cat *.coverage.txt | gzip --fast > coverage.txt.gz
    cat *.median_base_umi_group_size.txt | gzip --fast > median_base_umi_group_size.txt.gz
    cat *.n_umis_unfiltered.txt | gzip --fast > n_umis_unfiltered.txt.gz
    cat *.n_umis_filtered.txt | gzip --fast > n_umis_filtered.txt.gz
    cat *.depth.txt | gzip --fast > depth.txt.gz
    mkdir tables 
    mv *.gz tables
    """

    stub:
    """
    mkdir tables
    """

}