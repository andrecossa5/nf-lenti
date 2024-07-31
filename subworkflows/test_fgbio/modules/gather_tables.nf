// GATHER_ALLELIC_TABLES module

nextflow.enable.dsl = 2

process GATHER_TABLES {

    tag "${sample_name}"

    input:
    tuple val(sample_name),
        val(cells), 
        path(A), 
        path(C), 
        path(T), 
        path(G), 
        path(median_base_umi_group_size),
        path(n_umis_unfiltered),
        path(n_umis_filtered),
        path(depth),
        path(coverage),
        path(median_filtered_base_consensus_error),
        path(median_filtered_read_quality),
        path(n_reads_filtered),
        path(n_reads_unfiltered)

    output:
    tuple val(sample_name), path(tables), emit: tables
    tuple val(sample_name), path("stats.csv"), emit: stats

    script:
    """ 
    # Gather
    cat *.A.txt | gzip --fast > A.txt.gz
    cat *.C.txt | gzip --fast > C.txt.gz
    cat *.T.txt | gzip --fast > T.txt.gz
    cat *.G.txt | gzip --fast > G.txt.gz
    cat *.depth.txt | gzip --fast > depth.txt.gz
    cat *.coverage.txt | gzip --fast > coverage.txt.gz
    mkdir tables 
    mv *.gz tables

    # Create stats.csv
    # ....                  DO things with other tables
    touch stats.csv
    """

    stub:
    """
    mkdir tables
    touch stats.csv
    """

}