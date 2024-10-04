// freebayes module

nextflow.enable.dsl = 2

//

process FREEBAYES {

    tag "${sample_name}: ${cell}"
    
    input:
    tuple val(sample_name), val(cell), path(bam)
    path(reference)

    output:
    tuple val(sample_name), path('filtered.tsv'), emit: calls

    script:
    """ 
    # ...
    """

    stub:
    """
    touch filtered.tsv
    """

}

//

process COLLAPSE_FREEBAYES {

    tag "${sample_name}"

    input:
    tuple val(sample_name), val(muts)

    output:
    tuple val(sample_name),
        path('AD.mtx'), 
        path('DP.mtx'),
        path('cells.txt'), emit: matrices

    script:
    """ 
    # ...
    """

    stub:
    """
    touch AD.mtx
    touch DP.mtx
    touch cells.txt
    """

}
