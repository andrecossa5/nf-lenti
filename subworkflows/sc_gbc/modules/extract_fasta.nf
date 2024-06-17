// EXTRACT_FASTA module

nextflow.enable.dsl = 2

//

process EXTRACT_FASTA {
    
    input:
    val(pattern)
    
    output:
    path("${pattern}.fa"), emit: fasta

    script:
    """
    echo ">${pattern}" > name.txt
    sed -n '/^>${pattern}/,/^>/ { /^>/! p }' ${params.ref}/*.fa > fasta.fa
    cat name.txt fasta.fa > ${pattern}.fa
    samtools faidx ${pattern}.fa
    samtools dict ${pattern}.fa > ${pattern}.dict

    """

    stub:
    """
    touch "${pattern}.fa"
    """

} 
 