// EXTRACT_FASTA module

nextflow.enable.dsl = 2

//

process EXTRACT_FASTA {
    
    input:
    val(pattern)
    
    output:
    tuple path("${pattern}.fa"), path("${pattern}.dict"),  path("${pattern}.fa.amb"),  path("${pattern}.fa.ann"),  path("${pattern}.fa.bwt"),  path("${pattern}.fa.fai"),  path("${pattern}.fa.pac"),  path("${pattern}.fa.sa"), emit: fasta

    script:
    """
    echo ">${pattern}" > name.txt
    sed -n '/^>${pattern}/,/^>/ { /^>/! p }' ${params.ref}/*.fa > fasta.fa
    cat name.txt fasta.fa > ${pattern}.fa
    samtools faidx ${pattern}.fa
    samtools dict ${pattern}.fa > ${pattern}.dict
    bwa index  ${pattern}.fa
    #creare cartella bwa qualcosa mv i tre file e fai uscire la cartella essere 
    """

    stub:
    """
    touch "${pattern}.fa"
    """

} 
 