// COLLAPSE_OUTPUT  module
 
nextflow.enable.dsl = 2
 
//


process COLLAPSE_TSV {
 
    tag "${cells}"
 
    input:
        path(files)
 
    output:
        path("cells.tsv.gz"), emit: elements
 
    script:
    """
    outfile="cells.tsv"
    files=(${files})
    cat "\${files[0]}" > \$outfile
    for f in "\${files[@]:1}"; do
        tail -n +2 "\$f" >> \$outfile
    done
    gzip --fast \$outfile
    """

    stub:
    """
    touch cells.tsv.gz
    """
 
}