// freebayes modules

nextflow.enable.dsl = 2

//

process FREEBAYES {

    tag "${sample_name}: ${cell}"

    input:
    tuple val(sample_name), val(cell), path(bam)
    tuple path(reference), path(ref_idx)

    output:
    tuple val(sample_name), val(cell), path("${cell}_filtered.tsv"), emit: calls

    script:
    """ 
    picard MarkDuplicates I=${bam} O=cell_dedup.bam M=deduplication_metrics.txt REMOVE_DUPLICATES=true
    freebayes -C 0 -F 0 --fasta-reference ${reference} cell_dedup.bam > cell.vcf.gz
    bcftools filter  -i 'QUAL>20' -Oz -o filtered_variants.vcf.gz cell.vcf.gz
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%AD]\t[%DP]\n' filtered_variants.vcf.gz > ${cell}_filtered.tsv
    """

    stub:
    """
    touch ${cell}_filtered.tsv
    """

}


//

process COLLAPSE_FREEBAYES {

    tag "${sample_name}"
    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), 
        val(cells), 
        path(calls)

    output:
    tuple val(sample_name),
        path('allele_table.csv.gz'), emit: allele_counts

    script:
    """ 
    python ${baseDir}/bin/benchmark/collapse_bulk_methods.py
    """

    stub:
    """
    touch allele_table.csv.gz
    """

}
