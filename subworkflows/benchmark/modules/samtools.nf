// samtools modules

nextflow.enable.dsl = 2

//

process SAMTOOLS {

    tag "${sample_name}: ${cell}"

    input:
    tuple val(sample_name), val(cell), path(bam)
    path(reference)

    output:
    tuple val(sample_name), val(cell), path("${cell}_filtered.tsv"), emit: calls

    script:
    """ 
    picard MarkDuplicates I=${bam} O=cell_dedup.bam M=deduplication_metrics.txt REMOVE_DUPLICATES=true
    bcftools mpileup -f ${reference} -a FORMAT/AD,FORMAT/DP -Q 30 -Ou cell_dedup.bam | bcftools call -mv -Oz -o cell.vcf.gz
    bcftools filter -i 'QUAL>20' -Oz -o filtered.vcf.gz cell.vcf.gz
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%AD]\t[%DP]\n' filtered.vcf.gz > ${cell}_filtered.tsv
    """

    stub:
    """
    touch ${cell}_filtered.tsv
    """

}

//

process COLLAPSE_SAMTOOLS {

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
