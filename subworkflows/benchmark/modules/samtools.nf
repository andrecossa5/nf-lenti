// samtools modules

nextflow.enable.dsl = 2

//

process SAMTOOLS {

    tag "${sample_name}: ${cell}"

    input:
    tuple val(sample_name), val(cell), path(bam)
    path(reference)

    output:
    tuple val(sample_name), path("${cell}_filtered.tsv"), emit: calls

    script:
    """ 
    picard MarkDuplicates I=${bam} O=cell_dedup.bam M=deduplication_metrics.txt REMOVE_DUPLICATES=true
    bcftools mpileup -f ${reference} -a FORMAT/AD,FORMAT/DP -Q 30 -Ou cell_dedup.bam | bcftools call -mv -Oz -o cell.vcf.gz
    bcftools filter  -i 'QUAL>20' -Oz -o filtered.vcf.gz cell.vcf.gz
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

    input:
    tuple val(sample_name), path(muts)

    output:
    tuple val(sample_name),
        path('AD.npz'), 
        path('DP.npz'),
        path('cells.txt'), emit: matrices

    script:
    """ 
    python ${baseDir}/bin/benchmark/collapse_samtools.py 
    """

    stub:
    """
    touch AD.npz
    touch DP.npz
    touch cells.txt
    """

}
