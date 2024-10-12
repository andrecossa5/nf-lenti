// MERGE module

nextflow.enable.dsl = 2

//

process INDEX_AND_MERGE {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bams)

  output:
  tuple val(sample_name), path("merged_bam.bam"), emit: bam

  script:
  """
  samtools index -@ ${task.cpus} -M ${bams}
  samtools merge -@ ${task.cpus} -o merged_bam.bam ${bams}
  """

  stub:
  """
  touch merged_bam.bam
  """

} 