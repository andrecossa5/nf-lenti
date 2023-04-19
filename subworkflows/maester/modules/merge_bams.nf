// MERGE module

nextflow.enable.dsl = 2

//

process MERGE {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam_1), path(bam_2)

  output:
  tuple val(sample_name), path("merged_mitobam.bam"), emit: mitobam

  script:
  """
  samtools index -@ ${task.cpus} ${bam_1}
  samtools index -@ ${task.cpus} ${bam_2}
  samtools merge -@ ${task.cpus} ./merged_mitobam.bam ${bam_1} ${bam_2}
  """

  stub:
  """
  touch merged_mitobam.bam
  """

} 