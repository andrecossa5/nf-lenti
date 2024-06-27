// MERGE module

nextflow.enable.dsl = 2

//

process MERGE_BAM {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam_tenx), path(bam_maester)

  output:
  tuple val(sample_name), path("merged_mitobam.bam"), emit: bam

  script:
  """
  samtools index -@ ${task.cpus} ${bam_tenx}
  samtools index -@ ${task.cpus} ${bam_maester}
  samtools merge -@ ${task.cpus} -o merged_mitobam.bam ${bam_tenx} ${bam_maester}
  """

  stub:
  """
  touch merged_mitobam.bam
  """

} 