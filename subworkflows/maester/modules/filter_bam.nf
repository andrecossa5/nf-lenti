// FILTER BAM modules

nextflow.enable.dsl = 2

//

process FILTER_10X_BAM {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("mitobam_I.bam"), emit: bam

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  samtools view ${bam} -b -@ ${task.cpus} chrM > mitobam_I.bam
  """

  stub:
  """
  touch mitobam_I.bam
  """

}

//

process FILTER_MAESTER_BAM {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("mitobam_II.bam"), emit: bam

  script:
  """
  samtools sort -l 1 -@ ${task.cpus} -o sorted.bam ${bam}
  samtools index -@ ${task.cpus} sorted.bam
  samtools view sorted.bam -b -@ ${task.cpus} chrM > mitobam_II.bam
  """

  stub:
  """
  touch mitobam_II.bam
  """

}