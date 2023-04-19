// INDEX module

nextflow.enable.dsl = 2

//

process INDEX {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("${bam}"), emit: bam
  tuple val(sample_name), path("${bam}.bai"), emit: index

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  """

  stub:
  """
  touch ${bam}
  touch ${bam}.bai
  """

}