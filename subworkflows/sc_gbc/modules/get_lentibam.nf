// GET_LENTIBAM module

nextflow.enable.dsl = 2

//

process GET_LENTIBAM {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("lentibam.bam"), emit: lentibam

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  samtools view -b ${bam} lentiCassette > lentibam.bam
  samtools index -@ ${task.cpus} lentibam.bam
  """

  stub:
  """
  echo ${sample_name} > sample
  touch lentibam.bam
  touch lentibam.bam.bai
  """

}