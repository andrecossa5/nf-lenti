// PREP_GBC module

nextflow.enable.dsl = 2

//

process PREP_GBC {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(R1_raw), path(R2_raw), path(reads_aligned)

  output:
  tuple val(sample_name), path("aligned_R1.fq.gz"), path("aligned_R2.fq.gz"), emit: reads

  script:
  """
  seqtk subseq \
  ${R1_raw} \
  ${reads_aligned} \
  | pigz --fast -p ${task.cpus} \
  > aligned_R1.fq.gz

  seqtk subseq \
  ${R2_raw} \
  ${reads_aligned} \
  | pigz --fast -p ${task.cpus} \
  > aligned_R2.fq.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch aligned_R1.fq.gz
  touch aligned_R2.fq.gz
  """

}