// PREP_GBC module

nextflow.enable.dsl = 2

//

process PREP_TRANSCRIPT {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(R1_raw), path(R2_raw), path(reads_transcriptomic)

  output:
  tuple val(sample_name), path("transcriptomic_R1.fq.gz"), path("transcriptomic_R2.fq.gz"), emit: reads

  script:
  """
  seqtk subseq \
  ${R1_raw} \
  ${reads_transcriptomic} \
  | pigz --fast -p ${task.cpus} \
  > transcriptomic_R1.fq.gz

  seqtk subseq \
  ${R2_raw} \
  ${reads_transcriptomic} \
  | pigz --fast -p ${task.cpus} \
  > transcriptomic_R2.fq.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch transcriptomic_R1.fq.gz
  touch transcriptomic_R2.fq.gz
  """

}