// GET_NAMES_NOT_ALIGNED module

nextflow.enable.dsl = 2

//

process GET_NAMES_NOT_ALIGNED {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), path(all_reads), path(reads_aligned)

  output:
  tuple val(sample_name), path("reads_transcriptomic.txt"), emit: names

  script:
  """
  LC_ALL=C comm -23 \
  ${all_reads} \
  ${reads_aligned} \
  > reads_transcriptomic.txt
  """

  stub:
  """
  echo ${sample_name} > sample
  touch reads_transcriptomic.txt
  """

}