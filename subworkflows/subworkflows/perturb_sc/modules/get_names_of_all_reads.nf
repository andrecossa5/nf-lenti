// GET_READS_NAMES module

nextflow.enable.dsl = 2

//

process GET_READS_NAMES {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), path(R1_raw)

  output:
  tuple val(sample_name), path("all_reads.txt"), emit: names

  script:
  """
  zcat ${R1_raw} \
  | awk 'NR%4==1 {print substr(\$1, 2)}' \
  > reads_all_unsorted.txt

  LC_ALL=C sort --parallel=${task.cpus} reads_all_unsorted.txt > all_reads.txt

  rm reads_all_unsorted.txt
  """

  stub:
  """
  echo ${sample_name} > sample
  touch all_reads.txt
  """

}