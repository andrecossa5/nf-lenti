// GET_NAMES_ALIGNED module

nextflow.enable.dsl = 2

//

process GET_NAMES_ALIGNED {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(R2_first_33nt_aligned)

  output:
  tuple val(sample_name), path("reads_aligned.txt"), emit: names

  script:
  """
  zcat ${R2_first_33nt_aligned} \
  | grep -P '\tconstruct' \
  | cut -f1 \
  > reads_aligned_unsorted.txt

  LC_ALL=C sort --parallel=4 reads_aligned_unsorted.txt > reads_aligned.txt

  rm reads_aligned_unsorted.txt
  """

  stub:
  """
  echo ${sample_name} > sample
  touch reads_aligned.txt
  """

}
