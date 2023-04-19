// ALIGN_33_R2 module

nextflow.enable.dsl = 2

//

process ALIGN_33_R2 {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(R2_first_33nt)
  path bowtie_index_GBC_pattern_folder

  output:
  tuple val(sample_name), path("R2_first_33_nt.fq3.gz"), emit: R2_aligned

  script:
  """
  bowtie2 \
  -N 1 \
  --norc \
  -x ${bowtie_index_GBC_pattern_folder}/GBC_pattern \
  -f ${R2_first_33nt} \
  | pigz --fast -p ${task.cpus} \
  > R2_first_33_nt.fq3.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch R2_first_33_nt.fq3.gz
  """

}