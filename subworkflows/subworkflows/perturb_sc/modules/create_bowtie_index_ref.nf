// BOWTIE_INDEX_REF module

nextflow.enable.dsl = 2

//

process BOWTIE_INDEX_REF {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC_reference_fa)

  output:
  tuple val(sample_name), path("custom_ref_genome"), emit: index

  script:
  """
  mkdir -p custom_ref_genome

  bowtie2-build \
  -f ${GBC_reference_fa} \
  custom_ref_genome/GBC_reference
  """

  stub:
  """
  echo ${sample_name} > sample
  touch custom_ref_genome
  """

}