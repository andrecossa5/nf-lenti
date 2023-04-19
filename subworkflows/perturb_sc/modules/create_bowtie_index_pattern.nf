// BOWTIE_INDEX_GBC_PATTERN module

nextflow.enable.dsl = 2

//

process BOWTIE_INDEX_GBC_PATTERN {

  output:
  path("GBC_pattern"), emit: index

  script:
  """
  mkdir -p GBC_pattern

  echo ">construct" > GBC_pattern.fa
  echo "${params.perturb_sc_pattern}" >> GBC_pattern.fa

  bowtie2-build \
  -f GBC_pattern.fa \
  GBC_pattern/GBC_pattern
  """

  stub:
  """
  touch GBC_pattern
  """

}
