// MAEGATK module

nextflow.enable.dsl = 2

//

process MAEGATK {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("tables"), emit: tables
 
  script:
  """
  # ...
  """

  stub:
  """
  mkdir tables
  """

}

//

process COLLAPSE_MAEGATK {

  tag "${sample_name}"

  input:
  tuple val(sample_name), val(tables)

  output:
  tuple val(sample_name), path("final"), emit: matrices
 
  script:
  """
  # ...
  """

  stub:
  """
  mkdir final
  """

}