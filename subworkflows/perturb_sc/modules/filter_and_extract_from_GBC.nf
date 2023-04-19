// GET_GBC_ELEMENTs module

nextflow.enable.dsl = 2

//

process GET_GBC_ELEMENTS {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), path(R1), path(R2), path(filtered)

  output:
  tuple val(sample_name),  path('filtered_R1.fq.gz'), path('filtered_R2.fq.gz'), emit: reads
  tuple val(sample_name),  path('CBCs_by_read.tsv'), path('UMIs_by_read.tsv'), path('GBCs_by_read.tsv'), emit: elements

  script:
  """
  python ${baseDir}/bin/filter_and_extact_from_GBC_reads.py ${R1} ${R2} ${filtered}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch filtered_R1.fq.gz
  touch filtered_R2.fq.gz
  touch CBCs_by_read.tsv
  touch UMIs_by_read.tsv
  touch GBCs_by_read.tsv
  """

}

