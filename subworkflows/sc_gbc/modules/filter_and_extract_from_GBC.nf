// GET_GBC_ELEMENTs module

nextflow.enable.dsl = 2

//

process GET_GBC_ELEMENTS {

  tag "${sample_name}"
   
  input:
  tuple val(sample_name), path(lentibam), path(filtered)

  output:
  tuple val(sample_name),  path('GBC_read_elements.tsv'), emit: elements

  script:
  """
  python ${baseDir}/bin/sc_gbc/get_CBC_GBC_UMI.py \
  ${lentibam} ${filtered} ${params.sc_gbc_anchor_sequence} ${params.sc_gbc_anchor_treshold}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_read_elements.tsv
  """

}

