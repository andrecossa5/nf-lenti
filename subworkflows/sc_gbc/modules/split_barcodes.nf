// SPLIT_BARCODES module

nextflow.enable.dsl = 2

//

process SPLIT_BARCODES {

  tag "${sample_name}" 

  input:
  tuple val(sample_name), path(filtered)

  output:
  tuple val(sample_name), path('barcodes_*.csv'), emit: barcodes 

  script:
  """
  python ${baseDir}/bin/sc_gbc/split_barcodes.py ${filtered}/barcodes.tsv.gz ${params.CBs_chunk_size}
  """

  stub:
  """
  touch barcodes_1.csv barcodes_2.csv barcodes_3.csv
  """

}
