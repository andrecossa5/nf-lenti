// SPLIT_BARCODES module

nextflow.enable.dsl = 2

//

process SPLIT_BARCODES {

  tag "${sample_name}" 

  input:
  tuple val(sample_name), path(bam), path(barcodes)

  output:
  tuple val(sample_name), path(bam), path('barcodes_*.csv'), emit: barcodes 

  script:
  """
  python ${baseDir}/bin/common/split_barcodes.py ${barcodes} ${params.CBs_chunk_size}
  """

  stub:
  """
  touch barcodes_1.csv barcodes_2.csv barcodes_3.csv
  """

}
