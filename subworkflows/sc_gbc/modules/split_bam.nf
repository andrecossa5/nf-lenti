// SPLIT_BAM module

nextflow.enable.dsl = 2

//

process SPLIT_BAM {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam), path(filtered)

  output:
  tuple val(sample_name), path('cell_bams/*'), emit: cell_bams

  script:
  """
  # touch cell_bams
  # python ${baseDir}/bin/sc_gbc/split_bam.py <--args vari> --> cell_bams
  # cell_bams: un folder con tutti i bam delle singole cellule. Solo reads pulite e ben
  # allineate e clippate.
  """

  stub:
  """
  mkdir cell_bams
  cd cell_bams
  mkdir AAAA BBBB
  echo AAAA > AAAA/AAAA.txt
  echo AAAA > BBBB/BBBB.txt
  """

}
