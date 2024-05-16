// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(cell_folder)

  output:
  tuple val(sample_name), val(cell), path('filtered_consensus.bam'), emit: filtered_consensus_bam

  script:
  """
  touch filtered_consensus.bam
  # fgbio stuff + arguments --> 
  # cell_bams: un folder con tutti i bam delle singole cellule. Solo reads pulite e ben
  # allineate e clippate.
  """

  stub:
  """
  touch filtered_consensus.bam
  """

}
