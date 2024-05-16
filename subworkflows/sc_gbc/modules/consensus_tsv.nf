// CONSENSUS_TSV module

nextflow.enable.dsl = 2

//

process CONSENSUS_TSV {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(consensus_bam)

  output:
  tuple val(sample_name), path('filtered_consensus.tsv'), emit: filtered_consensus_tsv

  script:
  """
  touch filtered_consensus.tsv
  # fgbio stuff + arguments --> 
  # cell_bams: un folder con tutti i bam delle singole cellule. Solo reads pulite e ben
  # allineate e clippate.
  """

  stub:
  """
  touch filtered_consensus.tsv
  """

}
