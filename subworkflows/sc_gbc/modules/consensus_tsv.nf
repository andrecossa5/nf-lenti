// CONSENSUS_TSV module

nextflow.enable.dsl = 2

//

process CONSENSUS_TSV {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(consensus_filtered_bam)

  output:
  tuple val(sample_name), val(cell), path("consensus_filtered.bam"), emit: consensus_filtered_tsv

  script:
  """
  python ${baseDir}/bin/sc_gbc/consensus_tsv.py ${consensus_bam} ${cell}
  """

  stub:
  """
  touch "${cell}_consensus.tsv"
  """

}
