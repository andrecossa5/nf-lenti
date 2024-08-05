// CONSENSUS_NANOPORE module

nextflow.enable.dsl = 2

//

process CONSENSUS_NANOPORE {

  tag "${cell}"

  input:
  tuple val(cell), path(bam)

  output:
  tuple val("prova"), val(cell), path("${cell}_consensus_stats.csv"), path("${cell}_allelic_table.csv"), emit: allelic_tables

  script:
  """
  # Create allelic tables
  python ${baseDir}/bin/test_fgbio/pileup_nanopore.py ${cell} ${bam} ${params.scm_seq_positions}
  """

  stub:
  """
  touch ${cell}_consensus_stats.csv
  touch ${cell}_allelic_table.csv
  """

}


