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
  python ${baseDir}/bin/test_fgbio/pileup_nanopore.py \
  ${cell} ${bam} ${params.scm_seq_positions} ${params.scm_seq_min_reads} ${params.scm_seq_base_consensus_error_thr} ${params.scm_seq_base_quality_thr}
  """

  stub:
  """
  touch ${cell}_consensus_stats.csv
  touch ${cell}_allelic_table.csv
  """

}


