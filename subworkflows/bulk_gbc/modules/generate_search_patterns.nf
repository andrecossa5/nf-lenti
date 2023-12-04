// SEARCH_PATTERNS module

nextflow.enable.dsl = 2

//

process SEARCH_PATTERNS {

  output:
  path "search_patterns.tsv", emit: search_patterns

  script:
  """
  python \
  ${baseDir}/bin/bulk_gbc/generate_search_patterns.py ${params.bulk_gbc_anchor_sequence}
  """

  stub:
  """
  touch search_patterns.tsv
  """

}
