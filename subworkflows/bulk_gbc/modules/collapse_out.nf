// collapse_output module

nextflow.enable.dsl = 2

//

process collapse_output {

  // Publish
  publishDir "${params.bulk_gbc_outdir}", mode: 'copy'

  input:
  val last

  output:
  path summary

  script:
  """
  python ${baseDir}/bin/bulk_gbc/collapse_outputs.py -i ${params.bulk_gbc_outdir}
  """

  stub:
  """
  echo "Collapsing output in ${params.bulk_gbc_outdir}..."
  mkdir summary
  """

}
