// collapse_output module

nextflow.enable.dsl = 2

//

process collapse_output {

  label 'scLT'

  // Publish
  publishDir "${params.outdir}", mode: 'copy'

  input:
  val last

  output:
  path summary

  script:
  """
  python ${baseDir}/bin/bulk_gbc/collapse_outputs.py -i ${params.outdir}
  """

  stub:
  """
  echo "Collapsing output in ${params.outdir}..."
  mkdir summary
  """

}
