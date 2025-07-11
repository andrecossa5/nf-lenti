// publish_bulk module

nextflow.enable.dsl = 2

//

process publish_bulk {

  label 'scLT'
  tag "${sample_name}"

  // Publish
  publishDir "${params.outdir}/${sample_name}/", mode: 'copy'

  input:
  tuple val(sample_name), path(raw_counts), path(corrected_counts), path(correction_df), path(run_summary)

  output:
  path raw_counts
  path corrected_counts
  path correction_df
  path run_summary
  val sample_name, emit: finish_flag

  script:
  """
  echo "Moving all necessary files to ${params.outdir}/${sample_name}/..."
  """

}
