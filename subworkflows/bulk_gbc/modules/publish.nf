// publish_bulk module

nextflow.enable.dsl = 2

//

process publish_bulk {

  tag "${sample_name}"

  // Publish
  publishDir "${params.bulk_gbc_outdir}/${sample_name}/", mode: 'copy'

  input:
  tuple val(sample_name), path(raw_counts)
  tuple val(sample_name), path(corrected_counts)
  tuple val(sample_name), path(correction_df)
  tuple val(sample_name), path(run_summary)

  output:
  path raw_counts
  path corrected_counts
  path correction_df
  path run_summary
  val sample_name, emit: finish_flag

  script:
  """
  echo "Moving all necessary files to ${params.bulk_gbc_outdir}/${sample_name}/..."
  """

}
