// generate_run_summary_bulk module

nextflow.enable.dsl = 2

//

// Process
process generate_run_summary_bulk {

  label 'scLT'
  tag "${sample_name}"

  input:
  tuple val(sample_name), 
  path(raw_counts), 
  path(corrected_counts), 
  path(correction_df), 
  path(before_after_plot), 
  path(corrected_dist_plot), 
  path(raw_dist_plot)

  output:
  tuple val(sample_name), path('run_summary.json'), emit: summary

  script:
  """
  python \
  ${baseDir}/bin/bulk_gbc/create_run_summary.py \
  --indir ${params.raw_data_input} \
  --outdir . \
  --params_outdir ${params.outdir} \
  --anchor_sequence ${params.bulk_gbc_anchor_sequence} \
  --sample ${sample_name} \
  --raw_counts ${raw_counts} \
  --corrected_counts ${corrected_counts} \
  --correction_df ${correction_df} \
  --min_n_reads ${params.bulk_gbc_min_n_reads} \
  --hamming_treshold ${params.bulk_gbc_graph_clustering_hamming_treshold} \
  --qc_images ${before_after_plot} ${corrected_dist_plot} ${raw_dist_plot}
  """

  stub:
  """
  touch run_summary.json
  """

}
