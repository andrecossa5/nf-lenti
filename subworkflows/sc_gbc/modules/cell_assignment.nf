// CELL_ASSIGNMENT module

nextflow.enable.dsl = 2

//

process CELL_ASSIGNMENT {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), path(elements)

  output:
  tuple val(sample_name), path("CBC_GBC_combos.tsv.gz"), emit: CBC_GBC_combos
  tuple val(sample_name), path("clones_summary_table.csv"), emit: clones_summary
  tuple val(sample_name), path("cells_summary_table.csv"), emit: cells_summary
  tuple val(sample_name), path("CBC_GBC_combo_status.png"), emit: combo_plot
  tuple val(sample_name), path("clone_calling_summary.txt"), emit: summary

  script: 
  """
  python ${baseDir}/bin/sc_gbc/cell_assignment.py \
  --sample ${sample_name} \
  --path_sc ${elements} \
  --bulk_correction_treshold ${params.bulk_correction_treshold} \
  --umi_treshold ${params.umi_treshold} \
  --p_treshold ${params.p_treshold} \
  --max_ratio_treshold ${params.max_ratio_treshold} \
  --normalized_abundance_treshold ${params.normalized_abundance_treshold} \
  --sample_params ${params.sample_params}\
  
  --path_bulk ${params.bulk_gbc_outdir} \
  --sample_map ${params.gbc_sample_map} \
  """

  stub:
  """
  touch CBC_GBC_combos.tsv.gz
  touch clones_summary_table.csv
  touch cells_summary_table.csv
  touch CBC_GBC_combo_status.png
  touch clone_calling_summary.txt
  """

}