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
  tuple val(sample_name), path("CBC_GBC_UMI_read_distribution.png"), emit: combo_dist
  tuple val(sample_name), path("counts.pickle"), emit: counts
  tuple val(sample_name), path("selected_UMIs.png"), emit: selected_umi_plot

  script: 
  """
  python ${baseDir}/bin/sc_gbc/cell_assignment.py \
  --sample ${sample_name} \
  --path_bulk ${params.bulk_gbc_outdir} \
  --path_sc ${elements} \
  --sample_map ${params.gbc_sample_map} \
  --ncores ${task.cpus} \
  --bulk_correction_treshold ${params.bulk_correction_treshold} \
  --sc_correction_treshold ${params.sc_correction_treshold} \
  --filtering_method  ${params.umi_filtering_method} \
  --coverage_treshold ${params.coverage_treshold} \
  --correction_type ${params.correction_type} \
  --umi_treshold ${params.umi_treshold} \
  --p_treshold ${params.p_treshold} \
  --ratio_to_most_abundant_treshold ${params.ratio_to_most_abundant_treshold} \
  --sample_params ${params.sample_params}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch CBC_GBC_combos.tsv.gz
  touch counts.pickle
  touch clones_summary_table.csv
  touch cells_summary_table.csv
  touch CBC_GBC_combo_status.png
  touch clone_calling_summary.txt
  touch CBC_GBC_UMI_read_distribution.png
  touch selected_UMIs.png
  """

}