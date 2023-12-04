// INFER_PREVALENCES module

nextflow.enable.dsl = 2

//

// Infere GBC clones frequencies from their read counts using ground truth spikins data
process INFER_PREVALENCES {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(counts)

  output:
  tuple val(sample_name), path('df_spikeins.csv'), emit: df_spikeins
  tuple val(sample_name), path('spikeins_fit.png'), emit: spikeins_plot
  tuple val(sample_name), path('clonal_prevalences.csv'), emit: stats_table
  tuple val(sample_name), path('final_prevalences.png'), emit: prevalences_plot
  tuple val(sample_name), path('GBC_1read_distributions.png'), emit: distributions_1read
  tuple val(sample_name), path("GBC_morereads_distributions.png"), emit: distributions_morereads

  script:
  """
  python \
  ${baseDir}/bin/bulk_gbc/infer_clone_prevalences.py \
  -i ${counts} \
  --path_spikeins_table ${params.bulk_gbc_spikeins_table} \
  --n_reads ${params.bulk_gbc_n_reads}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch clonal_prevalences.csv
  touch df_spikeins.csv
  touch spikeins_fit.png
  touch final_prevalences.png
  touch GBC_1read_distributions.png
  touch GBC_morereads_distributions.png
  """

}
