// PLOT_QC_BULK module

nextflow.enable.dsl = 2

//

process PLOT_QC_BULK {
    tag { sample_name }
    label 'scLT'

    publishDir "${params.outdir}/${sample_name}/qc", mode: 'copy', overwrite: true

    input:
    tuple val(sample_name), path(raw_counts), path(corrected_counts), path(correction_df)

    output:
    tuple val(sample_name), path('before_vs_after_correction.png'), emit: before_after_plot
    tuple val(sample_name), path('corrected_counts_distribution.png'), emit: corrected_dist_plot
    tuple val(sample_name), path('raw_counts_distribution.png'), emit: raw_dist_plot

    script:
    """
    python ${baseDir}/bin/bulk_gbc/qc_plots_bulk.py \
        --raw_counts ${raw_counts} \
        --corrected_counts ${corrected_counts} \
        --correction_df ${correction_df} \
        --outdir . \
        --n_reads ${params.bulk_gbc_min_n_reads}
    """
}
