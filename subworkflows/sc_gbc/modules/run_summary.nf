// generate_run_summary_sc module

nextflow.enable.dsl = 2

//

process generate_run_summary_sc {

    label 'scLT'
    tag "${sample_name}"

    input:
    tuple val(sample_name),
        path(GBCs),       // cells.tsv.gz
        path(cells_summary), //cells_summary_table.csv
        path(clones_summary), //clones_summary_table.csv
        path(clone_summary_txt), //clone_calling_summary.txt
        path(combo_plot), //CBC_GBC_combo_status.png
        path(umi_dist), //umi_distribution.png
        path(moi_dist), //moi_distribution.png
        path(clone_sz), //clone_size_distribution.png
        path(umi_dist_interactive),
        path(moi_dist_interactive),
        path(clone_sz_interactive)

    output:
    tuple val(sample_name), path("run_summary.json"), emit: summary_json

    script:
    """
    python ${baseDir}/bin/sc_gbc/make_run_summary_json.py \
        --sample ${sample_name} \
        --clone_summary_txt ${clone_summary_txt} \
        --cells_summary ${cells_summary} \
        --clones ${clones_summary} \
        --gbcs ${GBCs} \
        --bulk_gbc ${params.bulk_gbc_outdir}/${sample_name}/GBC_counts_corrected.csv \
        --combo_plot ${combo_plot} \
        --umi_dist ${umi_dist} \
        --moi_dist ${moi_dist} \
        --clone_sz ${clone_sz} \
        --umi_dist_interactive ${umi_dist_interactive} \
        --moi_dist_interactive ${moi_dist_interactive} \
        --clone_sz_interactive ${clone_sz_interactive} \
        --raw_data_input ${params.raw_data_input} \
        --raw_data_input_type ${params.raw_data_input_type} \
        --sc_outdir ${params.sc_outdir} \
        --pattern ${params.sc_gbc_anchor_sequence} \
        --ref ${params.ref} \
        --out_json run_summary.json
    """
}
