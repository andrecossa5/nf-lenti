// publish_sc module

nextflow.enable.dsl = 2

//

process publish_sc_gbc {

    label 'scLT'
    tag "${sample_name}"
    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), 
          path(CBC_GBC),
          path(CBC_GBC_plot),
          path(umi_dist),
          path(moi_dist),
          path(clone_sz),
          path(cells_summary), 
          path(clones_summary), 
          path(cell_assignment_summary),
          path(run_summary_json),
          path(umi_dist_interactive),
          path(moi_dist_interactive),
          path(clone_sz_interactive)
          //path(filter_summary)

    output:
    path CBC_GBC
    path CBC_GBC_plot
    path umi_dist
    path moi_dist
    path clone_sz
    path cells_summary
    path clones_summary
    path run_summary_json
    path cell_assignment_summary
    path umi_dist_interactive
    path moi_dist_interactive
    path clone_sz_interactive
    //path filter_summary

    script:
    """
    echo "Moving all output files to ${params.sc_outdir}/${sample_name}/..."
    """

}
