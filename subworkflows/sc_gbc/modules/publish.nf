// publish_sc module

nextflow.enable.dsl = 2

//

process publish_sc_gbc {

    tag "${sample_name}"
    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), 
          path(CBC_GBC),
          path(CBC_GBC_plot),
          path(cells_summary), 
          path(clones_summary), 
          path(cell_assignment_summary),
          path(run_summary)
          //path(filter_summary)

    output:
    path CBC_GBC
    path CBC_GBC_plot
    path cells_summary
    path clones_summary
    path run_summary
    path cell_assignment_summary
    //path filter_summary

    script:
    """
    echo "Moving all output files to ${params.sc_outdir}/${sample_name}/..."
    """

}
