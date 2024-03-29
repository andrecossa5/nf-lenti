// publish_sc module

nextflow.enable.dsl = 2

//

process publish_sc {

    tag "${sample_name}"
    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), 
          path(CBC_GBC), 
          path(CBC_GBC_plot), 
          path(CBC_GBC_UMI_plot), 
          path(cells_summary), 
          path(clones_summary), 
          path(cell_assignment_summary),
          path(bam), 
          path(stats),
          path(alignment_summary), 
          path(filtered), 
          path(raw),
          path(run_summary),
          path(counts),
          path(selected_umi_plot)

    output:
    path CBC_GBC
    path CBC_GBC_plot
    path CBC_GBC_UMI_plot
    path cells_summary
    path clones_summary
    path bam
    path stats
    path raw
    path filtered
    path alignment_summary
    path run_summary
    path cell_assignment_summary
    path counts
    path selected_umi_plot 

    script:
    """
    echo "Moving all output files to ${params.sc_outdir}/${sample_name}/..."
    """

}
