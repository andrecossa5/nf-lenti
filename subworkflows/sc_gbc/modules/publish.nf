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
          path(occurrences),
          path(cells_summary), 
          path(clones_summary), 
          path(bam), 
          path(stats),
          path(summary), 
          path(filtered), 
          path(raw), 
          path(run_summary)

    output:
    path raw
    path CBC_GBC
    path CBC_GBC_plot
    path cells_summary
    path clones_summary
    path bam
    path stats
    path summary
    path filtered
    path run_summary
    path occurrences

    script:
    """
    echo "Moving all output files to ${params.sc_outdir}/${sample_name}/..."
    """

}
