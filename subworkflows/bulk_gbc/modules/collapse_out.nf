// collapse_output module

nextflow.enable.dsl = 2

//

process collapse_output {
    label 'scLT'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val last

    output:
    path "summary", emit: summary

    script:
    """
    python ${baseDir}/bin/bulk_gbc/collapse_outputs.py \
        -i ${params.outdir} \
        --template ${baseDir}/subworkflows/templates/report_template.html
    """
}


