// collapse_output_sc module

nextflow.enable.dsl = 2

//

process collapse_output_sc {
   label 'scLT'
   publishDir "${params.sc_outdir}", mode: 'copy'

   input:
   path run_jsons

   output:
   path "sc_summary/run_report.html"
   path "sc_summary/run_summary.json"
   path "sc_summary/*.png"
   path "sc_summary/*_interactive.html"

   script:
   """
   mkdir -p sc_summary
   python ${baseDir}/bin/sc_gbc/collapse_outputs_sc.py \
        --input ${params.sc_outdir} \
        --template ${baseDir}/subworkflows/templates/report_template.html \
        --output sc_summary
   """
}
