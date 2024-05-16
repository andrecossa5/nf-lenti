// generate_run_summary_sc

nextflow.enable.dsl = 2

//

process generate_run_summary_sc {

    tag "${sample_name}"
 
    input:
    tuple val(sample_name), 
        path(R1), 
        path(GBCs), 
        path(filtered), 
        path(cells_summary),
        path(clones_summary)
  
    output:
    tuple val(sample_name), path("run_summary.txt"), emit: summary

    script:
    """
    echo "Summary Step 2, sample ${sample_name}" > run_summary.txt
    echo "-------------------------------------" >> run_summary.txt
    echo "" >> run_summary.txt
    echo "Overview" >> run_summary.txt
    echo "- Date of analysis:  \$(date)" >> run_summary.txt
    echo "- User:              ${USER}" >> run_summary.txt
    echo "- Working directory: ${PWD}" >> run_summary.txt
    echo "" >> run_summary.txt
    echo "Parameters" >> run_summary.txt
    echo "--sc_tenx_indir:                         ${params.sc_tenx_indir}" >> run_summary.txt
    echo "--sc_gbc_indir:                         ${params.sc_gbc_indir}" >> run_summary.txt
    echo "--sc_outdir:                        ${params.sc_outdir}" >> run_summary.txt
    echo "--pattern:                          ${params.sc_gbc_anchor_sequence}" >> run_summary.txt
    echo "--ref:                              ${params.ref}" >> run_summary.txt
    echo "Numbers" >> run_summary.txt
    echo "- Total reads:            \$(zcat ${R1} | awk 'END{print NR/4}' | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC in bulk reference:         \$(cat ${params.bulk_gbc_outdir}/${sample_name}/GBC_counts_corrected.csv| wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC found in sc: \$(cat ${GBCs} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Putative cell n (Solo cell-calling): \$(zcat ${filtered}/barcodes.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Total number of transcripts:     \$(zcat ${filtered}/features.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- n clones:     \$(cat ${clones_summary} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- n cells confidently assigned to GBC clones: \$(cat ${cells_summary} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    """

}