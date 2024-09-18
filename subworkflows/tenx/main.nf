// sc_pp workflow

// Include here
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"

// 

process publish_tenx {

    tag "${sample_name}"

    // Publish
    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name),
          path (raw),
          path (filtered),
          path (stats), 
          path (summary),
          path (bam)

    output:
    path raw
    path filtered
    path stats
    path summary

    script:
    """
    echo moving everything to ${params.sc_outdir}
    """

}

// 


//----------------------------------------------------------------------------//
// tenx subworkflow
//----------------------------------------------------------------------------//


workflow tenx {

    take:
        ch_input

    main:
    
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
        publish_input = SOLO.out.raw
            .combine(SOLO.out.filtered, by:0)
            .combine(SOLO.out.stats, by:0)
            .combine(SOLO.out.summary, by:0)
            .combine(SOLO.out.bam, by:0)
        publish_tenx(publish_input)

    emit:
    
        filtered = SOLO.out.filtered
        bam = SOLO.out.bam

}


//----------------------------------------------------------------------------//