// sc subworkflow

nextflow.enable.dsl = 2

// Include here
include { MERGE_FQ } from "./modules/merge_fq.nf"
include { SOLO } from "./modules/Solo.nf"
include { GET_GBC_ELEMENTS } from "./modules/filter_and_extract_from_GBC.nf"
include { GET_LENTIBAM } from "./modules/get_lentibam.nf"
include { CELL_ASSIGNMENT } from "./modules/cell_assignment.nf"
include { generate_run_summary_sc } from "./modules/run_summary.nf"
include { publish_sc } from "./modules/publish.nf"

 
//

workflow sc_gbc {
    
    take:
        ch_sc

    main:
 
        // Merge reads
        MERGE_FQ(ch_sc)

        // STARSolo
        SOLO(MERGE_FQ.out.reads)
 
        // Assign cells to clones
        GET_LENTIBAM(SOLO.out.bam)
        GET_GBC_ELEMENTS(GET_LENTIBAM.out.lentibam.combine(SOLO.out.filtered, by:0))
        CELL_ASSIGNMENT(GET_GBC_ELEMENTS.out.elements)

        // Summary
        summary_input = MERGE_FQ.out.reads.map{ it -> tuple(it[0], it[1]) }
            .combine(GET_GBC_ELEMENTS.out.elements, by:0)
            .combine(SOLO.out.filtered, by:0)
            .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
            .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
        generate_run_summary_sc(summary_input)

        // Publishing
        publish_input = CELL_ASSIGNMENT.out.CBC_GBC_combos
            .combine(CELL_ASSIGNMENT.out.plot, by:0)
            .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
            .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
            .combine(SOLO.out.bam, by:0)
            .combine(SOLO.out.stats, by:0)
            .combine(SOLO.out.summary, by:0)
            .combine(SOLO.out.filtered, by:0)
            .combine(SOLO.out.raw, by:0)
            .combine(generate_run_summary_sc.out.summary, by:0)
        publish_sc(publish_input)

    emit:
        summary = generate_run_summary_sc.out.summary
        filtered = SOLO.out.filtered
        bam = SOLO.out.bam

}