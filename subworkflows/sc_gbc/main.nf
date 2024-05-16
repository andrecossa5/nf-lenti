// sc subworkflow
nextflow.enable.dsl = 2

// Include here
include { MERGE_R1 } from "../maester/modules/merge_R1.nf"
include { MERGE_R2 } from "../maester/modules/merge_R2.nf"
include { SOLO } from "./modules/Solo.nf"
include { SPLIT_BAM } from "./modules/split_bam.nf"
include { CONSENSUS_BAM } from "./modules/consensus_bam.nf"
include { CONSENSUS_TSV } from "./modules/consensus_tsv.nf"
include { COLLAPSE_TSV } from "./modules/collapse_tsv.nf"
include { CELL_ASSIGNMENT } from "./modules/cell_assignment.nf"
include { generate_run_summary_sc } from "./modules/run_summary.nf"
include { publish_sc } from "./modules/publish.nf"

 
//


workflow sc_gbc {
    
    take:
        ch_sc_gbc
        ch_filtered

    main:
 
        // Merge reads and Solo
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
 
        // Consensus UMIs
        // SPLIT_BAM(SOLO.out.bam.combine(ch_filtered, by:0))
        // ch_cell_bams = Channel
        //     .fromPath("${SPLIT_BAM.out.cell_bams}/*", type:'dir')
        //     .map{ tuple(it.getName(), it) }
        // CONSENSUS_BAM(ch_cell_bams)

        // Cell assignment
        // CONSENSUS_TSV(CONSENSUS_BAM.out.filtered_consensus_bam)
        // COLLAPSE_TSV(CONSENSUS_TSV.out.filtered_consensus_tsv.map{ it[1] }.collect())
        // CELL_ASSIGNMENT(COLLAPSE_TSV.out.elements)

        // Summary
        // summary_input = MERGE_R1.out.R1.map{ it -> tuple(it[0], it[1]) }
        //     .combine(GET_GBC_ELEMENTS.out.elements, by:0)
        //     .combine(SOLO.out.filtered, by:0)
        //     .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
        //     .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
        // generate_run_summary_sc(summary_input)

        // Publishing
        // publish_input = CELL_ASSIGNMENT.out.CBC_GBC_combos
        //     .combine(CELL_ASSIGNMENT.out.combo_plot, by:0)
        //     .combine(CELL_ASSIGNMENT.out.combo_dist, by:0)
        //     .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
        //     .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
        //     .combine(CELL_ASSIGNMENT.out.summary, by:0)
        //     .combine(SOLO.out.bam, by:0)
        //     .combine(SOLO.out.stats, by:0)
        //     .combine(SOLO.out.summary, by:0)
        //     .combine(SOLO.out.filtered, by:0)
        //     .combine(SOLO.out.raw, by:0)
        //     .combine(generate_run_summary_sc.out.summary, by:0)
        //     .combine(CELL_ASSIGNMENT.out.counts, by:0)
        //     .combine(CELL_ASSIGNMENT.out.selected_umi_plot, by:0)
        // publish_sc(publish_input)

    emit:

        ch_test = SOLO.out.bam.combine(ch_filtered, by:0)
        // summary = generate_run_summary_sc.out.summary
        // filtered = SOLO.out.filtered
        // bam = SOLO.out.bam

}