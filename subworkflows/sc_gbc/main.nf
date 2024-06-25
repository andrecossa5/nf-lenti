// sc subworkflow
nextflow.enable.dsl = 2

// Include here
include { MERGE_R1 } from "../maester/modules/merge_R1.nf"
include { MERGE_R2 } from "../maester/modules/merge_R2.nf"
include { SOLO } from "./modules/Solo.nf"
include { FILTER_LENTIBAM } from "./modules/filter_lentibam.nf"
include { SPLIT_BAM } from "./modules/split_bam.nf"
include { EXTRACT_FASTA } from "./modules/extract_fasta.nf"
include { CONSENSUS_LENTI } from "./modules/consensus_lenti.nf"
include { COLLAPSE_TSV } from "./modules/collapse_tsv.nf"
include { CELL_ASSIGNMENT } from "./modules/cell_assignment.nf"
include { generate_run_summary_sc } from "./modules/run_summary.nf"
include { publish_sc_gbc } from "./modules/publish.nf"

 
//


workflow sc_gbc {
    
    take:
        ch_sc_gbc
        ch_filtered

    main:
 
        // Merge reads and alignment
        MERGE_R1(ch_sc_gbc)
        MERGE_R2(ch_sc_gbc)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
 
        // Merge mitobams and split into cell-specific bams
        FILTER_LENTIBAM(SOLO.out.bam.combine(ch_filtered, by:0))
        SPLIT_BAM(FILTER_LENTIBAM.out.filtered_lentibam)
        ch_cell_bams = SPLIT_BAM.out.cell_bams
            .map { it ->
                def sample = it[0]
                def paths = it[1]      
                return paths.collect { cell_path ->
                    def path_splitted = cell_path.toString().split('/')
                    def cell = path_splitted[-1].toString().split('\\.')[0]
                    return [sample, cell, cell_path]
                }
            } 
            .flatMap { it } 
        EXTRACT_FASTA(params.string_lentiviral)

        // Create consensus reads and cell assignment
        CONSENSUS_LENTI(ch_cell_bams, EXTRACT_FASTA.out.fasta)
        ch_collapse = CONSENSUS_LENTI.out.consensus_filtered_tsv 
            .map { it -> tuple(it[0], it[2]) }
            .groupTuple(by: 0)
        COLLAPSE_TSV(ch_collapse)
        CELL_ASSIGNMENT(COLLAPSE_TSV.out.elements)

        // Summary
        summary_input = MERGE_R1.out.R1
            .combine(COLLAPSE_TSV.out.elements, by:0)
            .combine(ch_filtered, by:0)
            .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
            .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
        generate_run_summary_sc(summary_input)

        // Publishing
        publish_ch = CELL_ASSIGNMENT.out.CBC_GBC_combos 
            .combine(CELL_ASSIGNMENT.out.combo_plot, by:0)
            .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
            .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
            .combine(CELL_ASSIGNMENT.out.summary, by:0)
            .combine(generate_run_summary_sc.out.summary, by:0)
        publish_sc_gbc(publish_ch)

    emit:
        summary = generate_run_summary_sc.out.summary

}