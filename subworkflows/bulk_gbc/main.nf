// Bulk workflow

nextflow.enable.dsl = 2

// Module includes
include { SEARCH_PATTERNS }           from "./modules/generate_search_patterns.nf"
include { EXTRACT_GBC }               from "./modules/extract_GBC.nf"
include { CORRECT_AND_COUNT }         from "./modules/correct_and_count.nf"
include { generate_run_summary_bulk } from "./modules/run_summary.nf"
include { publish_bulk }              from "./modules/publish.nf"
include { collapse_output }           from "./modules/collapse_out.nf"
include { PLOT_QC_BULK }               from "./modules/plot_qc_bulk.nf"

workflow bulk_gbc {
    take:
        ch_input

    main:
        SEARCH_PATTERNS()
        EXTRACT_GBC(ch_input, SEARCH_PATTERNS.out.search_patterns)
        CORRECT_AND_COUNT(EXTRACT_GBC.out.GBC)
        PLOT_QC_BULK(CORRECT_AND_COUNT.out.gbc_counts)
        
        
        summary_input = CORRECT_AND_COUNT.out.gbc_counts
            .combine(PLOT_QC_BULK.out.before_after_plot,         by: 0)
            .combine(PLOT_QC_BULK.out.corrected_dist_plot,       by: 0)
            .combine(PLOT_QC_BULK.out.raw_dist_plot,             by: 0)

        generate_run_summary_bulk(summary_input)
        ch_counts  = CORRECT_AND_COUNT.out.gbc_counts
        ch_summary = generate_run_summary_bulk.out.summary

        ch_paired = ch_counts.join(ch_summary, by: 0)
        publish_bulk(ch_paired)

        // collapse_output now creates summary/ + run_report.html
        collapse_output(publish_bulk.out.finish_flag.collect().last())

    emit:
        flags = publish_bulk.out.finish_flag.collect()

}