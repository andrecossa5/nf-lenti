// Bulk workflow

// Include here
nextflow.enable.dsl = 2
include { SEARCH_PATTERNS } from "./modules/generate_search_patterns.nf"
include { EXTRACT_GBC } from "./modules/extract_GBC.nf"
include { CORRECT_AND_COUNT } from "./modules/correct_and_count.nf"
include { INFER_PREVALENCES } from "./modules/infer_clone_prevalences.nf"
include { generate_run_summary_bulk } from "./modules/run_summary.nf"
include { publish_bulk } from "./modules/publish.nf"
include { collapse_output } from "./modules/collapse_out.nf"

 
//


workflow bulk_gbc {

  take:
      ch_input
 
  main:
      SEARCH_PATTERNS()
      EXTRACT_GBC(ch_input, SEARCH_PATTERNS.out.search_patterns)
      CORRECT_AND_COUNT(EXTRACT_GBC.out.GBC)
  
      // Summary and cleanup 
      generate_run_summary_bulk(
        CORRECT_AND_COUNT.out.raw_counts,
        CORRECT_AND_COUNT.out.corrected_counts,
        CORRECT_AND_COUNT.out.correction_df,
      )
      publish_bulk(
        CORRECT_AND_COUNT.out.raw_counts,
        CORRECT_AND_COUNT.out.corrected_counts,
        CORRECT_AND_COUNT.out.correction_df,
        generate_run_summary_bulk.out.summary
      )
      collapse_output(
        publish_bulk.out.finish_flag.collect().last() // Need to fire only at the end
      )

  emit:
      flags = publish_bulk.out.finish_flag.collect()
}