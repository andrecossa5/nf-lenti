// sc subworkflow

nextflow.enable.dsl = 2

// Include here
include { MERGE_R1 } from "./modules/merge_R1.nf"
include { MERGE_R2 } from "./modules/merge_R2.nf"
include { EXTRACT_R2 } from "./modules/extract_first_33_R2.nf"
include { BOWTIE_INDEX_GBC_PATTERN } from "./modules/create_bowtie_index_pattern.nf"
include { BOWTIE_INDEX_REF } from "./modules/create_bowtie_index_ref.nf"
include { ALIGN_33_R2 } from "./modules/align_first_33_R2.nf"
include { GET_READS_NAMES } from "./modules/get_names_of_all_reads.nf"
include { GET_NAMES_ALIGNED } from "./modules/get_names_aligned.nf"
include { GET_NAMES_NOT_ALIGNED } from "./modules/get_names_not_aligned.nf"
include { PREP_GBC } from "./modules/prep_GBC.nf"
include { PREP_TRANSCRIPT } from "./modules/prep_transcript.nf"
include { SOLO } from "./modules/Solo.nf"
include { FASTA_FROM_REF } from "./modules/fasta_from_ref.nf"
include { GET_GBC_ELEMENTS } from "./modules/filter_and_extract_from_GBC.nf"
include { GBC_TO_FASTA } from "./modules/gbc_to_fasta.nf"
include { ALIGN_GBC } from "./modules/align_GBC_to_ref.nf"
include { CELL_ASSIGNMENT } from "./modules/cell_assignment.nf"

//

process generate_run_summary_sc {

    tag "${sample_name}"
 
    input:
    tuple val(sample_name), path(all_reads), path(reads_transcript), path(reads_aligned), path(GBCs), path(filtered)
  
    output:
    tuple val(sample_name), path("run_summary.txt"), emit: run_summary

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
    echo "--indir:                ${params.perturb_sc_indir}" >> run_summary.txt
    echo "--outdir:               ${params.perturb_sc_outdir}" >> run_summary.txt
    echo "${sample_name} specific I/O: ${params.perturb_sc_indir}/${sample_name}, ${params.perturb_sc_outdir}/${sample_name}" >> run_summary.txt
    echo "--step_1_out:           ${params.perturb_bulk_out}" >> run_summary.txt
    echo "--pattern:              ${params.perturb_sc_pattern}" >> run_summary.txt
    echo "--ref:                  ${params.ref}" >> run_summary.txt
    echo "Numbers" >> run_summary.txt
    echo "- Reads in input:                  \$(cat ${all_reads} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Transcriptomic reads:            \$(cat ${reads_transcript} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- GBC-containing reads:            \$(cat ${reads_aligned} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC in reference:         \$(cat ${params.perturb_bulk_out}/${sample_name}/read_count_by_GBC_corrected.tsv | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC found in this sample: \$(cat ${GBCs} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Putative cell n (Solo cell-calling): \$(zcat ${filtered}/barcodes.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Total number of transcripts:     \$(zcat ${filtered}/features.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    """

}

//

process publish_sc {

    tag "${sample_name}"
    publishDir "${params.perturb_sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), 
          path(CBC_GBC), 
          path(CBC_GBC_plot), 
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

    script:
    """
    echo "Moving all output files to ${params.perturb_sc_outdir}/${sample_name}/..."
    """

}

//----------------------------------------------------------------------------//
// perturb_sc subworkflow
//----------------------------------------------------------------------------//

workflow perturb_sc {
    
    take:
        ch_input

    main:

        // Prep all reads, independently
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        EXTRACT_R2(MERGE_R2.out.R2)
        BOWTIE_INDEX_GBC_PATTERN()
        ALIGN_33_R2(EXTRACT_R2.out.first_33, BOWTIE_INDEX_GBC_PATTERN.out.index)
        GET_READS_NAMES(MERGE_R1.out.R1)
        GET_NAMES_ALIGNED(ALIGN_33_R2.out.R2_aligned)

        // Get (for each sample) raw fastqs, read and aligned read names. Separate fastqs
        GET_NAMES_NOT_ALIGNED(GET_READS_NAMES.out.names.combine(GET_NAMES_ALIGNED.out.names, by:0))
        PREP_GBC(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0).combine(GET_NAMES_ALIGNED.out.names, by:0))
        PREP_TRANSCRIPT(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0).combine(GET_NAMES_NOT_ALIGNED.out.names, by:0))
        
        // STARSolo
        SOLO(PREP_TRANSCRIPT.out.reads)

        // Assign cells to clones
        FASTA_FROM_REF(ch_input)
        BOWTIE_INDEX_REF(FASTA_FROM_REF.out.fasta)
        GET_GBC_ELEMENTS(PREP_GBC.out.reads.combine(SOLO.out.filtered, by:0))
        GBC_TO_FASTA(GET_GBC_ELEMENTS.out.elements.map{ it -> tuple(it[0], it[3]) })
        ALIGN_GBC(BOWTIE_INDEX_REF.out.index.combine(GBC_TO_FASTA.out.fasta, by:0))
        CELL_ASSIGNMENT(GET_GBC_ELEMENTS.out.elements.combine(ALIGN_GBC.out.names, by:0))

        // Summary
        summary_input = GET_READS_NAMES.out.names
            .combine(GET_NAMES_NOT_ALIGNED.out.names, by:0)
            .combine(GET_NAMES_ALIGNED.out.names, by:0)
            .combine(GET_GBC_ELEMENTS.out.elements.map{ it -> tuple(it[0], it[3]) }, by:0)
            .combine(SOLO.out.filtered, by:0)
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
            .combine(generate_run_summary_sc.out.run_summary, by:0)
        publish_sc(publish_input)

    emit:
        filtered = SOLO.out.filtered
        bam = SOLO.out.bam
        

}

//----------------------------------------------------------------------------//