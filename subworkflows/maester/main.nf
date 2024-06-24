// measter subworkflow

// Include here
nextflow.enable.dsl = 2
include { MERGE_R1 } from "./modules/merge_R1.nf"
include { MERGE_R2 } from "./modules/merge_R2.nf"
include { ASSEMBLE_FQ } from "./modules/assemble_fastq.nf"
include { STAR } from "./modules/STAR.nf"
include { FILTER_10X_BAM } from "./modules/filter_bam.nf"
include { FILTER_MAESTER_BAM } from "./modules/filter_bam.nf"
include { FIX_TAGS } from "./modules/fix_tags.nf"
include { MERGE_BAM } from "./modules/merge_bams.nf"
include { FILTER_MITOBAM } from "./modules/filter_mitobam.nf"
include { SPLIT_BAM } from "../sc_gbc/modules/split_bam.nf"
include { EXTRACT_FASTA } from "../sc_gbc/modules/extract_fasta.nf"
include { CONSENSUS_BAM } from "../sc_gbc/modules/consensus_bam.nf"
include { ALLELIC_TABLES } from "./modules/allelic_tables.nf"
include { GATHER_TABLES } from "./modules/gather_allelic_tables.nf"
// include { TO_H5AD } from "./modules/to_h5ad.nf"

// 

process publish_maester {

    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), 
          path(bam),
          path(index),
          path(maegatk_out),
          path(afm)

    output:
    path bam
    path index
    path maegatk_out
    path afm

    script:
    """
    echo moving everything to ${params.sc_outdir}
    """

}
 
// 

//----------------------------------------------------------------------------//
// maester_pp subworkflow
//----------------------------------------------------------------------------//

workflow maester {
     
    take:
        ch_input
        filtered
        not_enriched_bam  

    main:

        // Get MT-reads from 10x and MAESTER libraries
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        ASSEMBLE_FQ(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
        STAR(ASSEMBLE_FQ.out.fq)
        FILTER_MAESTER_BAM(STAR.out.bam)
        FIX_TAGS(FILTER_MAESTER_BAM.out.bam)
        FILTER_10X_BAM(not_enriched_bam)
        MERGE_BAM(FILTER_10X_BAM.out.bam.combine(FIX_TAGS.out.bam, by:0))

        // Filter reads from good cells only, split into multiple bams and make consensus reads
        FILTER_MITOBAM(MERGE_BAM.out.bam.combine(filtered, by:0))
        SPLIT_BAM(FILTER_MITOBAM.out.filtered_mitobam)
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
        EXTRACT_FASTA(params.string_MT)
        CONSENSUS_BAM(ch_cell_bams, params.fgbio_min_reads_maester, EXTRACT_FASTA.out.fasta)

        // Create and aggregate cells allelic tables
        ALLELIC_TABLES(CONSENSUS_BAM.out.consensus_filtered_bam)
        ch_collapse = ALLELIC_TABLES.out.allelic_tables
            .map { it -> tuple(it[0], it[2], it[3], it[4], it[5], it[6]) }
            .groupTuple(by: 0)

        GATHER_TABLES(ch_collapse)

        //TO_H5AD(GATHER_TABLES.out.tables,EXTRACT_FASTA.out.ref_txt)
        // 
        // // Publish
        // publish_input = MERGE_BAM.out.bam
        //     .combine(MAEGATK.out.output, by:0)
        //     .combine(TO_H5AD.out.afm, by:0)
        // publish_maester(publish_input)

    emit:
        afm = GATHER_TABLES.out.tables

}

//----------------------------------------------------------------------------//