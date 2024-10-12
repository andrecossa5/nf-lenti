// maester subworkflow

// Include here
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"
include { FILTER_10X_BAM } from "./modules/filter_bam.nf"
include { FILTER_MAESTER_BAM } from "./modules/filter_bam.nf"
include { MERGE_BAM } from "./modules/merge_bams.nf"
include { EXTRACT_FASTA } from "../common/modules/extract_fasta.nf"
include { SPLIT_BARCODES } from "../common/modules/split_barcodes.nf"
include { FILTER_BAM_CB } from "../common/modules/filter_bam_cb.nf"
include { SPLIT_BAM } from "../common/modules/split_bam.nf"
include { CONSENSUS_MITO } from "./modules/consensus_mito.nf"
include { GATHER_TABLES } from "./modules/gather_allelic_tables.nf"

// 

process publish_maester {

    publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name),  
        path(bam),
        path(tables),
        path(fasta)

    output:
    path(bam)
    path(tables)
    path(afm)

    script:
    """
    echo moving everything to ${params.sc_outdir}
    """

}
 
// 

def processCellBams(cell_bams) {
    return cell_bams
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
}

//

//----------------------------------------------------------------------------//
// maester subworkflow
//----------------------------------------------------------------------------//

workflow maester {
     
    take:
        ch_input
        ch_filtered
        not_enriched_bam  

    main:

        // Get MT-reads from 10x and MAESTER libraries
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
        FILTER_MAESTER_BAM(SOLO.out.bam)
        FILTER_10X_BAM(not_enriched_bam)
        MERGE_BAM(FILTER_10X_BAM.out.bam.combine(FILTER_MAESTER_BAM.out.bam, by:0))

        // Filter reads from good cells only, split into multiple bams 
        SPLIT_BARCODES(ch_filtered)
        ch_barcodes = SPLIT_BARCODES.out.barcodes.flatMap { 
            sample_name, file_paths ->
            if (file_paths instanceof List) {
                file_paths.collect { file_path -> [sample_name, file_path] }
            } else {
                [[sample_name, file_paths]]
            }
        }
        FILTER_BAM_CB(MERGE_BAM.out.bam.combine(ch_barcodes, by:0))
        SPLIT_BAM(FILTER_BAM_CB.out.bam)
        ch_cell_bams = processCellBams(SPLIT_BAM.out.cell_bams)
        EXTRACT_FASTA(params.string_MT)

        // Make consensus reads, create and aggregate cells allelic tables
        CONSENSUS_MITO(ch_cell_bams, EXTRACT_FASTA.out.fasta)
        GATHER_TABLES(CONSENSUS_MITO.out.allelic_tables.groupTuple(by: 0))
        
        // Publish
        ch_pubb = MERGE_BAM.out.bam
                  .combine(GATHER_TABLES.out.tables, by:0)
                  .combine(EXTRACT_FASTA.out.fasta.map{it->it[0]})
        publish_maester(ch_pubb)

    emit:

        allelic_tables = CONSENSUS_MITO.out.allelic_tables

}

//----------------------------------------------------------------------------//