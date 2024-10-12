// measter subworkflow

// Include here
nextflow.enable.dsl = 2
include { SPLIT_BARCODES } from "./modules/split_barcodes.nf"
include { FILTER_BAM_CB } from "../common/modules/filter_bam_cb.nf"
include { SPLIT_BAM } from "../common/modules/split_bam.nf"
include { EXTRACT_FASTA } from "../common/modules/extract_fasta.nf"
include { SAMTOOLS } from "./modules/samtools.nf"
include { COLLAPSE_SAMTOOLS } from "./modules/samtools.nf"
include { INDEX_AND_MERGE } from "./modules/index_and_merge.nf"
include { CELLSNP } from "./modules/cellsnp_lite.nf"
include { FREEBAYES } from "./modules/freebayes.nf"
include { COLLAPSE_FREEBAYES } from "./modules/freebayes.nf"
include { MAEGATK } from "./modules/maegatk.nf"
include { COLLAPSE_MAEGATK } from "./modules/maegatk.nf"

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
// benchmark subworkflow
//----------------------------------------------------------------------------//

workflow benchmark {
     
    take:
        ch_bams  

    main:

        EXTRACT_FASTA(params.string_MT)
        SPLIT_BARCODES(ch_bams)
        ch_barcodes = SPLIT_BARCODES.out.barcodes.flatMap { 
            sample_name, bam, file_paths ->
            if (file_paths instanceof List) {
                file_paths.collect { file_path -> [sample_name, bam, file_path] }
            } else {
                [[sample_name, bam, file_paths]]
            }
        }
        FILTER_BAM_CB(ch_barcodes)

        if (params.pp_method == "samtools") {

            SPLIT_BAM(FILTER_BAM_CB.out.bam)
            ch_cell_bams = processCellBams(SPLIT_BAM.out.cell_bams)
            SAMTOOLS(ch_cell_bams, EXTRACT_FASTA.out.fasta.map{it->it[0]})
            COLLAPSE_SAMTOOLS(SAMTOOLS.out.calls.groupTuple(by:0))
            ch_output = COLLAPSE_SAMTOOLS.out.ch_output

        } else if (params.pp_method == "cellsnp-lite") {

            INDEX_AND_MERGE(FILTER_BAM_CB.out.bam.groupTuple(by:0))
            ch_cellsnp = INDEX_AND_MERGE.out.bam
                        .combine(ch_bams.map{it->tuple(it[0],it[2])}, by:0)
            CELLSNP(ch_cellsnp)
            ch_output = CELLSNP.out.ch_output

        } else if (params.pp_method == "freebayes") {

            SPLIT_BAM(FILTER_BAM_CB.out.bam)
            ch_cell_bams = processCellBams(SPLIT_BAM.out.cell_bams) 
            FREEBAYES(ch_cell_bams, EXTRACT_FASTA.out.fasta.map{it->tuple(it[0],it[5])})
            COLLAPSE_FREEBAYES(FREEBAYES.out.calls.groupTuple(by:0))
            ch_output = COLLAPSE_FREEBAYES.out.ch_output

        } else if (params.pp_method == "maegatk") {

            SPLIT_BAM(FILTER_BAM_CB.out.bam)
            ch_cell_bams = processCellBams(SPLIT_BAM.out.cell_bams)
            MAEGATK(ch_cell_bams, EXTRACT_FASTA.out.fasta)
            COLLAPSE_MAEGATK(MAEGATK.out.tables.groupTuple(by:0)) 
            ch_output = COLLAPSE_MAEGATK.out.ch_output
        
        } else {
            
            println('Current benchmarking include: maegatk, samtools, cellsnp-lite, freebayes.')
        
        }

    emit:

        ch_output = ch_output

}

//----------------------------------------------------------------------------//