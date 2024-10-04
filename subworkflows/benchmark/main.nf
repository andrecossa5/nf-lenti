// measter subworkflow

// Include here
nextflow.enable.dsl = 2
include { SPLIT_BAM } from "../common/modules/split_bam.nf"
include { SAMTOOLS } from "./modules/samtools.nf"
include { COLLAPSE_SAMTOOLS } from "./modules/samtools.nf"
include { CELLSNP } from "./modules/cellsnp_lite.nf"
include { FREEBAYES } from "./modules/freebayes.nf"
include { COLLAPSE_FREEBAYES } from "./modules/freebayes.nf"
include { MAEGATK } from"./modules/maegatk.nf"
include { COLLAPSE_MAEGATK } from"./modules/maegatk.nf"

// 

//----------------------------------------------------------------------------//
// benchmark subworkflow
//----------------------------------------------------------------------------//

workflow benchmark {
     
    take:
        ch_bams  

    main:

        if (params.pp_method == "samtools") {

            matrices = SPLIT_BAM(ch_bams.map{it->tuple(it[0],it[1])})
            SAMTOOLS(SPLIT_BAM.out.cell_bams)
            matrices = COLLAPSE_SAMTOOLS(SAMTOOLS.out.calls)

        } else if (params.pp_method == "cellsnp-lite") {

            matrices = CELLSNP(ch_bams)

        } else if (params.pp_method == "freebayes") {

            SPLIT_BAM(ch_bams.map{it->tuple(it[0],it[1])})
            FREEBAYES(SPLIT_BAM.out.cell_bams)
            matrices = COLLAPSE_FREEBAYES(FREEBAYES.out.calls)

        } else if (params.pp_method == "maegatk") {

            SPLIT_BAM(ch_bams.map{it->tuple(it[0],it[1])})
            MAEGATK(SPLIT_BAM.out.cell_bams)
            matrices = COLLAPSE_MAEGATK(MAEGATK.out.tables)
        
        } else {
            
            println('Current benchmarking include: maegatk, samtools, cellsnp-lite, freebayes.')
        
        }

    emit:

        matrices = matrices

}

//----------------------------------------------------------------------------//