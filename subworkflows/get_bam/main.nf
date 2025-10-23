// get_gbc_bam subworkflows

// Modules
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"

 
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

//----------------------------------------------------------------------------//
// get_gbc_bam subworkflow
//----------------------------------------------------------------------------//

workflow get_tenx_gbc_bam {
     
    take:
        ch_fastqs
        cell_barcodes

    main:

        MERGE_R1(ch_fastqs)
        MERGE_R2(ch_fastqs)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))

    emit:
        tenx_bam = SOLO.out.bam.combine(cell_barcodes, by: 0)


}

workflow get_gbc_bam {

    take:
        ch_fastqs
        cell_barcodes

    main:

        MERGE_R1(ch_fastqs)
        MERGE_R2(ch_fastqs)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))

    emit:
        gbc_bam = SOLO.out.bam.combine(cell_barcodes, by: 0)

}