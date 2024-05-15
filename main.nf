// Main pipeline mi_to_preprocessing
nextflow.enable.dsl = 2
include { bulk_gbc } from "./subworkflows/bulk_gbc/main"
include { tenx } from "./subworkflows/tenx/main"
include { sc_gbc } from "./subworkflows/sc_gbc/main"
include { maester } from "./subworkflows/maester/main"

//

// Bulk DNA target GBC sequencing
ch_bulk_gbc = Channel
    .fromPath("${params.bulk_gbc_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }
    
// 10x and GBC data, single-cell
ch_sc = Channel
    .fromPath("${params.sc_indir}/*", type:'dir')
    .map{ tuple(it.getName(), it) }

// MAESTER
ch_maester = Channel
    .fromPath("${params.sc_maester_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }


//


//----------------------------------------------------------------------------//
// mito_preprocessing pipeline
//----------------------------------------------------------------------------//

//

workflow TENX {

    tenx(ch_sc)
    tenx.out.filtered.view()

}

//

workflow TENX_MITO {

    tenx(ch_sc)
    maester(ch_maester, tenx.out.filtered, tenx.out.bam)
    maester.out.afm.view()

}

//

workflow BULK_GBC {
 
    bulk_gbc(ch_bulk_gbc)
    bulk_gbc.out.flags.view()

}

//

workflow SC_GBC {

    sc_gbc(ch_sc)
    sc_gbc.out.summary.view()

}

//

workflow SC_GBC_MITO {

    sc_gbc(ch_sc)
    maester(ch_maester, sc_gbc.out.filtered, sc_gbc.out.bam)
    maester.out.afm.view()

}

//

// Mock
workflow {
    
    Channel.of(1,2,3,4) | view

}
