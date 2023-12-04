// mi_to_preprocessing
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

// 10x expression data
ch_tenx = Channel
    .fromPath("${params.sc_tenx_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

// 10x sub-library, from target GBC enrichment
ch_sc_gbc = Channel
    .fromPath("${params.sc_gbc_indir}/*", type:'dir')
    .map{ tuple(it.getName(), it) }

// MAESTER
ch_maester = Channel
    .fromPath("${params.sc_maester_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }


//

//----------------------------------------------------------------------------//
// mito_preprocessing pipeline: new version
//----------------------------------------------------------------------------//

//

workflow TENX {

    tenx(ch_tenx)
    tenx.out.filtered.view()

}

//

workflow TENX_MITO {

    tenx(ch_tenx)
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

    sc_gbc(ch_tenx, ch_sc_gbc)
    sc_gbc.out.summary.view()

}

//

workflow SC_GBC_MITO {

    sc_gbc(ch_tenx, ch_sc_gbc)
    maester(ch_maester, sc_gbc.out.filtered, sc_gbc.out.bam)
    maester.out.afm.view()

}

//

// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}