// nf-lenti
nextflow.enable.dsl = 2

// Include here
include { tenx } from "./subworkflows/tenx/main"
include { bulk_gbc } from "./subworkflows/bulk_gbc/main"
include { sc_gbc } from "./subworkflows/sc_gbc/main"

//


//----------------------------------------------------------------------------//
// mito_preprocessing pipeline
//----------------------------------------------------------------------------//

//

workflow TENX {

    // 10x GEX library
    ch_tenx = Channel
        .fromPath("${params.sc_tenx_indir}/*", type:'dir')
        .map{ tuple(it.getName(), it) }

    tenx(ch_tenx)

}

//

workflow BULK_GBC {
 
    // (Bulk DNA) targeted DNA sequencing of GBC
    ch_bulk_gbc = Channel
        .fromPath(params.bulk_gbc_indir) 
        .splitCsv(header : true)
        .map{ row -> 
            [row.ID_we_want, "${row.path_bulk}/${row.folder_name_bulk}"]}
    bulk_gbc(ch_bulk_gbc)

}

//

workflow TENX_GBC {

    // GBC enrichment from 10x library
    ch_sc_gbc = Channel
        .fromPath("${params.sc_gbc_indir}/*", type:'dir')
        .map{ tuple(it.getName(), it) }
    
    // 10x GEX library
    ch_tenx = Channel
        .fromPath("${params.sc_tenx_indir}/*", type:'dir')
        .map{ tuple(it.getName(), it) }

    tenx(ch_tenx)
    sc_gbc(ch_sc_gbc, tenx.out.filtered)
    sc_gbc.out.summary.view()

}

//

workflow {
    
    Channel.of(1,2,3,4) | view

}

//

