// nf-lenti
nextflow.enable.dsl = 2

// Include here
include { bulk_gbc } from "./subworkflows/bulk_gbc/main"
include { tenx } from "./subworkflows/tenx/main"
include { sc_gbc } from "./subworkflows/sc_gbc/main"

//


//----------------------------------------------------------------------------//
// nf-lenti pipeline
//----------------------------------------------------------------------------//

//

workflow BULK_GBC {
 
    ch = Channel.fromPath(params.raw_data_input)
        .splitCsv(header: true)
        .map { row -> [ row.ID_we_want, "${row.path_bulk}/${row.folder_name_bulk}"] }
    bulk_gbc(ch)

}

//

workflow TENX_GBC {

    // FIX AS IN nf-MiTo !!

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

