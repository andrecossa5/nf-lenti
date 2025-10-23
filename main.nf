// nf-lenti, tested
nextflow.enable.dsl = 2

// Include here
include { bulk_gbc } from "./subworkflows/bulk_gbc/main"
include { tenx_gbc } from "./subworkflows/tenx_gbc/main"
include { createPreprocessingChannel } from "./subworkflows/tenx_gbc/main"

//


//----------------------------------------------------------------------------//
// nf-lenti pipeline entrypoints
//----------------------------------------------------------------------------//

//

workflow BULK_GBC {
 
    ch_preprocessing = Channel.fromPath(params.raw_data_input)
        .splitCsv(header: true)
        .map { row -> [ row.ID_we_want, "${row.path_bulk}/${row.folder_name_bulk}"] }
    bulk_gbc(ch_preprocessing)

}

//

workflow TENX_GBC {

    ch_preprocessing = createPreprocessingChannel()
    tenx_gbc(ch_preprocessing)

}

//

workflow {
    
    Channel.of(1,2,3,4) | view

}

//

