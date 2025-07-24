// nf-lenti
nextflow.enable.dsl = 2

// Include here
include { bulk_gbc } from "./subworkflows/bulk_gbc/main"
include { tenx } from "./subworkflows/tenx/main"
include { sc_gbc } from "./subworkflows/sc_gbc/main"
include { get_tenx_gbc_bam } from "./subworkflows/get_bam/main"
include { get_gbc_bam } from "./subworkflows/get_bam/main"

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

def createPreprocessingChannel() {

    if (params.raw_data_input_type == "fastq") {

        // From raw reads, unaligned
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.fastq_folder, row.library ] }
        
    } else if (params.raw_data_input_type == "fastq,GBC") {

        // From raw reads, unaligned (GBC) and a .txt file of valid 10x barcodes
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.fastq_folder, row.cell_barcodes ] }
        
    } else if (params.raw_data_input_type == "bam") {

        // From aligned reads, and a .txt file of valid 10x barcodes
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.bam, row.cell_barcodes ] }
    }
    else {
        error 'Unsupported raw_data_input_type: ${params.raw_data_input_type}. Available: "fastq", "fastq,GBC", or "bam".'
    }
    return ch
}


workflow tenx_gbc {

    take:
        ch

    main:
    
        if (params.raw_data_input_type == "fastq") {

            // All 10x and GBC reads
            tenx_fastqs = ch.filter{it->it[2]=='TENX'}.map{it->tuple(it[0],it[1])}
            gbc_fastqs = ch.filter{it->it[2]=='GBC'}.map{it->tuple(it[0],it[1])}
            tenx(tenx_fastqs)
            get_tenx_gbc_bam(gbc_fastqs, tenx.out.cell_barcodes)
            ch_filtered = get_tenx_gbc_bam.out.tenx_bam

        } else if (params.raw_data_input_type == "fastq,GBC") {

            // GBC reads
            gbc_fastqs = ch.map{it->tuple(it[0],it[1])}
            cell_barcodes = ch.map{it->tuple(it[0],it[2])}
            get_gbc_bam(gbc_fastqs, cell_barcodes)
            ch_filtered = get_gbc_bam.out.gbc_bam

        } else if (params.raw_data_input_type == "bam") {

            // Previously aligned reads
            ch_filtered = ch

        }
        sc_gbc(ch_filtered)
}

workflow TENX_GBC {

    ch_preprocessing = createPreprocessingChannel()
    tenx_gbc(ch_preprocessing)

}


//

workflow {
    
    Channel.of(1,2,3,4) | view

}

//

