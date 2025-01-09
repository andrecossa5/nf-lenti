// mi_to_preprocessing
nextflow.enable.dsl = 2

// Include here
include { bulk_gbc } from "./subworkflows/bulk_gbc/main"
include { tenx } from "./subworkflows/tenx/main"
include { sc_gbc } from "./subworkflows/sc_gbc/main"
include { maester } from "./subworkflows/maester/main"
include { benchmark } from "./subworkflows/benchmark/main"    


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

workflow TENX_MITO {

    // 10x GEX library
    ch_tenx = Channel
        .fromPath("${params.sc_tenx_indir}/*", type:'dir')
        .map{ tuple(it.getName(), it) }

    // MAESTER library
    ch_maester = Channel
        .fromPath("${params.sc_maester_indir}/*", type:'dir') 
        .map{ tuple(it.getName(), it) }

    tenx(ch_tenx)
    maester(ch_maester, tenx.out.filtered, tenx.out.bam)

} 

//

workflow BULK_GBC {
 
    // (Bulk DNA) targeted DNA sequencing of GBC
    ch_bulk_gbc = Channel
        .fromPath(params.bulk_gbc_indir) 
        .splitCsv(header : true)
        .map{ row -> 
            //[row.ID_we_want, "${row.path_bulk}/${row.folder_name_bulk}"]}
            [row.id_want, "${row.path}/${row.path_sample}"]}
    //ch_bulk_gbc.view()
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

workflow TENX_GBC_MITO {

    // GBC enrichment from 10x library
    ch_sc_gbc = Channel
        .fromPath("${params.sc_gbc_indir}/*", type:'dir')
        .map{ tuple(it.getName(), it) }

    // 10x GEX library
    ch_tenx = Channel
        .fromPath("${params.sc_tenx_indir}/*", type:'dir')
        .map{ tuple(it.getName(), it) }

    // MAESTER library
    ch_maester = Channel
        .fromPath("${params.sc_maester_indir}/*", type:'dir') 
        .map{ tuple(it.getName(), it) }

    tenx(ch_tenx)
    sc_gbc(ch_sc_gbc, tenx.out.filtered)
    maester(ch_maester, tenx.out.filtered, tenx.out.bam)
    maester.out.allelic_tables.view()

}

//

workflow BENCH {

    // Bench
    ch_bams = ch_jobs = Channel.fromPath(params.bam_file)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.bam, row.barcodes ]}

    benchmark(ch_bams)
    benchmark.out.ch_output.view()

}

//

workflow {
    
    Channel.of(1,2,3,4) | view

}

//

