// test_fgbio

// Include here
nextflow.enable.dsl = 2
include { FGBIO_CUSTOM_FILTER } from "./modules/fgbio_custom_filter.nf"
include { FGBIO_READ_FILTER } from "./modules/fgbio_read_filter.nf"
include { FGBIO_NO_FILTER } from "./modules/fgbio_no_filter.nf"
include { GATHER_TABLES } from "./modules/gather_tables.nf"
include { COLLAPSE_NANOPORE } from "./modules/collapse_nanopore.nf"
include { CONSENSUS_NANOPORE } from "./modules/consensus_nanopore.nf"
include { EXTRACT_FASTA } from "../common/modules/extract_fasta.nf"
include { TO_H5AD } from "../maester/modules/to_h5ad.nf"


// 


//----------------------------------------------------------------------------//
// test_fgbio subworkflow
//----------------------------------------------------------------------------//

//

process publish_MITO {

    publishDir "${params.test_outdir}/", mode: 'copy'

    input:
    tuple val(sample_name),  
        path(stats),
        path(afm)

    output:
    path(stats)
    path(afm)

    script:
    """
    echo moving everything to ${params.test_outdir}
    """

}

//


workflow test_fgbio {

    take:
        ch_bams

    main:

        // Prep channel
        ch_cell_bams = ch_bams
        .map { it ->
            def path_splitted = it.toString().split('/')
            def cell = path_splitted[-1].toString().split('\\.')[0]
                return [cell, it]
        }

        if (params.library == "MITO") {
        
            EXTRACT_FASTA(params.string_MT)

            // Make consensus reads, create and aggregate cells allelic tables. Compute testing stats
            if (params.test_fgbio_version == "fgbio_custom_filter") {
                out = FGBIO_CUSTOM_FILTER(ch_cell_bams, EXTRACT_FASTA.out.fasta)
            } else if (params.test_fgbio_version == "fgbio_read_filter") {
                out = FGBIO_READ_FILTER(ch_cell_bams, EXTRACT_FASTA.out.fasta)
            } else {
                out = FGBIO_NO_FILTER(ch_cell_bams, EXTRACT_FASTA.out.fasta)
            }
            GATHER_TABLES(out.allelic_tables.groupTuple(by: 0))
            TO_H5AD(GATHER_TABLES.out.tables, EXTRACT_FASTA.out.fasta.map{it -> it[0]})
            counts = TO_H5AD.out.afm
            publish_MITO(GATHER_TABLES.out.stats.combine(TO_H5AD.out.afm, by:0))

        } else if (params.library == "SCM-seq") {

            out = CONSENSUS_NANOPORE(ch_cell_bams)
            COLLAPSE_NANOPORE(out.allelic_tables.groupTuple(by: 0))
            counts = COLLAPSE_NANOPORE.out.tables

        }

    emit:
    
        results = counts

}


//----------------------------------------------------------------------------//