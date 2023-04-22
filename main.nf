// MI_TO pipeline
nextflow.enable.dsl = 2
include { perturb_sc } from "./subworkflows/perturb_sc/main"
include { maester } from "./subworkflows/maester/main"
include { tenx } from "./subworkflows/tenx/main"

// Perturb-seq input sc_fastqs 
ch_perturb_sc = Channel
    .fromPath("${params.perturb_sc_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

//

// MAESTER input sc_fastq
ch_maester = Channel
    .fromPath("${params.maester_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

// 

// 10x only input sc_fastqs
ch_tenx = Channel
    .fromPath("${params.tenx_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

//

//----------------------------------------------------------------------------//
// Perturb-seq pipeline entry points and (mock) main workflow
//----------------------------------------------------------------------------//

//

workflow tenx_only {

    tenx(ch_tenx)
    tenx.out.filtered.view()
    tenx.out.bam.view()

}

//

workflow perturbseq_only {

    perturb_sc(ch_perturb_sc)
    perturb_sc.out.filtered.view()
    perturb_sc.out.bam.view()

}

//

workflow tenx_mito {

    tenx(ch_tenx)
    maester(ch_maester, tenx.out.filtered, tenx.out.bam)
    maester.out.afm.view()

}

//

workflow gbc_mito {

    perturb_sc(ch_perturb_sc)
    maester(ch_maester, perturb_sc.out.filtered, perturb_sc.out.bam)
    maester.out.afm.view()

}

//

// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}