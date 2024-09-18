// CONSENSUS_MITO module

nextflow.enable.dsl = 2

//

process FGBIO_NO_FILTER {

  tag "${cell}"

  input:
  tuple val(cell), path(bam)
  tuple path(ref), 
      path(ref_dict),  
      path(ref_fa_amb),  
      path(ref_fa_ann),  
      path(ref_fa_bwt),  
      path(ref_fa_fai),  
      path(ref_fa_pac),  
      path(ref_fa_sa)

  output:
  tuple val("prova"),
    val(cell), 
    path("${cell}.A.txt"),
    path("${cell}.C.txt"),
    path("${cell}.T.txt"),
    path("${cell}.G.txt"), 
    path("${cell}.median_filtered_base_umi_group_size.txt"), 
    path("${cell}.n_umis_unfiltered.txt"), 
    path("${cell}.n_umis_filtered.txt"), 
    path("${cell}.depth.txt"), 
    path("${cell}.coverage.txt"), 
    path("${cell}.median_filtered_base_consensus_error.txt"),
    path("${cell}.median_filtered_read_quality.txt"),
    path("${cell}.n_reads_filtered.txt"),
    path("${cell}.n_reads_unfiltered.txt"), emit: allelic_tables

  script:
  """
  # fgbio consensus pipeline
  fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi \
	  --input ${bam} \
    --strategy ${params.fgbio_UMI_consensus_mode} \
    --edits ${params.fgbio_UMI_consensus_edits} \
	  --output grouped.bam \
	  -t UB \
	  -T MI

  fgbio -Xmx4g --compression 1 CallMolecularConsensusReads \
    --input grouped.bam \
    --output cons_unmapped.bam \
    --min-reads ${params.fgbio_min_reads_mito} \
    --min-input-base-quality ${params.fgbio_base_quality}

  samtools fastq cons_unmapped.bam \
  | bwa mem -t 16 -p -K 150000000 -Y ${ref} - \
  | fgbio -Xmx4g --compression 1 --async-io ZipperBams \
    --unmapped cons_unmapped.bam \
    --ref ${ref} \
    --tags-to-reverse Consensus \
    --tags-to-revcomp Consensus \
    --output cons_mapped.bam 

  # Create allelic tables
  python ${baseDir}/bin/test_fgbio/make_allelic_tables_test.py \
  --consensus_bam cons_mapped.bam \
  --cell ${cell} \
  --min_base_qual ${params.fgbio_base_quality} \
  --min_alignment_quality ${params.fgbio_min_alignment_quality} \
  --min_base_depth ${params.fgbio_min_reads_mito} \
  --min_base_consensus_error ${params.fgbio_base_error_rate_mito} \
  --filtering no_filter
  """

  stub:
  """
  touch ${cell}.A.txt
  touch ${cell}.C.txt
  touch ${cell}.T.txt
  touch ${cell}.G.txt
  touch ${cell}.median_filtered_base_umi_group_size.txt
  touch ${cell}.n_umis_unfiltered.txt
  touch ${cell}.n_umis_filtered.txt 
  touch ${cell}.depth.txt
  touch ${cell}.coverage.txt
  touch ${cell}.median_filtered_base_consensus_error.txt
  touch ${cell}.median_filtered_read_quality.txt
  touch ${cell}.n_reads_filtered.txt
  touch ${cell}.n_reads_unfiltered.txt
  """

}


