// CONSENSUS_LENTI module

nextflow.enable.dsl = 2

//

process CONSENSUS_LENTI {

  label 'scLT'
  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), 
      val(cell), 
      path(bam)
  tuple path(ref), 
      path(ref_dict),  
      path(ref_fa_amb),  
      path(ref_fa_ann),  
      path(ref_fa_bwt),  
      path(ref_fa_fai),  
      path(ref_fa_pac),  
      path(ref_fa_sa)

  output:
  tuple val(sample_name), val(cell), path("${cell}_consensus_filtered.tsv"), emit: consensus_filtered_tsv

  script:
  """
  # Samtools sort and extract reads from lentiviral cassette
  samtools sort -@ 1 ${bam} --write-index -o sorted.bam 
  samtools view -b sorted.bam ${params.string_lentiviral} > filtered.bam

  ##

  # fgbio consensus pipeline
  fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi \
	  --input filtered.bam \
    --strategy ${params.fgbio_UMI_consensus_mode} \
    --edits ${params.fgbio_UMI_consensus_edits} \
	  --output grouped.bam \
	  -t UB \
	  -T MI

  fgbio -Xmx4g --compression 1 CallMolecularConsensusReads \
    --input grouped.bam \
    --output cons_unmapped.bam \
    --min-reads ${params.fgbio_min_reads_gbc} \
    --min-input-base-quality ${params.fgbio_base_quality}

  samtools fastq cons_unmapped.bam \
  | bwa mem -t 16 -p -K 150000000 -Y ${ref} - \
  | fgbio -Xmx4g --compression 1 --async-io ZipperBams \
    --unmapped cons_unmapped.bam \
    --ref ${ref} \
    --tags-to-reverse Consensus \
    --tags-to-revcomp Consensus \
    --output cons_mapped.bam 

  fgbio -Xmx8g --compression 0 FilterConsensusReads \
    --input cons_mapped.bam \
    --output /dev/stdout \
    --ref ${ref} \
    --min-reads ${params.fgbio_min_reads_gbc} \
    --min-base-quality ${params.fgbio_base_quality} \
    --max-base-error-rate ${params.fgbio_base_error_rate_gbc} \
    | samtools sort -@ 1 -o consensus_filtered_mapped.bam --write-index
  
  
  ##


  # Create tables
  python ${baseDir}/bin/sc_gbc/consensus_tsv.py consensus_filtered_mapped.bam ${cell}
  """

  stub:
  """
  touch ${cell}_consensus_filtered.tsv
  """

}

