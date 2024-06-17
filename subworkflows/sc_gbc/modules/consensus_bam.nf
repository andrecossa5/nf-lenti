// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(bam)
  val(min_reads)
  path(ref)

  output:
  tuple val(sample_name), val(cell), path("${cell}_consensus_filtered_mapped.bam"), emit: consensus_filtered_bam

  script:
  """
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
    --min-reads ${min_reads} \
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
    --min-reads ${min_reads} \
    --min-base-quality ${params.fgbio_base_quality} \
    --max-base-error-rate ${params.fgbio_base_error_rate} \
    | samtools sort -@ 1 -o "${cell}_consensus_filtered_mapped.bam" --write-index

  """

  stub:
  """
  touch ${cell}_consensus_filtered_mapped.bam
  """

}

