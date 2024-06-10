// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(bam)

  output:
  tuple val(sample_name), val(cell), path("${cell}_consensus_filtered_mapped.bam"), emit: consensus_filtered_bam

  script:
  """
  
  fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi \
    --input ${bam} \
    --strategy ${params.fgbio_UMI_consensus_mode} \
    --edits ${params.fgbio_UMI_consensus_edits} \
    --output /dev/stdout \
    -t UB \
    -T MI \
  | fgbio -Xmx4g --compression 1 CallMolecularConsensusReads \
    --input /dev/stdin \
    --output /dev/stdout \
    --min-reads ${params.fgbio_min_reads} \
    --min-input-base-quality ${params.fgbio_base_quality} \
  | samtools fastq - \
  | bwa mem -t 16 -p -K 150000000 -Y ${params.ref}/cassette_up.fa - \
  | fgbio -Xmx4g --compression 1 --async-io ZipperBams \
      --unmapped /dev/stdin \
      --ref ${params.ref}/cassette_up.fa \
      --tags-to-reverse Consensus \
      --tags-to-revcomp Consensus \
      --output /dev/stdout \
  | samtools sort -n -u - \
  | fgbio -Xmx8g --compression 0 FilterConsensusReads \
    --input /dev/stdin \
    --output /dev/stdout \
    --ref ${params.ref}/cassette_up.fa \
    --min-reads ${params.fgbio_min_reads} \
    --min-base-quality ${params.fgbio_base_quality} \
    --max-base-error-rate ${params.fgbio_base_error_rate} \
  | samtools sort -@ 1 -o ${cell}_consensus_filtered_mapped.bam --write-index

  """

  stub:
  """
  touch ${cell}_consensus_filtered_mapped.bam
  """

}

