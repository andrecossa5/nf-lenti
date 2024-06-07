// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(cell_folder)

  output:
  tuple val(sample_name), val(cell), path("${cell_folder}_consensus_filtered_mapped.bam"), emit: consensus_filtered_bam

  script:
  """

  #fgbio GroupReadsByUmi -i ${cell_folder}/${cell}.bam -o${cell_folder}/grouped.bam -s ${params.fgbio_UMI_consensus_mode}  -e ${params.fgbio_UMI_consensus_edits} -t UB -T MI
#
  #fgbio CallMolecularConsensusReads -t UB -i ${cell_folder}/grouped.bam  -o ${cell_folder}/consensus_unmapped.bam   -M ${params.fgbio_min_reads}
#
  #samtools fastq ${cell_folder}/consensus_unmapped.bam \
  #  | bwa mem -t 16 -p -K 150000000 -Y ${params.ref}/cassette_up.fa - \
  #  | fgbio --compression 1 --async-io ZipperBams --unmapped ${cell_folder}/consensus_unmapped.bam --ref ${params.ref}/cassette_up.fa  --tags-to-reverse Consensus --tags-to-revcomp Consensus --output ${cell_folder}/consensus_mapped.bam 
#
  #fgbio -Xmx8g --compression 0 FilterConsensusReads \
  #  --input ${cell_folder}/consensus_mapped.bam \
  #  --output /dev/stdout \
  #  --ref ${params.ref}/cassette_up.fa  \
  #  --min-reads ${params.fgbio_min_reads} \
  #  --min-base-quality ${params.fgbio_base_quality} \
  #  --max-base-error-rate ${params.fgbio_base_error_rate} \
  #| samtools sort --threads 1 -o ${cell_folder}/consensus_filtered_mapped.bam --write-index
#


  fgbio GroupReadsByUmi -i ${cell}.bam -o ${cell}_grouped.bam -s ${params.fgbio_UMI_consensus_mode}  -e ${params.fgbio_UMI_consensus_edits} -t UB -T MI

  fgbio CallMolecularConsensusReads -t UB -i ${cell}_grouped.bam  -o ${cell_folder}_consensus_unmapped.bam   -M ${params.fgbio_min_reads}

  samtools fastq ${cell_folder}_consensus_unmapped.bam \
    | bwa mem -t 16 -p -K 150000000 -Y ${params.ref}/cassette_up.fa - \
    | fgbio --compression 1 --async-io ZipperBams --unmapped ${cell_folder}_consensus_unmapped.bam  --ref ${params.ref}/cassette_up.fa  --tags-to-reverse Consensus --tags-to-revcomp Consensus --output ${cell_folder}_consensus_mapped.bam 

  fgbio -Xmx8g --compression 0 FilterConsensusReads \
    --input ${cell_folder}_consensus_mapped.bam \
    --output /dev/stdout \
    --ref ${params.ref}/cassette_up.fa  \
    --min-reads ${params.fgbio_min_reads} \
    --min-base-quality ${params.fgbio_base_quality} \
    --max-base-error-rate ${params.fgbio_base_error_rate} \
  | samtools sort --threads 1 -o ${cell_folder}_consensus_filtered_mapped.bam --write-index

  """

  stub:
  """
  touch ${cell}_consensus.bam
  """

}

