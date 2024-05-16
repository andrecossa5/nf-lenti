// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(cell_folder)

  output:
  tuple val(sample_name), val(cell), path("${cell}_consensus.bam"), emit: consensus_bam

  script:
  """
  # Group reads by UMI
  fgbio GroupReadsByUmi \
      -i ${cell_folder}/${cell}/${cell}.bam \
      -o grouped.bam \
      -s ${params.fgbio_UMI_consensus_mode} \
      -e ${params.fgbio_UMI_consensus_edits} \
      -t UB \
      -T MI    
  
  # Call Molecular consensus reads
  fgbio CallMolecularConsensusReads -t UB -i grouped.bam -o consensus.bam -M ${params.fgbio_min_reads} 

  # Filter high-quality consensus reads
  # fgbio FilterConsensusReads \
  #     -i consensus.bam \
  #     -o ${cell}_consensus_filtered.bam \
  #     -r ${params.ref}/new_genome_masked.fa \
  #     -M ${params.fgbio_min_reads} 
  #     -e 0.1 \                        # NB: need to be tuned from params.
  #     -N 20 \
  #     -E 0.025   

  # Index consensus_bam
  samtools index ${cell}_consensus.bam
  """

  stub:
  """
  touch ${cell}_consensus.bam
  """

}
