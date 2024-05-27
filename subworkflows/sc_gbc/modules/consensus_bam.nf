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


  java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar GroupReadsByUmi -i ${cell_folder}/${cell}.bam -o ${cell_folder}/grouped.bam -s ${params.fgbio_UMI_consensus_mode} -e ${params.fgbio_UMI_consensus_edits} -t UB -T MI    
  
  java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar CallMolecularConsensusReads -t UB -i ${cell_folder}/grouped.bam -o ${cell_folder}/consensus.bam -M ${params.fgbio_min_reads} 
  
  # Filter high-quality consensus reads
  python ${baseDir}/bin/sc_gbc/filter_consensus.py \
  --input ${cell_folder}/consensus.bam \
  --output filtered.bam \
  --base_quality 30 \
  --mean_quality_th 30 \
  --E 0.025 \
  --e 0.1 \
  --mask_th 0.2
  """

  stub:
  """
  touch ${cell}_consensus.bam
  """

}

