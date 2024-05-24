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
  fgbio="java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar"

  $fgbio GroupReadsByUmi \
      -i ${cell_folder}/${cell}.bam \ 
      -o ${cell_folder}/grouped.bam \
      -s ${params.fgbio_UMI_consensus_mode} \
      -e ${params.fgbio_UMI_consensus_edits} \
      -t UB \
      -T MI    
  
  $fgbio CallMolecularConsensusReads -t UB -i ${cell_folder}/grouped.bam -o ${cell_folder}/consensus.bam -M ${params.fgbio_min_reads} 
  
  # Filter high-quality consensus reads
  python ${baseDir}/bin/sc_gbc/filter_consensus.py \
  --input consensus.bam \
  --output filtered.bam \
  --base_quality 30 \
  --mean_quality_th 30 \
  --E 0.025 \
  --e 0.1 \
  --mask_th 0.2
  #$fgbio FilterConsensusReads -i TGTAACGCATTCTCTA/consensus.bam  -o TGTAACGCATTCTCTA/TGTAACGCATTCTCTA_consensus_filtered.bam -r $ref -M 3 -e 0.1 -N 20 -E 0.025



  #attenzione a params ref
  # fgbio FilterConsensusReads \
  #     -i consensus.bam \
  #     -o ${cell}_consensus_filtered.bam \
  #     -r ${params.ref}/new_genome_masked.fa \
  #     -M ${params.fgbio_min_reads} 
  #     -e 0.1 \
  #     -N 20 \
  #     -E 0.025   

  # Index consensus_bam
  #samtools index ${cell}_consensus.bam
  """

  stub:
  """
  touch ${cell}_consensus.bam
  """

}


//fgbio FilterConsensusReads -i consensus.bam -o consensus_filtered.bam -r $ref -M 3 -e 0.1 -N 20 -E 0.025 