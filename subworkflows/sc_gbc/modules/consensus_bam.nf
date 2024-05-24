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
  fgbio="java -Xmx 4000m -jar + /maegatk/maegatk /bin/fgbio.jar" 
  fgbio GroupReadsByUmi \
      -i ${cell_folder}/${cell}/${cell}.bam \
      -o grouped.bam \

  fgbio="java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar"

  $fgbio GroupReadsByUmi \
      -i ${cell_folder}/${cell}.bam \ 
      -o ${cell_folder}/grouped.bam \
      -s ${params.fgbio_UMI_consensus_mode} \
      -e ${params.fgbio_UMI_consensus_edits} \
      -t UB \
      -T MI    

  # 2.5) Modify the UB tag
  #samtools view -H ${cell_folder}/grouped.bam > ${cell_folder}/grouped_temp.bam 
  #samtools view ${cell_folder}/grouped.bam | awk \'OFS="\t" {$13=$13""$4; print $0}\' >> ${cell_folder}/grouped_temp.sam
  #samtools view -b ${cell_folder}/grouped_temp.sam > ${cell_folder}/grouped_temp.bam 
  #os.system('echo "'+samtoolscall+'"')
  #os.system(samtoolscall)

  # 2.5) Modify the UB tag
  #samtools view -H TGTAACGCATTCTCTA/grouped.bam > TGTAACGCATTCTCTA/grouped_temp.bam 
  #samtools view TGTAACGCATTCTCTA/grouped.bam | awk \'OFS="\t" {$13=$13""$4; print $0}\' >> TGTAACGCATTCTCTA/grouped_temp.sam
  #samtools view -b TGTAACGCATTCTCTA/grouped_temp.sam > TGTAACGCATTCTCTA/grouped_temp.bam 
  
  # Call Molecular consensus reads
  # -i ${cell_folder}/grouped_temp.bam
  # -o {cell_folder}/consensus.bam
  $fgbio CallMolecularConsensusReads -t UB -i ${cell_folder}/grouped.bam -o ${cell_folder}/consensus.bam -M ${params.fgbio_min_reads} 
  

  # $fgbio CallMolecularConsensusReads -t UB -i TGTAACGCATTCTCTA/grouped.bam -o TGTAACGCATTCTCTA/consensus.bam -M 3


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