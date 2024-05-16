// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(cell_folder)

  output:
  tuple val(sample_name), val(cell), path('filtered_consensus.bam'), emit: filtered_consensus_bam

  //script aggiunto da Guido da riguardare se fila con input e output
  script:
  """
  touch filtered_consensus.bam
  # fgbio stuff + arguments --> 
  # cell_bams: un folder con tutti i bam delle singole cellule. Solo reads pulite e ben
  # allineate e clippate.


  fgbio --compression 1 FastqToBam -i ${cell_folder}/*R1.fastq ${cell_folder}/*R2.fastq \
        -r 16B12M +T -o unmapped.bam --sample MDA_clones --library UMI
  samtools fastq unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y ${params.ref} - | \
        fgbio -Xmx4g --compression 1 --async-io ZipperBams --unmapped unmapped.bam \
        --ref ${params.ref} --output mapped.bam
  fgbio GroupReadsByUmi -s Identity -e 0 -i mapped.bam -o grouped.bam --raw-tag RX -T MI    
  fgbio CallMolecularConsensusReads -t RX -i grouped.bam -o ${cell}_consensus.bam -M ${params.min_reads} 
  fgbio FilterConsensusReads -i ${cell}_consensus.bam -o ${cell}_consensus_filtered.bam -r ${params.ref} \
        -M ${params.min_reads} -e 0.1 -N 20 -E 0.025
  samtools index ${cell}_consensus_filtered.bam
  """

  stub:
  """
  touch filtered_consensus.bam
  """

}
