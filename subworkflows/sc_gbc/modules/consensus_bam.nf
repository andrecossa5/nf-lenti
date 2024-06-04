// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(cell_folder)

  output:
  tuple val(sample_name), val(cell), path("${cell_folder}/consensus_filtered.bam"), emit: consensus_filtered_bam
  tuple val(sample_name), path("seq_qual_err.csv"), emit: filter_summary

  script:
  """


  java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar GroupReadsByUmi -i ${cell_folder}/${cell}.bam -o ${cell_folder}/grouped.bam -s ${params.fgbio_UMI_consensus_mode} -e ${params.fgbio_UMI_consensus_edits} -t UB -T MI    
  
  java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar CallMolecularConsensusReads -t UB -i ${cell_folder}/grouped.bam -o ${cell_folder}/consensus_unmapped.bam -M ${params.fgbio_min_reads} 
  

  samtools fastq ${cell_folder}/consensus_unmapped.bam \
    | bwa mem -t 16 -p -K 150000000 -Y ${ref} - \
    | java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar -Xmx4g --compression 1 --async-io ZipperBams \
        --unmapped ${cell_folder}/consensus_unmapped.bam \
        --ref ${params.ref} \
        --tags-to-reverse Consensus \
        --tags-to-revcomp Consensus \
        --output ${cell_folder}/consensus_mapped.bam 

  java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar -Xmx8g --compression 0 FilterConsensusReads \
    --input ${cell_folder}/consensus_mapped.bam \
    --output /dev/stdout \
    --ref ${params.ref} \
    --min-reads 3 \
    --min-base-quality 45 \
    --max-base-error-rate 0.2 \
    | samtools sort --threads 8 -o ${cell_folder}/consensus_filtered_mapped.bam --write-index

  """

  stub:
  """
  touch ${cell}_consensus.bam
  """

}

