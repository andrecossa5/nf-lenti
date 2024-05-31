// CONSENSUS_BAM module

nextflow.enable.dsl = 2

//

process CONSENSUS_BAM {

  tag "${sample_name}: ${cell}

  input:
  tuple val(sample_name), val(cell), path(cell_folder)

  output:
  tuple val(sample_name), val(cell), path("${cell_folder}/consensus_filtered.bam"), emit: consensus_filtered_bam
  tuple val(sample_name), path("seq_qual_err.csv"), emit: filter_summary

  script:
  """


  java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar GroupReadsByUmi -i ${cell_folder}/${cell}.bam -o ${cell_folder}/grouped.bam -s ${params.fgbio_UMI_consensus_mode} -e ${params.fgbio_UMI_consensus_edits} -t UB -T MI    
  
  java -Xmx4000m -jar /maegatk/maegatk/bin/fgbio.jar CallMolecularConsensusReads -t UB -i ${cell_folder}/grouped.bam -o ${cell_folder}/consensus.bam -M ${params.fgbio_min_reads} 
  
  # Filter high-quality consensus reads
  python ${baseDir}/bin/sc_gbc/filter_consensus.py \
  --input ${cell_folder}/consensus.bam \
  --output ${cell_folder}/consensus_filtered.bam \
  --min_quality 30 \
  --read_max_consensus_error 0.025 \
  --base_consensus_error 0.1 \
  --read_max_N 0.2
  --consensus_filter_mode 'GBC'
  --GBC_max_N 


  # Params
min_quality = args.min_quality
read_max_consensus_error = args.read_max_consensus_error
base_consensus_error = args.base_consensus_error
read_max_N = args.read_max_N
consensus_filter_mode = args.consensus_filter_mode
GBC_max_N =args.GBC_max_N
  """

  stub:
  """
  touch ${cell}_consensus.bam
  """

}

