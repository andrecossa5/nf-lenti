// MERGE_GBC module

nextflow.enable.dsl = 2

//

process MERGE_GBC {

  tag "${sample_name}"

  input:
  tuple val(sample_name), val(in_folder)

  output:
  tuple val(sample_name), path("R1_gbc.fastq.gz"), path("R2_gbc.fastq.gz"), emit: reads

  script:
  """
  zcat ${in_folder}/*R1*.fastq.gz \
  | awk '{if(NR%4==1){print "@"(NR%1?c+1:++c)} else {print \$0}}' \
  | pigz --fast -p ${task.cpus}  \
  > R1_gbc.fastq.gz

  zcat ${in_folder}/*R2*.fastq.gz \
  | awk '{if(NR%4==1){print "@"(NR%1?c+1:++c)} else {print \$0}}' \
  | pigz --fast -p ${task.cpus}  \
  > R2_gbc.fastq.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch R1_gbc.fastq.gz
  touch R2_gbc.fastq.gz
  """

}