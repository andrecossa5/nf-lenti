// MERGE module

nextflow.enable.dsl = 2

//

process MERGE_FQ {

  tag "${sample_name}"

  input:
  tuple val(sample_name), val(in_folder)

  output:
  tuple val(sample_name), path("R1.fastq.gz"), path("R2.fastq.gz"), emit: reads

  script:
  """
  zcat ${in_folder}/*R1*.fastq.gz \
  | awk '{if(NR%4==1){print "@"(NR%1?c+1:++c)} else {print \$0}}' \
  | pigz --fast -p ${task.cpus}  \
  > R1.fastq.gz

  zcat ${in_folder}/*R2*.fastq.gz \
  | awk '{if(NR%4==1){print "@"(NR%1?c+1:++c)} else {print \$0}}' \
  | pigz --fast -p ${task.cpus}  \
  > R2.fastq.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch R1.fastq.gz
  touch R2.fastq.gz
  """

}