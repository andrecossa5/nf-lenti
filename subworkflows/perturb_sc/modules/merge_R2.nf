// MERGE_R2 module

nextflow.enable.dsl = 2

//

process MERGE_R2 {

  tag "${sample_name}"

  input:
  tuple val(sample_name), val(in_folder)

  output:
  tuple val(sample_name), path("R2_raw.fastq.gz"), emit: R2

  script:
  """
  zcat ${in_folder}/*R2*.fastq.gz \
  | awk '{if(NR%4==1){print "@"(NR%1?c+1:++c)} else {print \$0}}' \
  | pigz --fast -p ${task.cpus}  \
  > R2_raw.fastq.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch R2_raw.fastq.gz
  """

}