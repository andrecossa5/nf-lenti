// MERGE_R1 module

nextflow.enable.dsl = 2

//

process MERGE_R1 {

  tag "${sample_name}"

  input:
  tuple val(sample_name), val(in_folder)

  output:
  tuple val(sample_name), path("R1_raw.fastq.gz"), emit: R1

  script:
  """
  zcat ${in_folder}/*R1*.fastq.gz | pigz --fast -p ${task.cpus} > R1_raw.fastq.gz
  """

  stub:
  """
  touch R1_raw.fastq.gz
  """

}%
