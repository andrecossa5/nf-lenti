// GET_LENTIBAM module

nextflow.enable.dsl = 2

//

process GET_LENTIBAM {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("grepped.txt.gz"), emit: lentiviral_records

  script:
  """
  samtools view -h ${bam} \
  | grep lentiCassette \
  | pigz --fast -p ${task.cpus} \
  > grepped.txt.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch grepped.txt.gz
  """

}