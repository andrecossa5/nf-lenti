// ALIGN_33_R2 module

nextflow.enable.dsl = 2

//

process ASSEMBLE_FQ {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(R1), path(R2)

  output:
  tuple val(sample_name), path("assembled_fastq.gz"), emit: fq

  script:
  """
  python ${baseDir}/bin/maester/assemble_trim_fastq.py ${R1} ${R2} ${task.cpus}
  """

  stub:
  """
  touch assembled_fastq.gz
  """

}