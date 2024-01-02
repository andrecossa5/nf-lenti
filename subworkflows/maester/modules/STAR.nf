// STAR module

nextflow.enable.dsl = 2

//

process STAR {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(fastq)

  output:
  tuple val(sample_name), path("Aligned.out.bam"), emit: bam

  script:
  """
  STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${params.ref} \
    --readFilesIn ${fastq} \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI nM AS
  """

  stub:
  """
  touch Aligned.out.bam
  """

}