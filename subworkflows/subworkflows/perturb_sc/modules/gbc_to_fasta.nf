// GBC_TO_FASTA module

nextflow.enable.dsl = 2

//

process GBC_TO_FASTA {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBCs)

  output:
  tuple val(sample_name), path("GBC_to_align.fa.gz"), emit: fasta

  script:
  """
  awk '{ gsub("@",">",\$1); print }' ${GBCs} \
  | tr ' ' '\n' \
  | pigz --fast -p ${task.cpus} \
  > GBC_to_align.fa.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_to_align.fa.gz
  """
}