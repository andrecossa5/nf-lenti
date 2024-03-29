// MAEGATK module

nextflow.enable.dsl = 2

//

process MAEGATK {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(mitobam), path(index), path(filtered)

  output:
  tuple val(sample_name), path("final"), emit: output
 
  script:
  """
  zcat ./filtered/barcodes.tsv.gz > ./barcodes.txt

  python ${baseDir}/bin/maester/maegatk_cli.py \
  ${params.maester_code_dir} \
  ${mitobam} \
  ${task.cpus} \
  ./barcodes.txt \
  ${params.maester_min_reads} \
  ${params.min_base_qual} \
  ${params.min_alignment_quality}
  """

  stub:
  """
  mkdir final
  """

}