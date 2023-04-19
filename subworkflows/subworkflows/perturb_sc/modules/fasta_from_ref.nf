// FASTA_FROM_REF module

nextflow.enable.dsl = 2

//

process FASTA_FROM_REF {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), val(in_folder)

  output:
  tuple val(sample_name), path("GBC_reference.fa"), emit: fasta

  script:
  """
  awk 'FNR > 1 {print ">"NR-1"\\n"\$1}' \
  ${params.perturb_bulk_out}/${sample_name}/read_count_by_GBC_corrected.tsv > GBC_reference.fa
  """

}