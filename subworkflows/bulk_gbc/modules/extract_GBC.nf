// EXTRACT_GBC module

nextflow.enable.dsl = 2

//

process EXTRACT_GBC {

  label 'scLT'
  tag "${sample_name}"

  input:
  tuple val(sample_name), val(in_folder)
  path search_patterns

  output:
  tuple val(sample_name), path('GBC_not_corrected.tsv'), emit: GBC

  script:
  """
  zcat ${in_folder}/*_R1_*.fastq.gz ${in_folder}/*_R2_*.fastq.gz | \
  awk 'NR % 4 == 2' | \
  egrep -f ${search_patterns} -o | \
  awk '{print substr(\$0, 23, 18);}' \
  > GBC_not_corrected.tsv
  """

  stub:
  """
  touch GBC_not_corrected.tsv
  """
  
}
