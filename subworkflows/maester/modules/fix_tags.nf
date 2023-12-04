// FIX_TAGS module

nextflow.enable.dsl = 2

//

process FIX_TAGS {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(mitobam_no_UB_CB)

  output:
  tuple val(sample_name), path("mitobam_fixed.bam"), emit: mitobam

  script:
  """
  python ${baseDir}/bin/maester/fix_tags.py ${mitobam_no_UB_CB}
  """

  stub:
  """
  touch mitobam_fixed.bam
  """

}