// FIX_TAGS module

nextflow.enable.dsl = 2

//

process FIX_TAGS {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("mitobam_fixed_tags.bam"), emit: bam

  script:
  """
  python ${baseDir}/bin/maester/fix_tags.py ${bam}
  """

  stub:
  """
  touch mitobam_fixed_tags.bam
  """

}