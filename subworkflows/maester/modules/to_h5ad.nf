// TO_H5AD module

nextflow.enable.dsl = 2

//

process TO_H5AD {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(output)

  output:
  tuple val(sample_name), path("AFM.h5ad"), emit: afm

  script:
  """
  python ${baseDir}/bin/maester/to_h5ad.py ${output}
  """
  
  stub:
  """
  touch AFM.h5ad
  """

}