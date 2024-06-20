// TO_H5AD module

nextflow.enable.dsl = 2

//

process TO_H5AD {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(tables) 
  path(ref_txt)

  output:
  tuple val(sample_name), path("AFM.h5ad"), emit: afm

  script:
  """
  tr -d '\n' < ${ref_txt} > new.txt
  sed 's/./& /g' new.txt > space.txt
  python ${baseDir}/bin/maester/to_h5ad.py ${tables} space.txt
  """
  
  stub:
  """
  touch AFM.h5ad
  """

}