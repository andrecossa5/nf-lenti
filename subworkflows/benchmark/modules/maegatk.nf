// MAEGATK module

nextflow.enable.dsl = 2

//

process MAEGATK {

  tag "${sample_name}: ${cell}"

  input:
  tuple val(sample_name), val(cell), path(bam)
  path(reference)

  output:
  tuple val(sample_name), 
    val(cell), 
    path("${cell}.A.txt"),
    path("${cell}.C.txt"),  
    path("${cell}.G.txt"),  
    path("${cell}.T.txt"),
    path("${cell}.coverage.txt"),  
    path("${cell}.depth.txt"), emit: tables
  
  script:
  """
  python ${baseDir}/bin/benchmark/oneSample_maegatk.py ${bam} ${cell} chrM.fa 
  rm *bam *bai *fastq *sam *log
  """

  stub:
  """
  touch ${cell}.A.txt  
  touch ${cell}.C.txt  
  touch ${cell}.G.txt  
  touch ${cell}.T.txt 
  touch ${cell}.coverage.txt  
  touch ${cell}.depth.txt
  """

}


//

process COLLAPSE_MAEGATK {

  tag "${sample_name}"
  publishDir "${params.sc_outdir}/${sample_name}/", mode: 'copy'

  input:
  tuple val(sample_name), 
      val(cells), 
      path(A), 
      path(C), 
      path(G), 
      path(T), 
      path(coverage),
      path(depth)

  output:
  tuple val(sample_name), path(tables), emit: ch_output
 
  script:
  """ 
  cat *.A.txt | gzip --fast > A.txt.gz
  cat *.C.txt | gzip --fast > C.txt.gz
  cat *.T.txt | gzip --fast > T.txt.gz
  cat *.G.txt | gzip --fast > G.txt.gz
  cat *.coverage.txt | gzip --fast > coverage.txt.gz
  cat *.depth.txt | gzip --fast > depth.txt.gz
  mkdir tables 
  mv *.gz tables
  """

  stub:
  """
  mkdir tables
  """

}
