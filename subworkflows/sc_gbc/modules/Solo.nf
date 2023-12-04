// SOLO module

nextflow.enable.dsl = 2

//

process SOLO {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(transcript_R1), path(transcript_R2)

  output:
  tuple val(sample_name), path('raw'), emit: raw 
  tuple val(sample_name), path('filtered'), emit: filtered 
  tuple val(sample_name), path('Features.stats'), emit: stats 
  tuple val(sample_name), path('Summary.csv'), emit: summary 
  tuple val(sample_name), path('Aligned.sortedByCoord.out.bam'), emit: bam 

  script:
  """
  STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${params.ref} \
    --readFilesIn ${transcript_R2} ${transcript_R1} \
    --readFilesCommand zcat \
    --outTmpDir tmp \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 50000000000 \
    --outSAMattributes NH HI nM AS CR UR CB UB \
    --soloType CB_UMI_Simple \
    --soloBarcodeReadLength 28 \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloCBwhitelist ${params.whitelist} \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter EmptyDrops_CR 
  mv Solo.out/Gene/* ./ &&
  rm -rf Solo.out &&
  gzip raw/*.tsv
  gzip raw/*.mtx
  gzip filtered/*.tsv
  gzip filtered/*.mtx
  """

  stub:
  """
  echo ${sample_name} > sample
  touch raw
  touch filtered
  touch Features.stats 
  touch Summary.csv
  touch Aligned.sortedByCoord.out.bam
  """

}