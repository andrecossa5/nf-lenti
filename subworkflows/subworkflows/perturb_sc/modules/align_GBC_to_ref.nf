// GBC_TO_FASTA module

nextflow.enable.dsl = 2

//

process ALIGN_GBC {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(custom_ref_genome), path(GBC_to_align_fa)

  output:
  tuple val(sample_name), path("reads_aligned_to_in_vitro_ref.tsv"), emit: names

  script:
  """
  bowtie2 \
  -N 0 \
  -L 9 \
  -i L,0,0.35 \
  -x ${custom_ref_genome}/GBC_reference \
  --end-to-end \
  --score-min L,0,0 \
  --n-ceil L,0,0 \
  -f ${GBC_to_align_fa} \
  -t \
  --no-1mm-upfront \
  --no-unal \
  --no-hd \
  --no-sq \
  --al-gz GBC_aligned.fa.gz \
  > GBC_aligned.sam

  rm GBC_aligned.sam

  zcat GBC_aligned.fa.gz \
  | awk '\$1 ~ /^>/' \
  | sed 's/>/@/' \
  > reads_aligned_to_in_vitro_ref.tsv

  rm GBC_aligned.fa.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch reads_aligned_to_in_vitro_ref.tsv
  """

}