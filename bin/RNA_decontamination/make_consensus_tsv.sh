#!/bin/sh

source ~/.bashrc


type_ex=$1
sample=$2
cell=$3
path=/hpcnfs/scratch/PGP/cossa_campani/results/$type_ex/$sample/
fasta_path=/hpcnfs/data/PGP/reference_genomes/custom/make_new_index/cassette_up.fa
min_reads=3

cp $fasta_path $path
fasta_path="${path}cassette_up.fa"
bwa index $fasta_path
samtools faidx $fasta_path
samtools dict -o /hpcnfs/scratch/PGP/cossa_campani/results/sc_gbc/MDA_clones/cassette_up.dict $fasta_path
cell_path="${path}cbc/${cell}/"

cd $cell_path

R1="reads_${cell}_R1.fastq"
R2="reads_${cell}_R2.fastq"
fgbio --compression 1 FastqToBam -i $R1 $R2 -r 16B12M +T -o unmapped.bam --sample MDA_clones --library UMI
samtools fastq unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y $fasta_path - | fgbio -Xmx4g --compression 1 --async-io ZipperBams --unmapped unmapped.bam --ref $fasta_path --output mapped.bam
fgbio GroupReadsByUmi -s Identity -e 0 -i mapped.bam -o grouped.bam --raw-tag RX -T MI
fgbio CallMolecularConsensusReads -t RX -i grouped.bam -o consensus.bam -M $min_reads
fgbio FilterConsensusReads -i consensus.bam -o consensus_filtered.bam -r $fasta_path -M $min_reads -e 0.1 -N 20 -E 0.025

samtools index consensus_filtered.bam

cd /hpcnfs/home/ieo6943/mito_preprocessing/bin/RNA_decontamination/
mamba activate Mi_To_G
python from_consensus_bam_to_tsv.py --input_path ${cell_path}consensus_filtered.bam --output_path ${cell_path}consensus_filtered.tsv --cbc $cell