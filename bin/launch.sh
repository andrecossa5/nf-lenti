#!/bin/sh

source ~/.bashrc

path_wd=/data/cossa2_dare/AML_MITO
path_pipeline=/data/cossa2/phylo_inference/
path_data=$path_wd/data/samples_no_NUMTs/
path_results=$path_wd/results/trial_simpler_and_final/

# For separate lineages
for lineage in barcodes lymphoid malignant; do
   
    nextflow run $path_pipeline/main.nf \
    -c $path_wd/nextflow.config \
    --cell_file $lineage.txt \
    --outdir $path_results/$lineage \
    -params-file $path_wd/params.json \
    -profile conda_garr \
    -entry phylo 

done