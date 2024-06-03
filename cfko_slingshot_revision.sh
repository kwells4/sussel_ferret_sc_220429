#! /usr/bin/env bash

#BSUB -J slinghot
#BSUB -o logs/cfko_slingshot_%J.out
#BSUB -e logs/cfko_slingshot_%J.err
#BSUB -R "select[mem>200] rusage[mem=200]" 
#BSUB -q rna

# Run in slingshot conda environment

script=src/scripts/All_combined_2_9/06_run_slingshot.R

data_path=results/All_combined_2_9/R_analysis/images/revision

pca_file=$data_path/CFKO_pca.tsv 

cluster_file=$data_path/CFKO_clusters.tsv

start_cluster="1"

save_name=$data_path/CFKO

seurat_object=results/All_combined_2_9/R_analysis/rda_obj/seurat_processed.rds

nknots="5"

set -o nounset -o pipefail -o errexit -x

Rscript --vanilla \
    $script \
    $pca_file \
    $cluster_file \
    $start_cluster \
    $seurat_object \
    $save_name \
    $nknots