#! /usr/bin/env bash

#BSUB -J slinghot
#BSUB -o logs/wt_slingshot_%J.out
#BSUB -e logs/wt_slingshot_%J.err
#BSUB -R "select[mem>200] rusage[mem=200]" 
#BSUB -q rna

script=src/scripts/All_combined/08_run_slingshot.R

data_path=results/All_combined/R_analysis/files/slingshot

pca_file=$data_path/WT_pca.tsv 

cluster_file=$data_path/WT_clusters.tsv

start_cluster="3"

save_name=$data_path/WT_second_cluster

seurat_object=results/All_combined/R_analysis/rda_obj/seurat_processed.rds

nknots="7"

set -o nounset -o pipefail -o errexit -x

Rscript --vanilla \
    $script \
    $pca_file \
    $cluster_file \
    $start_cluster \
    $seurat_object \
    $save_name \
    $nknots
