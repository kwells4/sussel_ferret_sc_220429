#! /usr/bin/env bash

#BSUB -J slinghot
#BSUB -o logs/wt_slingshot_%J.out
#BSUB -e logs/wt_slingshot_%J.err
#BSUB -R "select[mem>200] rusage[mem=200]" 
#BSUB -q rna

script=src/scripts/All_combined/08_run_slingshot.R

data_path=results/All_combined/R_analysis/files/slingshot

pca_file=$data_path/WT_pca.tsv 

pca_file_all=$data_path/WT_pca_all.tsv

cluster_file=$data_path/WT_clusters.tsv

start_cluster="4"

save_name=$data_path/WT

save_name_all=$data_path/WT_all

seurat_object=results/All_combined/R_analysis/rda_obj/seurat_processed.rds

nknots="7"

nknots_all="8"

set -o nounset -o pipefail -o errexit -x

Rscript --vanilla \
    $script \
    $pca_file \
    $cluster_file \
    $start_cluster \
    $seurat_object \
    $save_name \
    $nknots

Rscript --vanilla \
    $script \
    $pca_file_all \
    $cluster_file \
    $start_cluster \
    $seurat_object \
    $save_name_all \
    $nknots_all