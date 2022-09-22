#! /usr/bin/env bash

#BSUB -J slinghot
#BSUB -o logs/wt_slingshot_%J.out
#BSUB -e logs/wt_slingshot_%J.err
#BSUB -R "select[mem>20] rusage[mem=20]" 
#BSUB -q rna

script=src/scripts/All_combined/08_run_slingshot.R

data_path=results/All_combined/R_analysis/files/slingshot

pca_file=$data_path/WT_pca.tsv 

cluster_file=$data_path/WT_clusters.tsv

start_cluster=WT_D2_1

save_name=$data_path/WT_res.rds

set -o nounset -o pipefail -o errexit -x

Rscript --vanilla \
    $script \
    $pca_file \
    $cluster_file \
    $start_cluster \
    $save_name