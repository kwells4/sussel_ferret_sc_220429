#! /usr/bin/env bash

#BSUB -J slinghot
#BSUB -o logs/cfko_slingshot_%J.out
#BSUB -e logs/cfko_slingshot_%J.err
#BSUB -R "select[mem>20] rusage[mem=20]" 
#BSUB -q rna

script=src/scripts/All_combined/08_run_slingshot.R

data_path=results/All_combined/R_analysis/files/slingshot

pca_file=$data_path/CFKO_pca.tsv 

cluster_file=$data_path/CFKO_clusters.tsv

start_cluster=CFKO_D2_0

save_name=$data_path/CFKO_res.rds

set -o nounset -o pipefail -o errexit -x

Rscript --vanilla \
    $script \
    $pca_file \
    $cluster_file \
    $start_cluster \
    $save_name