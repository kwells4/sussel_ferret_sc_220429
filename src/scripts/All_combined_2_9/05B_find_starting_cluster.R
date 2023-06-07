# Document information
# This document uses visualizations to find a likely starting cluster based on 
# the clusters that are enriched in the D2 sample.

library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined_2_9"

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

seurat_wt <- readRDS(file.path(sample_dir, "rda_obj/seurat_wt.rds"))
seurat_cfko <- readRDS(file.path(sample_dir, "rda_obj/seurat_cfko.rds"))

# Add cluster info from individual analysis
wt_cluster_info <- seurat_wt[[]] %>%
  dplyr::select(RNA_cluster) %>%
  dplyr::rename(wt_cluster = RNA_cluster)

cfko_cluster_info <- seurat_cfko[[]] %>%
  dplyr::select(RNA_cluster) %>%
  dplyr::rename(cfko_cluster = RNA_cluster)


merged_seurat <- AddMetaData(merged_seurat, metadata = wt_cluster_info)
merged_seurat <- AddMetaData(merged_seurat, metadata = cfko_cluster_info)

# Find clusters that are enriched in the D2 sample
cM <- confusionMatrix(merged_seurat$cfko_cluster, merged_seurat$sample)
cM <- cM / rowSums(cM)
pheatmap::pheatmap(cM)


cM <- confusionMatrix(merged_seurat$wt_cluster, merged_seurat$sample)
cM <- cM / rowSums(cM)
pheatmap::pheatmap(cM)

# Plot UMAPs to identify the cluster that is furthest from the D5-9 samples
plotDimRed(merged_seurat, "wt_cluster", plot_type = "rna.umap")

plotDimRed(merged_seurat, "wt_cluster", plot_type = "rna.umap",
           highlight_group = TRUE, meta_data_col = "wt_cluster",
           group = 4)

plotDimRed(merged_seurat, "cfko_cluster", plot_type = "rna.umap")

plotDimRed(merged_seurat, "cfko_cluster", plot_type = "rna.umap",
           highlight_group = TRUE, meta_data_col = "cfko_cluster",
           group = 2)

plotDimRed(merged_seurat, "cfko_cluster", plot_type = "rna.umap",
           highlight_group = TRUE, meta_data_col = "cfko_cluster",
           group = 9)


