library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))

plotDimRed(merged_seurat, "wt_cluster", plot_type = "rna.umap")

plotDimRed(merged_seurat, "wt_cluster", plot_type = "rna.umap",
           highlight_group = TRUE, meta_data_col = "wt_cluster",
           group = 3)

plotDimRed(merged_seurat, "cfko_cluster", plot_type = "rna.umap")

plotDimRed(merged_seurat, "cfko_cluster", plot_type = "rna.umap",
           highlight_group = TRUE, meta_data_col = "cfko_cluster",
           group = 19)

plotDimRed(merged_seurat, "cfko_cluster", plot_type = "rna.umap",
           highlight_group = TRUE, meta_data_col = "cfko_cluster",
           group = 9)

# 
# "wt_cluster"
# "cfko_cluster"
# wt_seurat <- readRDS(file.path("results/WT_combined/R_analysis/rda_obj",
#                                "seurat_processed.rds"))
# 
# # Make sure barcodes match
# cluster_info <- wt_seurat[[]] %>%
#   dplyr::select(sample) %>%
#   dplyr::mutate(barcode = gsub("_[0-9]+", "", rownames(.))) %>%
#   dplyr::mutate(barcode_sample = paste(barcode, sample, sep = "_")) 
# 
# merged_seurat_info <- merged_seurat[[]] %>%
#   dplyr::mutate(barcode = gsub("_[0-9]+", "", rownames(.))) %>%
#   dplyr::filter(genotype == "WT") %>%
#   dplyr::select(barcode, sample, RNA_combined_celltype) %>%
#   dplyr::mutate(barcode_sample = paste(barcode, sample, sep = "_"))
# 
# merged_seurat_info <- merged_seurat_info[order(
#   match(merged_seurat_info$barcode_sample,
#         cluster_info$barcode_sample)),]
# 
# if(!identical(cluster_info$barcode_sample,
#               merged_seurat_info$barcode_sample)){
#   stop("Something went wrong when reordering!")
# }
# 
# wt_seurat$RNA_combined_celltype_combined <- merged_seurat_info$RNA_combined_celltype
# 
# # Test new clusters
# RNA_pcs <- 30
# 
# wt_seurat <- FindNeighbors(wt_seurat, dims = 1:15)
# wt_seurat <- FindClusters(wt_seurat, resolution = c(1, 1.2, 1.6, 2, 3, 4, 5, 6))
# 
# 
# confusion_matrix <- confusionMatrix(wt_seurat$RNA_snn_res.3,
#                                     wt_seurat$RNA_combined_celltype_combined)
# 
# confusion_matrix <- confusion_matrix[!is.na(rownames(confusion_matrix)),]
# 
# confusion_matrix <- confusion_matrix / Matrix::rowSums(confusion_matrix)
# 
# pheatmap::pheatmap(confusion_matrix)
# 
# 
# max_val <- apply(confusion_matrix, 1, max)
# 
# length(max_val)
# 
# length(max_val[max_val < 0.75])
# 
# 
# 
# confusion_matrix <- confusionMatrix(wt_seurat$RNA_combined_celltype,
#                                     wt_seurat$RNA_combined_celltype_combined)
# 
# confusion_matrix <- confusion_matrix[!is.na(rownames(confusion_matrix)),]
# 
# confusion_matrix <- confusion_matrix / Matrix::rowSums(confusion_matrix)
# 
# pheatmap::pheatmap(confusion_matrix)
# 
# 
# apply(confusion_matrix, 1, max)
