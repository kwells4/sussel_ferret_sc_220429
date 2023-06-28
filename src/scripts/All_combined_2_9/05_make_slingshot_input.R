# Document information
# This document pulls out cluster and embedding information to pass to
# slingshot. Both joint embeddings from the full dataset and embeddings from
# each genotype alone are pulled. Clusters are pulled from each genotype
# alone.

# Note cluster resolution of 1.2 was pretty good for both
# Note cluster resolution of 0.8 was good for WT bad for CFKO
# Note cluster resolution of 1.5 was similar to 1.2
# The CFKO problem is I can't get all cells that were in the original Lineage
# 1 together

library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(slingshot)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined_2_9"

sample_dir <- here("results", sample, "R_analysis")

write_data <- function(merged_seurat, sample_seurat, sample_group,
                     cluster_name, sample_dir){
  
  subset_seurat <- subset(merged_seurat, subset = genotype == sample_group)
  
  # Pull out same number of PCs as used to make the UMAP
  dim_red_all <- Embeddings(subset_seurat, reduction = "pca")[ , 1:30]
  
  # Pull out cluster information
  clusters <- sample_seurat[[cluster_name]]

  dim_red_sample <- Embeddings(sample_seurat, reduction = "pca")[ , 1:30]
  
  barcode_mapping <- lapply(c(sample_seurat, subset_seurat), function(x){
    barcode_mapping_one <- x[[]] %>%
      dplyr::mutate(barcode = gsub("_[0-9]*", "", rownames(.))) %>%
      dplyr::select(sample, barcode) %>%
      dplyr::mutate(sample_barcode = paste(sample, barcode, sep  = "_"))
    return(barcode_mapping_one)
  })

  barcode_mapping[[2]] <- barcode_mapping[[2]][
    order(match(barcode_mapping[[2]]$sample_barcode,
                barcode_mapping[[1]]$sample_barcode)),]
  
  if(!identical(barcode_mapping[[1]]$sample_barcode,
                barcode_mapping[[2]]$sample_barcode)){
    stop("barcodes don't match")
  } else if (!identical(rownames(barcode_mapping[[1]]),
                        rownames(dim_red_sample))){
    stop("barcodes don't match")
  }
  
  rownames(dim_red_sample) <- rownames(barcode_mapping[[2]])
  
  # Save to a new directory
  save_dir <- file.path(sample_dir, "files", "slingshot")
  
  ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)
  
  write.table(dim_red_all, sep = "\t", quote = FALSE, row.names = TRUE,
              file = file.path(save_dir, paste0(sample_group, "_pca_all.tsv")))
  
  write.table(clusters, sep = "\t", quote = FALSE, row.names = TRUE,
              file = file.path(save_dir, paste0(sample_group, "_clusters.tsv")))
  
  write.table(dim_red_sample, sep = "\t", quote = FALSE, row.names = TRUE,
              file = file.path(save_dir, paste0(sample_group, "_pca.tsv")))
  
}

# Read in seurat objects
merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


seurat_wt <- readRDS(file.path(sample_dir, "rda_obj/seurat_wt.rds"))
seurat_cfko <- readRDS(file.path(sample_dir, "rda_obj/seurat_cfko.rds"))

# Get clusters for slingshot, lower resolution than the cell types
seurat_cfko <- FindClusters(seurat_cfko, resolution = 0.7)
seurat_cfko$cfko_slingshot_cluster <- seurat_cfko$RNA_snn_res.0.7

all_plots <- lapply(unique(seurat_cfko$cfko_slingshot_cluster), function(x){
  plotDimRed(seurat_cfko, col_by = "cfko_slingshot_cluster",
             plot_type = "rna.umap",
             highlight_group = TRUE, group = x, 
             meta_data_col = "cfko_slingshot_cluster")
})

names(all_plots) <- unique(seurat_cfko$cfko_slingshot_cluster)

seurat_cfko_sce <- as.SingleCellExperiment(seurat_cfko)

# Get a sense of parameters
test_lineages <- getLineages(seurat_cfko_sce, 
                             clusterLabels = "cfko_slingshot_cluster",
                             reducedDim = "PCA",
                             start.clus = 0, use.median = TRUE,
                             dist.method = "mnn")
cfko_results <- SlingshotDataSet(test_lineages)

# Was 1.5
seurat_wt <- FindClusters(seurat_wt, resolution = 1.5)
seurat_wt$wt_slingshot_cluster <- seurat_wt$RNA_snn_res.1.5

plotDimRed(seurat_wt, col_by = "wt_slingshot_cluster",
           plot_type = "rna.umap")

wt_cluster_info <- seurat_wt[[]] %>%
  dplyr::select(wt_slingshot_cluster, RNA_cluster) %>%
  dplyr::rename(wt_cluster = RNA_cluster)

cfko_cluster_info <- seurat_cfko[[]] %>%
  dplyr::select(cfko_slingshot_cluster, RNA_cluster) %>%
  dplyr::rename(cfko_cluster = RNA_cluster)


merged_seurat <- AddMetaData(merged_seurat, metadata = wt_cluster_info)
merged_seurat <- AddMetaData(merged_seurat, metadata = cfko_cluster_info)


all_plots <- lapply(unique(merged_seurat$wt_slingshot_cluster), function(x){
  plotDimRed(merged_seurat, col_by = "wt_slingshot_cluster",
             plot_type = "rna.umap",
             highlight_group = TRUE, group = x, 
             meta_data_col = "wt_slingshot_cluster")
})


names(all_plots) <- unique(merged_seurat$wt_slingshot_cluster)

seurat_wt_sce <- as.SingleCellExperiment(seurat_wt)

# Get a sense of parameters
test_lineages <- getLineages(seurat_wt_sce, 
                             clusterLabels = "wt_slingshot_cluster",
                             reducedDim = "PCA",
                             start.clus = 4, use.median = TRUE,
                             dist.method = "mnn")

wt_results <- SlingshotDataSet(test_lineages)
pheatmap::pheatmap(wt_results@adjacency)

write_data(merged_seurat, seurat_wt, "WT", "wt_slingshot_cluster",
           sample_dir)

write_data(merged_seurat, seurat_cfko, "CFKO", "cfko_slingshot_cluster",
           sample_dir)

saveRDS(seurat_cfko, file.path(sample_dir, "rda_obj/seurat_cfko.rds"))
saveRDS(seurat_wt, file.path(sample_dir, "rda_obj/seurat_wt.rds"))