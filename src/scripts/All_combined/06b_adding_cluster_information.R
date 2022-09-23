library(slingshot)
library(scAnalysisR)
library(Seurat)
library(here)
library(tidyverse)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

all_samples <- "All_combined"

all_sample_dir <- here("results", all_samples, "R_analysis")

merged_seurat <- readRDS(file.path(all_sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


add_cluster_data <- function(sample, cluster_name, subset_name,
                             merged_seurat){
  sample_dir <- here("results", sample, "R_analysis")
  
  sample_seurat <- cfko_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                                    "seurat_processed.rds"))
  
  
  # Make sure barcodes match
  cluster_info <- sample_seurat[[]] %>%
    dplyr::select(sample, uncorrected_cluster) %>%
    dplyr::mutate(barcode = gsub("_[0-9]+", "", rownames(.))) %>%
    dplyr::mutate(barcode_sample = paste(barcode, sample, sep = "_")) 
  
  merged_seurat_info <- merged_seurat[[]] %>%
    dplyr::mutate(barcode = gsub("_[0-9]+", "", rownames(.))) %>%
    dplyr::filter(genotype == subset_name) %>%
    dplyr::select(barcode, sample) %>%
    dplyr::mutate(barcode_sample = paste(barcode, sample, sep = "_"))
  
  if(nrow(cluster_info) != nrow(merged_seurat_info)){
    stop("different number of rows for your two objects")
  }
  
  cluster_info <- cluster_info[order(
    match(cluster_info$barcode_sample,
          merged_seurat_info$barcode_sample)),]
  
  if(!identical(cluster_info$barcode_sample,
                merged_seurat_info$barcode_sample)){
    stop("Something went wrong when reordering!")
  }
  
  # Change rownames and select only the cluster column
  rownames(cluster_info) <- rownames(merged_seurat_info)
  
  cluster_info <- cluster_info %>%
    dplyr::select(uncorrected_cluster)
  
  colnames(cluster_info) <- cluster_name
  cluster_info[[paste0(subset_name, "_barcodes")]] <- rownames(cluster_info)
  
  merged_seurat <- AddMetaData(merged_seurat, metadata = cluster_info)
  
  return(merged_seurat)
}

merged_seurat <- add_cluster_data(sample = "CFKO_combined",
                                  cluster_name = "cfko_cluster",
                                  subset_name = "CFKO",
                                  merged_seurat = merged_seurat)

merged_seurat <- add_cluster_data(sample = "WT_combined",
                                  cluster_name = "wt_cluster",
                                  subset_name = "WT",
                                  merged_seurat = merged_seurat)

plotDimRed(merged_seurat, "cfko_cluster", plot_type = "rna.umap")
plotDimRed(merged_seurat, "wt_cluster", plot_type = "rna.umap")


# Okay, so the clusters here look good...

# Still to do -->
# Repeat analysis with these new clusters --> Try with both the full PCA
# And the PCA for the genotypes individually... should be similar.

# cluster_simalarity
cm <- confusionMatrix(merged_seurat$uncorrected_cluster,
                      merged_seurat$sample_cluster)

pheatmap::pheatmap(log(cm + 1))

# Think about how to show percent - ie what percent of cells in the
# original cluster fall into the new cluster
# do this with dplyr --> 
seurat_data <- merged_seurat[[]] %>%
  dplyr::select(uncorrected_cluster, sample_cluster) %>%
  dplyr::group_by(uncorrected_cluster) %>%
  dplyr::add_count(name = "uncorrected_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_cluster) %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(all_data = paste(uncorrected_cluster, sample_cluster,
                                 sep = "_")) %>%
  dplyr::group_by(all_data) %>%
  dplyr::add_count(name = "all_count") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_sample = all_count / sample_count,
                percent_all = all_count / uncorrected_count)

max_vals <- seurat_data %>%
  dplyr::group_by(sample_cluster) %>%
  dplyr::summarise(max = max(percent_sample))


# All except 24 clusters have 75% plus cells mapping to one of the large
# clusters

seurat_data <- merged_seurat[[]] %>%
  dplyr::select(uncorrected_cluster, RNA_combined_celltype) %>%
  dplyr::group_by(uncorrected_cluster) %>%
  dplyr::add_count(name = "uncorrected_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(RNA_combined_celltype) %>%
  dplyr::add_count(name = "celltype_count") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(all_data = paste(uncorrected_cluster, RNA_combined_celltype,
                                 sep = "_")) %>%
  dplyr::group_by(all_data) %>%
  dplyr::add_count(name = "all_count") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_celltype = all_count / celltype_count,
                percent_all = all_count / uncorrected_count)

max_vals <- seurat_data %>%
  dplyr::group_by(uncorrected_cluster) %>%
  dplyr::summarise(max = max(percent_all))


# This is less good... summary - basically, the larger clusters are neatly
# subset into smaller, higher resolution clusters with the individual datasets
# being analyzed alone. Sometimes, this resolution allowed us to call cell
# types that we would not have otherwise called.

saveRDS(merged_seurat, file.path(all_sample_dir,
                                 "rda_obj/seurat_processed.rds"))
