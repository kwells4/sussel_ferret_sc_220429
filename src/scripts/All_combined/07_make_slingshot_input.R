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


write_data <- function(merged_seurat, input_dir, sample_group,
                     cluster_name){
  
  subset_seurat <- subset(merged_seurat, subset = genotype == sample_group)
  
  # Pull out same number of PCs as used to make the UMAP
  dim_red_all <- Embeddings(subset_seurat, reduction = "pca")[ , 1:30]
  
  # Pull out cluster information
  clusters <- subset_seurat[[cluster_name]]
  
  specific_seurat <- here("results", sample, "R_analysis")
  
  sample_seurat <- readRDS(here("results", input_dir, "R_analysis",
                                "rda_obj", "seurat_processed.rds"))
  
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

write_data(merged_seurat, "WT_combined", "WT", "wt_cluster")
write_data(merged_seurat, "CFKO_combined", "CFKO", "cfko_cluster")
