# Document information
# This document makes all of the figures that are seen in the manuscript.

library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(viridis)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

counts <- GetAssayData(seurat_data, assay = "RNA", slot = "data")
dimensionality_reduction <- Embeddings(seurat_data, reduction = "rna.umap")

meta_data <- seurat_data[[]]

meta_data$cluster <- paste(meta_data$genotype, meta_data$ind_cluster,
                           sep = "_")

keep_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
                  "S.Score", "G2M.Score", "Phase", "sample", "day", "genotype",
                  "cluster", "celltype")

save_meta <- meta_data[ , keep_columns]

lineage_columns <- colnames(meta_data)[grepl("Lineage", colnames(meta_data))]

save_lineage <- meta_data[, lineage_columns]

write_dir <- file.path(save_dir, "files", "publish_files")
ifelse(!dir.exists(write_dir), dir.create(write_dir), FALSE)

write.csv(counts, 
          file.path(write_dir, "normalized_counts.csv"))


write.csv(save_meta, 
          file.path(write_dir, "merged_metadata.csv"))


write.csv(save_lineage, 
          file.path(write_dir, "lineage_scores.csv"))
