library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "All_combined"

wt_samples <- c("WT_D2", "WT_D5", "WT_D7", "WT_D9", "WT_D14")

cfko_samples <- c("CFKO_D2", "CFKO_D5", "CFKO_D7", "CFKO_D9", "CFKO_D14")


sample_dir <- here("results", sample, "R_analysis")

merged_seurat <- readRDS(file.path(sample_dir, "rda_obj",
                                   "seurat_processed.rds"))


wt_seurat <- subset(merged_seurat, subset = sample %in% wt_samples)

cfko_seurat <- subset(merged_seurat, subset = sample %in% cfko_samples)

# Pull out same number of PCs as used to make the UMAP
wt_dim_red <- Embeddings(wt_seurat, reduction = "pca")[ , 1:32]

cfko_dim_red <- Embeddings(cfko_seurat, reduction = "pca")[, 1:32]


# Pull out cluster information
wt_clusters <- wt_seurat[["sample_cluster"]]

cfko_clusters <- cfko_seurat[["sample_cluster"]]

# Save to a new directory
save_dir <- file.path(sample_dir, "files", "slingshot")

ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

write.table(wt_dim_red, sep = "\t", quote = FALSE, row.names = TRUE,
            file = file.path(save_dir, "WT_pca.tsv"))

write.table(cfko_dim_red, sep = "\t", quote = FALSE, row.names = TRUE,
            file = file.path(save_dir, "CFKO_pca.tsv"))

write.table(wt_clusters, sep = "\t", quote = FALSE, row.names = TRUE,
            file = file.path(save_dir, "WT_clusters.tsv"))

write.table(cfko_clusters, sep = "\t", quote = FALSE, row.names = TRUE,
            file = file.path(save_dir, "CFKO_clusters.tsv"))