# Document information
# This document makes all of the figures that are seen in the manuscript.

library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(viridis)
library(clustree)
library(harmony)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "All_combined_2_9"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

all_colors2 <- c("acinar" = "#D4405B",
                 "ductal" = "#A5903E",
                 "Prolif_acinar" = "#55A470",
                 "Prolif_ductal" = "#767FC9",
                 "centroacinar_progenitor" = "#297878",
                 "centroacinar" = "#78295D")

# Batch correction -------------------------------------------------------------
seurat_data <- RunHarmony(seurat_data, c("day"),
                          plot_convergence = TRUE,
                          theta = 3)

harmony_plots <- plotDimRed(seurat_data, col_by = "sample",
                            plot_type = "harmony")

# UMAP -------------------------------------------------------------------------

RNA_pcs <- 30
ADT_pcs <- 8

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}


set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = "RNA", HTO = FALSE,
                         reduction = "harmony")

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data <- FindClusters(seurat_data, resolution = c(0.5, 0.8, 1, 1.2,
                                                        1.5, 2.0, 3.0))
clustree(seurat_data)

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 3.0, assay = "RNA", HTO = FALSE,
                         reduction = "harmony")


seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data$full_corrected_cluster <- seurat_data$RNA_cluster

cm <- confusionMatrix(seurat_data$full_corrected_cluster, seurat_data$celltype )
cm <- cm / rowSums(cm)

pdf(file.path(save_dir, "images", "revision",
              "full_batch_corrected_cluster_vs_celltype.pdf"))
print(pheatmap::pheatmap(cm))

dev.off()

pdf(file.path(save_dir, "images", "revision", 
              "celltype_on_corrected_umap.pdf"))
print(plotDimRed(seurat_data, col_by = "celltype", plot_type = "harmony.umap",
           color = all_colors2))

print(plotDimRed(seurat_data, col_by = "sample", plot_type = "harmony.umap"))

dev.off()

pdf(file.path(save_dir, "images", "revision", "slingshot_lineages.pdf"))
# Called lineage 1 in paper
print(plotDimRed(seurat_data, col_by = "wt_Lineage1", plot_type = "harmony.umap")[[1]] +
  ggplot2::ggtitle("wt_Lineage1"))

# Called lineage 2 in paper
print(plotDimRed(seurat_data, col_by = "wt_Lineage9", plot_type = "harmony.umap")[[1]] +
        ggplot2::ggtitle("wt_Lineage2"))

# Called lineage 3 in paper
print(plotDimRed(seurat_data, col_by = "wt_Lineage7", plot_type = "harmony.umap")[[1]] +
        ggplot2::ggtitle("wt_Lineage3"))

# Called lineage 1 in paper
print(plotDimRed(seurat_data, col_by = "cfko_Lineage1", plot_type = "harmony.umap")[[1]] +
        ggplot2::ggtitle("cfko_Lineage1"))

# Called lineage 2 in paper
print(plotDimRed(seurat_data, col_by = "cfko_Lineage4", plot_type = "harmony.umap")[[1]] +
        ggplot2::ggtitle("cfko_Lineage2"))

# Called lineage 3 in paper
print(plotDimRed(seurat_data, col_by = "cfko_Lineage2", plot_type = "harmony.umap")[[1]] +
        ggplot2::ggtitle("cfko_Lineage3"))

dev.off()

# Slingshot --------------------------------------------------------------------
# Look for starting cluster
plotDimRed(seurat_data, col_by = "genotype_cluster_celltype", plot_type = "rna.umap")

plotDimRed(seurat_data, col_by = "genotype_cluster_celltype", plot_type = "rna.umap",
           highlight_group = TRUE, group = "WT_4_ductal",
           meta_data_col = "genotype_cluster_celltype")


plotDimRed(seurat_data, col_by = "genotype_cluster_celltype", plot_type = "rna.umap",
           highlight_group = TRUE, group = "CFKO_2_ductal",
           meta_data_col = "genotype_cluster_celltype")

cm <- confusionMatrix(seurat_data$genotype_cluster_celltype,
                      seurat_data$full_corrected_cluster)

cm <- cm / rowSums(cm)
pheatmap::pheatmap(cm)

# Start cluster = cluster 7 for WT 1 for CFKO
plotDimRed(seurat_data, col_by = "full_corrected_cluster", plot_type = "harmony.umap",
           highlight_group = TRUE, group = "7",
           meta_data_col = "full_corrected_cluster")
plotDimRed(seurat_data, col_by = "full_corrected_cluster", plot_type = "harmony.umap",
           highlight_group = TRUE, group = "1",
           meta_data_col = "full_corrected_cluster")


# Make input for slingshot
# PCA 1:30


make_slingshot_info <- function(seurat_data, sample_group, save_dir){
  # Subset to just the genotype of interest
  subset_seruat <- subset(seurat_data, subset = genotype == sample_group)

  # Want only cells that are in the genotype to be part of the clustering
  # UMAP of gene expression
  set.seed(0)
  umap_data <- group_cells(subset_seruat, sample, save_dir = NULL,
                           nPCs = RNA_pcs,
                           resolution = 3.0, assay = "RNA", HTO = FALSE,
                           reduction = "harmony")
  
  subset_seruat <- umap_data$object
  
  dim_red <- Embeddings(subset_seruat, reduction = "harmony")[ , 1:30]
  
  clusters <- subset_seruat$RNA_cluster
  

  write.table(clusters, sep = "\t", quote = FALSE, row.names = TRUE,
              file = file.path(save_dir, paste0(sample_group, "_clusters.tsv")))
  
  write.table(dim_red, sep = "\t", quote = FALSE, row.names = TRUE,
              file = file.path(save_dir, paste0(sample_group, "_pca.tsv")))
  
  return(subset_seruat)
  
}

wt <- make_slingshot_info(seurat_data, "WT", 
                    save_dir = file.path(save_dir, "images", "revision"))

cfko <- make_slingshot_info(seurat_data, "CFKO", 
                    save_dir = file.path(save_dir, "images", "revision"))


# Find starting cluster
plotDimRed(wt, "wt_Lineage1", plot_type = "harmony.umap")
plotDimRed(wt, "RNA_cluster", plot_type = "harmony.umap")
plotDimRed(wt, "RNA_cluster", plot_type = "harmony.umap",
           highlight_group = TRUE, group = "26", meta_data_col = "RNA_cluster")

plotDimRed(wt, "RNA_cluster", plot_type = "harmony.umap",
           highlight_group = TRUE, group = "13", meta_data_col = "RNA_cluster")

featDistPlot(wt, geneset = "CRABP2", combine = FALSE, sep_by = "RNA_cluster",
             col_by = "RNA_cluster")
featDistPlot(wt, geneset = c("wt_Lineage1", "wt_Lineage2"),
             combine = FALSE, sep_by = "RNA_cluster",
             col_by = "RNA_cluster")


plotDimRed(cfko, "cfko_Lineage3", plot_type = "harmony.umap")
featDistPlot(cfko, geneset = c("cfko_Lineage1", "cfko_Lineage2"),
             combine = FALSE, sep_by = "RNA_cluster",
             col_by = "RNA_cluster")
plotDimRed(cfko, "RNA_cluster", plot_type = "harmony.umap")
plotDimRed(cfko, "RNA_cluster", plot_type = "harmony.umap",
           highlight_group = TRUE, group = "18", meta_data_col = "RNA_cluster")

featDistPlot(cfko, geneset = "CRABP2", combine = FALSE, sep_by = "RNA_cluster",
             col_by = "RNA_cluster")
