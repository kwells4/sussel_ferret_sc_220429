library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

source("src/scripts/functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "CFKO_D7"

cell_types <- "RNA_combined_celltype"
clusters <- "RNA_cluster"
celltype_two <- "RNA_tabula_muris_celltype"

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

rownames(seurat_data)[grepl("SPINK", rownames(seurat_data))]

plotDimRed(seurat_data, col_by = cell_types, plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = celltype_two, plot_type = "rna.umap")

# Ductal markers ---------------------------------------------------------------
gene <- "KRT19"
plotDimRed(seurat_data, col_by = gene, plot_type = "rna.umap")
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype",
                  col_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype_tm",
                  col_by = celltype_two))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = celltype_two))

gene <- "SOX9"
plotDimRed(seurat_data, col_by = gene, plot_type = "rna.umap")
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype",
                  col_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype_tm",
                  col_by = celltype_two))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = celltype_two))

# Epithelial markers -----------------------------------------------------------
gene <- "CDH1"
plotDimRed(seurat_data, col_by = gene, plot_type = "rna.umap")
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype",
                  col_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype_tm",
                  col_by = celltype_two))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = celltype_two))

# Acinar markers ---------------------------------------------------------------
gene <- "SPINK1"
plotDimRed(seurat_data, col_by = gene, plot_type = "rna.umap")
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype",
                  col_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = "cluster_celltype_tm",
                  col_by = celltype_two))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = cell_types))
plot(featDistPlot(seurat_data, geneset = gene, sep_by = celltype_two))
