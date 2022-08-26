library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "CFKO_D5"

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

cell_cycle_genes <- readRDS(here("files/cell_cycle_genes.rds"))

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

seurat_data <- CellCycleScoring(seurat_data,
                                s.features = cell_cycle_genes$s_genes$gene_id,
                                g2m.features = cell_cycle_genes$g2m_genes$gene_id)

plotDimRed(seurat_data, "Phase", plot_type = "rna.umap")

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
